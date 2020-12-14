from pyGenicPipeline.utils import errors as ec
from pyGenicPipeline.utils import misc as mc
from pyGenicPipeline.core.Input import Input

import numpy as np
import time


class FilterSnps(Input):
    def __init__(self, args):
        super().__init__(args)
        # todo: Counts are not considering already filtered snps - change sum to be the sum of filter?
        self._filter_error_dict = {"Filter case": "Count", f"{self.freq} Discrepancy": 0, "MAF": 0,
                                   "Monomorphic": 0, "In Long Range LD": 0}

    def filter_snps(self, sm_dict, ref, chromosome, validation=None):
        """
        Large numbers of snps and individuals can lead to significant memory issues. This will filter the snps in chunks
        vus allowing it to run with less memory
        """
        t0 = self._assert_filter_snps()

        # Extract the information we need for filtering and then filter our references
        snp_list, chunks = self.chunked_snp_names(sm_dict, chunk_return=True)
        bp_positions = np.array_split(mc.variant_array(self.bp_position.lower(), sm_dict[self.sm_variants]), chunks)
        freqs = np.array_split(sm_dict[self.freq], chunks)

        accepted_snps = []
        for index, (snps, f, bp) in enumerate(zip(snp_list, freqs, bp_positions), start=1):
            print(f"Filtering chunk {index} out of {len(snp_list)}: {mc.terminal_time()}")
            # Setup the base filter from our extract chunks
            filter_dict = {self.snp_id: snps, self.freq: f, self.bp_position: bp,
                           self.filter_key: np.full(len(snps), True)}

            # Filter out undesirable snps, then store the new snps into a list
            filter_dict = self.filter_snp_chunk(ref, filter_dict, chromosome)

            # If we also have a validation file, filter on that as well. This will likely have a different std/freq so
            # may filter out additional snps
            if validation:
                filter_dict = self.filter_snp_chunk(validation, filter_dict, chromosome)

            # Store any remaining snps to a list
            accepted_snps.append(filter_dict[self.filter_key])

        # Return the filter summary dict
        t1 = mc.error_dict_to_terminal(self._filter_error_dict)

        print(f"Cleaned summary stats for Chromosome {chromosome} in {round(t1 - t0, 2)} Seconds\n")
        return mc.filter_array(sm_dict, mc.flatten(accepted_snps))

    def filter_snp_chunk(self, gen_file, filter_dict, chromosome):
        """
        From our cleaned summary statistics we can construct a set of information for our validation and reference
        genetic samples. These samples need there raw snps to be loaded, and then we can use these to calculate genetic
        frequencies and standard deviations.

        If we have information on summary frequencies, the user has specified a maf_min or wants to filter snps in long
        range LD then we can filter those out. This process will always filter out monomorphic snps.
        """
        # Load the raw snps from the .bed or .bgen file and use it to isolate statistical information
        raw_snps = self.isolate_raw_snps(gen_file, filter_dict[self.snp_id])
        filter_dict[self.f_std] = np.std(raw_snps, 1, dtype='float32')
        filter_dict[self.f_freq] = np.sum(raw_snps, 1, dtype='float32') / (2 * float(gen_file.iid_count))

        raw_snps = None
        print(f"Loaded raw snps {len(filter_dict[self.snp_id])} row removing from memory: raw_snps = {raw_snps}")

        # If the frequencies in the summary stats are not just a list of -1 errors then screen genetic snp frequencies
        if (self.freq_discrepancy < 1) and (np.sum(filter_dict[self.freq] == -1) != len(filter_dict[self.freq])):
            self._summary_frequencies(filter_dict)

        # Filter minor allele frequency SNPs.
        if self.maf_min > 0:
            self._maf_filter(filter_dict)

        # Filter any Monomorphic snps
        self._monomorphic_filter(filter_dict)

        # Filter long range LD if set
        if self.lr_ld_path:
            self._long_range_ld_filter(filter_dict, chromosome)

        return filter_dict

    def combined_filter_chunks(self):
        raise NotImplementedError("Not done yet!")

    def _summary_frequencies(self, filter_dict):
        """
        If the summary frequencies existed in the summary stats then cross check them with our genetic frequencies,
        removing anything that is outside of the frequency discrepancy allow by the user.
        """
        freq_filter = np.logical_or(
            np.absolute(filter_dict[self.freq] - filter_dict[self.f_freq]) < self.freq_discrepancy,
            np.absolute(filter_dict[self.freq] + (filter_dict[self.f_freq] - 1)) < self.freq_discrepancy)

        # Invalid frequencies from summary stats where coded as -1 so these should be removed
        freq_filter = np.logical_or(freq_filter, filter_dict[self.freq] <= 0)
        self._filter_error_dict[f"{self.freq} Discrepancy"] += len(freq_filter) - np.sum(freq_filter)

        # Log failures to the filter
        filter_dict[self.filter_key] = filter_dict[self.filter_key] * freq_filter

    def _maf_filter(self, filter_dict):
        """
        If the maf filter is greater than zero, then use this to remove any maf frequencies less that the frequency
        provided or greater than 1 - frequency provided.
        """
        maf_filter = (filter_dict[self.f_freq] > self.maf_min) * (filter_dict[self.f_freq] < (1 - self.maf_min))
        self._filter_error_dict["MAF"] += len(maf_filter) - np.sum(maf_filter)

        # Log filter to dict
        filter_dict[self.filter_key] = filter_dict[self.filter_key] * maf_filter

    def _monomorphic_filter(self, filter_dict):
        """Remove any Monomorphic snps, those with no variation"""
        monomorphic_filter = filter_dict[self.f_std] > 0
        self._filter_error_dict["Monomorphic"] += len(monomorphic_filter) - np.sum(monomorphic_filter)

        # Log filter to dict
        filter_dict[self.filter_key] = filter_dict[self.filter_key] * monomorphic_filter

    def _long_range_ld_filter(self, filter_dict, chromosome):
        """
        If we have filtering of long range LD in place then we load a dict of long range ld from Price et all 2008 as
        was done within, and data taken from, LdPred. Then we isolate the long range ld regions for this chromosome,
        and filter any in long range LD.
        """
        lr_pos = self.load_lr_ld_dict()[chromosome]
        if len(lr_pos) != 0:
            for key in lr_pos.keys():
                position = filter_dict[self.bp_position]
                long_filter = np.where((lr_pos[key]["start_pos"] < position) & (position < lr_pos[key]["end_pos"]),
                                       False, True)
                self._filter_error_dict["In Long Range LD"] += len(long_filter) - np.sum(long_filter)
                filter_dict[self.filter_key] = filter_dict[self.filter_key] * long_filter

    def _assert_filter_snps(self):
        """Different files require different load type operations, so load type must be set"""
        assert self.gen_type, ec.missing_arg(self.operation, "Load_Type")
        return time.time()
