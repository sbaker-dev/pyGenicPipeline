from pyGenicPipeline.utils import errors as ec
from pyGenicPipeline.utils import misc as mc
from pyGenicPipeline.core.Input import Input

from miscSupports import terminal_time, flatten
import numpy as np
import time


class FilterSnps(Input):
    def __init__(self, args):
        super().__init__(args)
        # todo: Counts are not considering already filtered snps - change sum to be the sum of filter?
        self._filter_error_dict = {"Filter case": "Count", f"{self.freq} Discrepancy": 0, "MAF": 0,
                                   "Monomorphic": 0, "In Long Range LD": 0}

    def pgs_filter_snps(self):
        """
        Large numbers of snps and individuals can lead to significant memory issues. This will filter the snps in chunks
        vus allowing it to run with less memory
        """

        # Construct the reference panel
        gen_file = self.construct_reference_panel()
        t0 = time.time()

        # Chunk the snps, freqs, and bp positions so we can load raw dosage data in a memory conscious way
        sm_dict = self.sm_dict_from_csv()
        snp_list, chunks = self.chunked_snp_names(sm_dict, chunk_return=True)
        bp_positions = np.array_split(mc.variant_array(self.bp_position.lower(), sm_dict[self.sm_variants]), chunks)
        freqs = np.array_split(sm_dict[self.freq], chunks)

        # Filter each chunk to clean any snps that may be probabilistic
        accepted_snps = [self.filter_snp_chunk(gen_file, snps, f, bp, index, len(snp_list))
                         for index, (snps, f, bp) in enumerate(zip(snp_list, freqs, bp_positions), start=1)]

        mc.filter_array(sm_dict, flatten(accepted_snps), "Filter")
        print(f"Found {len(sm_dict[self.sm_variants])} Snps that passed filtering")

    def filter_snp_chunk(self, gen_file, filter_dict, chromosome):
        """
        From our cleaned summary statistics we can construct a set of information for our validation and reference
        genetic samples. These samples need there raw snps to be loaded, and then we can use these to calculate genetic
        frequencies and standard deviations.

        If we have information on summary frequencies, the user has specified a maf_min or wants to filter snps in long
        range LD then we can filter those out. This process will always filter out monomorphic snps.
        """
        print(f"Filtering chunk {index} out of {total}: {terminal_time()}")
        filter_dict = self._set_filter_dict(gen_file, snps, freqs, base_positions)

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

    def _set_filter_dict(self, gen_file, snps, freqs, base_positions):
        """
        This sets up the filtered dict, which acts similar to our sm_dict, but just for our filter variables. We
        construct this in a method call so that the memory of the raw_snps, that are not required after this is
        constructed, can be garbage collected.
        """
        # Setup the base filter of all snps are valid from our extract chunks
        filter_dict = {self.snp_id: snps, self.freq: freqs, self.bp_position: base_positions,
                       self.filter_key: np.full(len(snps), True)}

        # Load the raw snps from the .bed or .bgen file and use it to isolate statistical information
        raw_snps = self.isolate_raw_snps(gen_file, filter_dict[self.snp_id])
        filter_dict[self.f_std] = np.std(raw_snps, 1, dtype='float32')
        filter_dict[self.f_freq] = np.sum(raw_snps, 1, dtype='float32') / (2 * float(gen_file.iid_count))
        return filter_dict

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

