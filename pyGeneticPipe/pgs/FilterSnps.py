from pyGeneticPipe.utils import error_codes as ec
from pyGeneticPipe.utils import misc as mc
from pyGeneticPipe.core.Input import Input
import numpy as np
import time


class FilterSnps(Input):
    def __init__(self, args):
        super().__init__(args)
        self._filter_error_dict = {"Filter case": "Count", self.freq: 0, "MAF": 0, "Monomorphic": 0, "Long_Range": 0}

    def filter_snps(self, gen_type, gen_file, sm_dict, chromosome):
        """
        From our cleaned summary statistics we can construct a set of information for our validation and reference
        genetic samples. These samples need there raw snps to be loaded, and then we can use these to calculate genetic
        frequencies and standard deviations.

        If we have information on summary frequencies, the user has specified a maf_min or wants to filter snps in long
        range LD then we can filter those out. This process will always filter out monomorphic snps.
        """

        t0 = self._assert_filter_snps()
        print(f"\nStarting the filtering for {gen_type}")

        # Load the raw snps from the .bed or .bgen file and use it to isolate statistical information
        raw_snps = self._isolate_raw_snps(gen_file, sm_dict)
        sm_dict[f"{gen_type}_{self.mean}"] = np.mean(raw_snps, 1, dtype='float32')
        sm_dict[f"{gen_type}_{self.stds}"] = np.std(raw_snps, 1, dtype='float32')
        sm_dict[f"{gen_type}_{self.freq}"] = np.sum(raw_snps, 1, dtype='float32') / (2 * float(gen_file.iid_count))
        sm_dict[f"{gen_type}_{self.raw_snps}"] = raw_snps

        # If the frequencies in the summary stats are not just a list of -1 errors then screen genetic snp frequencies
        if (self.freq_discrepancy < 1) and (np.sum(sm_dict[self.freq] == -1) != len(sm_dict[self.freq])):
            if not self._summary_frequencies(sm_dict, gen_type):
                return None

        # Filter minor allele frequency SNPs.
        if self.maf_min > 0:
            if not self._maf_filter(sm_dict, gen_type):
                return None

        # Filter any Monomorphic snps
        if not self._monomorphic_filter(sm_dict, gen_type):
            return None

        # Filter long range LD if set
        if self.lr_ld_path:
            sm_dict = self._long_range_ld_filter(sm_dict, chromosome)
            if not sm_dict:
                return None

        # Now that the genetic information has been cleaned and filtered, use this information to create a normalised
        # snp and clean up sm_dict of anything we know long need
        self._normalise_snps(sm_dict, gen_type)

        t1 = mc.error_dict_to_terminal(self._filter_error_dict)
        print(f"Cleaned summary stats for Chromosome {chromosome} in {round(t1 - t0, 2)} Seconds\n")
        return sm_dict

    def _isolate_raw_snps(self, gen_file, sm_dict):
        """
        This will isolate the raw snps for a given bed or bgen file

        :param gen_file: Genetic file you wish to load from
        :param sm_dict: dict of clean information
        :return: raw snps
        """
        ordered_common = gen_file[:, gen_file.sid_to_index(self._extract_variant_name(sm_dict))].read().val

        # bed returns 2, 1, 0 rather than 0, 1, 2 although it says its 0, 1, 2; so this inverts it
        if self.load_type == ".bed":
            return np.array([abs(snp - 2) for snp in ordered_common.T])

        # We have a [1, 0, 0], [0, 1, 0], [0, 0, 1] array return for 0, 1, 2 respectively. So if we multiple the arrays
        # by their index position and then sum them we get [0, 1, 2]
        elif self.load_type == ".bgen":
            return sum(np.array([snp * i for i, snp in enumerate(ordered_common.T)]))

        else:
            raise Exception(f"Critical Error: Unknown load type {self.load_type} found in _isolate_dosage")

    def _extract_variant_name(self, sm_dict):
        """
        Different file types have different naming standards.

        .bed: ["rs123", "rs124", ... "rsN"]
        .bgen: ["rs123,rs123", "rs124,rs124", ... "rsN,rsN"]

        This will standardise the names to be a list of type equivalent to bed
        :param sm_dict: dict of clean information
        :return: list of snp names
        """
        if self.load_type == ".bed":
            return [variant.variant_id for variant in sm_dict[self.sm_variants]]
        elif self.load_type == ".bgen":
            print("Bgen load type, so need to restructure return type ... will take a bit longer longer!")
            return [variant.bgen_variant_id() for variant in sm_dict[self.sm_variants]]
        else:
            raise Exception(f"Critical Error: Unknown load type {self.load_type} found in _isolate_dosage")

    def _summary_frequencies(self, sm_dict, gen_type):
        """
        If the summary frequencies existed in the summary stats then cross check them with our genetic frequencies,
        removing anything that is outside of the frequency discrepancy allow by the user.
        """
        freq_filter = np.logical_or(
            np.absolute(sm_dict[self.freq] - sm_dict[f"{gen_type}_{self.freq}"]) < self.freq_discrepancy,
            np.absolute(sm_dict[self.freq] + (sm_dict[f"{gen_type}_{self.freq}"] - 1)) < self.freq_discrepancy)

        # Invalid frequencies from summary stats where coded as -1 so these should be removed
        freq_filter = np.logical_or(freq_filter, sm_dict[self.freq] <= 0)
        self._filter_error_dict[self.freq] = len(freq_filter) - np.sum(freq_filter)
        return mc.filter_array(sm_dict, freq_filter)

    def _maf_filter(self, sm_dict, gen_type):
        """
        If the maf filter is greater than zero, then use this to remove any maf frequencies less that the frequency
        provided or greater than 1 - frequency provided.
        """
        maf_filter = (sm_dict[f"{gen_type}_{self.freq}"] > self.maf_min) * \
                     (sm_dict[f"{gen_type}_{self.freq}"] < (1 - self.maf_min))

        self._filter_error_dict["MAF"] = len(maf_filter) - np.sum(maf_filter)
        return mc.filter_array(sm_dict, maf_filter)

    def _monomorphic_filter(self, sm_dict, gen_type):
        """Remove any Monomorphic snps, those with no variation"""
        monomorphic_filter = sm_dict[f"{gen_type}_{self.stds}"] > 0
        self._filter_error_dict["Monomorphic"] = len(monomorphic_filter) - np.sum(monomorphic_filter)
        return mc.filter_array(sm_dict, monomorphic_filter)

    def _long_range_ld_filter(self, sm_dict, chromosome):
        """
        If we have filtering of long range LD in place then we load a dict of long range ld from Price et all 2008 as
        was done within, and data taken from, LdPred. Then we isolate the long range ld regions for this chromosome,
        and filter any in long range LD.
        """
        lr_pos = self.load_lr_ld_dict()[chromosome]
        if len(lr_pos) != 0:
            for key in lr_pos.keys():
                position = mc.variant_array(self.bp_position.lower(), sm_dict[self.sm_variants])
                long_filter = np.where((lr_pos[key]["start_pos"] < position) & (position < lr_pos[key]["end_pos"]),
                                       False, True)
                self._filter_error_dict["Long_Range"] += len(long_filter) - np.sum(long_filter)
                if not mc.filter_array(sm_dict, long_filter):
                    return None

        return sm_dict

    def _normalise_snps(self, sm_dict, gen_type):
        """For gibbs we use normalised snps, this process will use the information we have filtered to construct it"""

        # Get the number of snps and individuals in the filtered dict
        n_snps, n_individuals = sm_dict[f"{gen_type}_{self.raw_snps}"].shape

        # Need to reformat the shape to construct the normalised snps
        raw_means = sm_dict[f"{gen_type}_{self.mean}"]
        raw_stds = sm_dict[f"{gen_type}_{self.stds}"]
        raw_means.shape = (n_snps, 1)
        raw_stds.shape = (n_snps, 1)

        # Use this information to construct a normalised snps
        normalised_snps = np.array((sm_dict[f"{gen_type}_{self.raw_snps}"] - raw_means) / raw_stds, dtype="float32")
        assert normalised_snps.shape == sm_dict[f"{gen_type}_{self.raw_snps}"].shape
        sm_dict[f"{gen_type}_Normalised_Snps"] = normalised_snps

    def _assert_filter_snps(self):
        """Different files require different load type operations, so load type must be set"""
        assert self.load_type, ec.missing_arg(self.operation, "Load_Type")
        return time.time()
