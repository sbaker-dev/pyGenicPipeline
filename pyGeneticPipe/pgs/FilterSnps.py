from pyGeneticPipe.core.Input import Input
from pyGeneticPipe.utils import misc as mc
import numpy as np


class FilterSnps(Input):
    def __init__(self, args):
        super().__init__(args)
        self._error_dict = {"Filter case": "Count", "Frequency": 0, "MAF": 0,
                            "Monomorphic": 0, "Long_Range": 0}

    def _filter_snps(self, gen_type, genetic_file, sm_dict, chromosome):
        # testing validation for comparision TEMP

        # Construct the genetic raw snps and genetics freqs
        sm_dict[f"{gen_type}_Raw_Snps"] = self._isolate_raw_snps(genetic_file, sm_dict)
        sm_dict[f"{gen_type}_Freqs"] = np.sum(sm_dict[f"{gen_type}_Raw_Snps"], 1, dtype='float32') / (2 * float(genetic_file.iid_count))

        # If the frequencies in the summary stats are not just a list of -1 errors then screen genetic snp frequencies
        if (self.freq_discrepancy < 1) and (np.sum(sm_dict[self.frequency] == -1) != len(sm_dict[self.frequency])):
            freq_filter = np.logical_or(
                np.absolute(sm_dict[self.frequency] - sm_dict[f"{gen_type}_Freqs"]) < self.freq_discrepancy,
                np.absolute(sm_dict[self.frequency] + (sm_dict[f"{gen_type}_Freqs"] - 1)) < self.freq_discrepancy)

            # Invalid frequencies from summary stats where coded as -1 so these will be removed
            freq_filter = np.logical_or(freq_filter, sm_dict[self.frequency] <= 0)
            self._error_dict["Filtered_Frequency"] = len(freq_filter) - np.sum(freq_filter)
            if not mc.filter_array(sm_dict, freq_filter):
                return None

        # Filter minor allele frequency SNPs.
        if self.maf_min > 0:
            maf_filter = (sm_dict[f"{gen_type}_Freqs"] > self.maf_min) * (sm_dict[f"{gen_type}_Freqs"] < (1 - self.maf_min))
            self._error_dict["Filtered_MAF"] = len(maf_filter) - np.sum(maf_filter)
            if not mc.filter_array(sm_dict, maf_filter):
                return None

        # Do the same for std
        stds = np.std(sm_dict[f"{gen_type}_Raw_Snps"], 1, dtype='float32')
        monomorphic_filter = stds > 0
        self._error_dict["Monomorphic"] = len(monomorphic_filter) - np.sum(monomorphic_filter)
        if not mc.filter_array(sm_dict, monomorphic_filter):
            return None

        # Filter long range LD if set
        if self.lr_ld_path:
            lr_pos = self.load_lr_ld_dict()[chromosome]
            if len(lr_pos) != 0:
                for key in lr_pos.keys():
                    position = mc.variant_array(self.bp_position.lower(), sm_dict[self.sm_variants])
                    long_filter = np.where((lr_pos[key]["start_pos"] < position) & (position < lr_pos[key]["end_pos"]),
                                           False, True)
                    self._error_dict["Long_Range"] += len(long_filter) - np.sum(long_filter)
                    if not mc.filter_array(sm_dict, long_filter):
                        return None

        # set accept snps for gibbs
        print(len(sm_dict[self.sm_variants]))
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
