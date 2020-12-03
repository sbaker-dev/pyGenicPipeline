from pyGeneticPipe.utils import misc as mc
from pyGeneticPipe.core.Input import Input
from csvObject import CsvObject
from pathlib import Path
import numpy as np


class LDHerit(Input):
    def __init__(self, args):
        super().__init__(args)

    def compute_ld_scores(self, sm_dict, snp_count, iid_count, ld_dict=False):
        """
        This will calculate the ld scores and create a dict of the ld reference for each snp in our normalised list of
        snps
        """
        ld_scores = np.ones(snp_count)

        if ld_dict:
            ld_dict = {}
        else:
            ld_dict = None

        norm_snps = sm_dict[f"{self.ref_prefix}_{self.norm_snps}"]
        for i, snp in enumerate(norm_snps):
            self._calculate_disequilibrium(i, snp, norm_snps, ld_scores, iid_count, snp_count, ld_dict)

        sm_dict[self.ld_scores] = ld_scores
        sm_dict[self.ld_dict] = ld_dict
        return np.mean(ld_scores)

    def _calculate_disequilibrium(self, snp_index, current_snp, norm_snps, ld_scores, iid_count, snp_count,
                                  ld_dict=None):
        """
        This will calculate the disequilibrium, for when we don't have a genetic map.
        """

        # Create a window of normalised snps around the current snp with a maximum length of (self.ld_radius * 2) + 1
        snps_in_ld = self.local_values(norm_snps, snp_index, snp_count)

        # Calculate the distance dot product
        distance_dp = np.dot(current_snp, snps_in_ld.T) / iid_count
        if isinstance(ld_dict, dict):
            ld_dict[snp_index] = mc.shrink_r2_matrix(distance_dp, iid_count)

        # calculate the ld score
        r2s = distance_dp ** 2
        ld_scores[snp_index] = np.sum(r2s - ((1 - r2s) / (iid_count - 2)), dtype="float32")

    def chromosome_heritability(self, sm_dict, average_ld, chromosome, iid_count, snp_count):
        """
        This will calculated the chromosome chi-squared lambda (maths from LDPred), and then take the maximum of 0.0001
        or the computed heritability of

        This chromosomes chi-sq lambda (or 1 if its less than 1 for reasons that are beyond me)
        ---------------------------------------------------------------------------------------
        Number of samples in the summary stats * (average ld score / number snps)

        :return: estimated heritability (h2 in ldpred)
        """
        if self.heritability_calculated:
            try:
                return self.heritability_calculated[chromosome]
            except KeyError:
                raise Exception("You have said you will provided pre-calculated heritability but failed to find it for"
                                f"chromosome {chromosome}")
        else:
            sum_beta_sq = np.sum(sm_dict[self.beta] ** 2)
            char_chi_sq_lambda = np.mean((iid_count * sum_beta_sq) / snp_count)

            return max(0.0001, (max(1.0, float(char_chi_sq_lambda)) - 1) / (iid_count * (average_ld / snp_count)))

    def _genome_wide_heritability(self):

        cumulative_ld = 0
        sum_sq_beta = 0
        total_snps = 0
        for file in [file for file in mc.directory_iterator(self.working_dir) if "inf" in file]:
            load_file = CsvObject(Path(self.working_dir, file),
                                  [int, int, str, str, str, float, float, float, str, float],
                                  set_columns=True)

            cumulative_ld += np.sum(load_file.column_data[self.w_ld_score])
            total_snps += load_file.column_length
            sum_sq_beta += np.sum(np.array(load_file.column_data[self.w_beta]) ** 2)

        average_gw_ld_score = cumulative_ld / float(total_snps)
        chi_square_lambda = np.mean(self.sample_size * sum_sq_beta / float(total_snps))
        gw_h2_ld_score_est = max(0.0001, (max(1.0, float(chi_square_lambda)) - 1.0) /
                                 (self.sample_size * (average_gw_ld_score / total_snps)))

        print(gw_h2_ld_score_est)
