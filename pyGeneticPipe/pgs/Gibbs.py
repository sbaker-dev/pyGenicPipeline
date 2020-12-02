from pyGeneticPipe.core.Input import Input
from pyGeneticPipe.utils import misc as mc
import numpy as np


class Gibbs(Input):
    def __init__(self, args):
        super().__init__(args)

    def construct_gibbs_weights(self, sm_dict, chromosome):

        # Calculate the mean ld score and construct a dict of ld reference
        average_ld = self.compute_ld_scores(sm_dict)

        # Calculate the estimate heritability for this chromosome
        estimated_herit = self._estimate_heritability(sm_dict, average_ld, chromosome)

        # Update the betas via infinitesimal shrinkage using ld information
        updated_betas = self._infinitesimal_betas(sm_dict, estimated_herit)

        for variant_fraction in self.gibbs_causal_fractions:
            print(variant_fraction)

        return

    def compute_ld_scores(self, sm_dict):
        """
        This will calculate the ld scores and create a dict of the ld reference for each snp in our normalised list of
        snps
        """
        ld_scores = np.ones(sm_dict[f"{self.ref_prefix}_{self.snp_count}"])
        ld_dict = {}

        # todo Allow for genetic map
        if "Genetic_Map" in sm_dict.keys():
            raise NotImplementedError("Genetic Map Not Yet Implemented")
        else:
            norm_snps = sm_dict[f"{self.ref_prefix}_{self.norm_snps}"]
            number_iid = sm_dict[f"{self.ref_prefix}_{self.iid_count}"]
            for i, snp in enumerate(norm_snps):
                self.calculate_disequilibrium(i, snp, norm_snps, sm_dict, ld_scores, ld_dict, number_iid)

        sm_dict[f"{self.ref_prefix}_{self.ld_dict}"] = ld_dict
        return np.mean(ld_scores)

    def calculate_disequilibrium(self, snp_index, current_snp, norm_snps, sm_dict, ld_scores, ld_dict, iid_count):
        """
        This will calculate the disequilibrium, for when we don't have a genetic map.
        """

        # Create a window of normalised snps around the current snp with a maximum length of (self.ld_radius * 2) + 1
        snps_in_ld = self.local_values(norm_snps, snp_index, sm_dict[f"{self.ref_prefix}_{self.snp_count}"])

        # Calculate the distance dot product
        distance_dp = np.dot(current_snp, snps_in_ld.T) / iid_count
        ld_dict[snp_index] = mc.shrink_r2_matrix(distance_dp, iid_count)

        # calculate the ld score
        r2s = distance_dp ** 2
        ld_scores[snp_index] = np.sum(r2s - ((1 - r2s) / (iid_count - 2)), dtype="float32")

    def _estimate_heritability(self, sm_dict, average_ld, chromosome):
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
            iid_count = sm_dict[f"{self.ref_prefix}_{self.iid_count}"]
            snp_count = float(sm_dict[f"{self.ref_prefix}_{self.snp_count}"])

            char_chi_sq_lambda = np.mean((iid_count * sum_beta_sq) / snp_count)

            return max(0.0001, (max(1.0, float(char_chi_sq_lambda)) - 1) / (iid_count * (average_ld / snp_count)))

    def _infinitesimal_betas(self, sm_dict, estimate_herit):
        """
        Apply the infinitesimal shrink w LD (which requires LD information), from LDPred directly (mostly).

        """
        ld_window = self.ld_radius * 2
        iid_count = sm_dict[f"{self.ref_prefix}_{self.iid_count}"]
        snp_count = sm_dict[f"{self.ref_prefix}_{self.snp_count}"]

        updated_betas = np.empty(snp_count)
        for wi in range(0, snp_count, ld_window):
            window_snps = mc.snps_in_window(sm_dict[f"{self.ref_prefix}_{self.norm_snps}"], wi, snp_count, ld_window)

            D = mc.shrink_r2_matrix(np.dot(window_snps, window_snps.T) / iid_count, iid_count)

            # numpy.eye is just an identity matrix
            A = np.array(((snp_count / estimate_herit) * np.eye(min(snp_count, (wi + (self.ld_radius * 2))) - wi))
                         + (iid_count / 1.0) * D)

            # Update the betas for this window
            start_i = wi
            stop_i = min(snp_count, wi + ld_window)
            updated_betas[start_i: stop_i] = np.dot(np.linalg.pinv(A) * iid_count, sm_dict[self.beta][start_i: stop_i])
        return updated_betas

    def local_values(self, values, snp_index, number_of_snps):
        """
        We want to construct a window of -r + r around each a given list of values where r is the radius. However, the
        first r and last N-r of the snps will not have r number of snps before or after them so we need to account for
        this by:

        Taking the maximum of (0, i-r) so that we never get a negative index
        Taking the minimum of (n_snps, (i + radius + 1)) to ensure we never get an index out of range

        :param values: A set of values to extract a local off
        :param snp_index: Index
        :param number_of_snps: total number of snps

        :return: An array of shape snps of a maximum of 'radius' number of snps surrounding the current snp accessed via
            index.
        """
        return values[max(0, snp_index - self.ld_radius): min(number_of_snps, (snp_index + self.ld_radius + 1))]
