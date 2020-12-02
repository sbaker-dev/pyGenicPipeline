from pyGeneticPipe.core.Input import Input
from pyGeneticPipe.utils import misc as mc
import numpy as np


class Gibbs(Input):
    def __init__(self, args):
        super().__init__(args)

    def construct_gibbs_weights(self, sm_dict):

        # Calculate the mean ld score and construct a dict of ld reference
        self.compute_ld_scores(sm_dict)


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
            norm_snps = sm_dict[f"{self.ref_prefix}_{self.normalised_snps}"]
            number_iid = sm_dict[f"{self.ref_prefix}_{self.iid_count}"]
            for i, snp in enumerate(norm_snps):
                self.calculate_disequilibrium(i, snp, norm_snps, sm_dict, ld_scores, ld_dict, number_iid)

        sm_dict[f"{self.ref_prefix}_{self.ld_scores}"] = np.mean(ld_scores)
        sm_dict[f"{self.ref_prefix}_{self.ld_dict}"] = ld_dict

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

    def local_values(self, values, snp_index, number_of_snps):
        """
        We want to construct a window of -r + r around each a given list of values where r is the radius. However, the first
        r and last N-r of the snps will not have r number of snps before or after them so we need to account for this by:

        Taking the maximum of (0, i-r) so that we never get a negative index
        Taking the minimum of (n_snps, (i + radius + 1)) to ensure we never get an index out of range

        :param values: A set of values to extract a local off
        :param snp_index: Index
        :param number_of_snps: total number of snps

        :return: An array of shape snps of a maximum of 'radius' number of snps surrounding the current snp accessed via
            index.
        """
        return values[max(0, snp_index - self.ld_radius): min(number_of_snps, (snp_index + self.ld_radius + 1))]
