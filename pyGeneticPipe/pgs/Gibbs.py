from pyGeneticPipe.core.Input import Input
from pyGeneticPipe.utils import misc as mc
from scipy import stats
import numpy as np


class Gibbs(Input):
    def __init__(self, args):
        super().__init__(args)

    def construct_gibbs_weights(self, sm_dict, chromosome):

        # Extract information on count data
        iid_count = sm_dict[f"{self.ref_prefix}_{self.iid_count}"]
        snp_count = sm_dict[f"{self.ref_prefix}_{self.snp_count}"]

        # Calculate the mean ld score and construct a dict of ld reference
        average_ld = self.compute_ld_scores(sm_dict, iid_count)

        # Calculate the estimate heritability for this chromosome
        estimated_herit = self._estimate_heritability(sm_dict, average_ld, chromosome, iid_count, snp_count)

        # Update the betas via infinitesimal shrinkage using ld information
        updated_betas = self._infinitesimal_betas(sm_dict, estimated_herit, iid_count, snp_count)

        for variant_fraction in self.gibbs_causal_fractions:
            # Run the LDPred gibbs processor to calculate a beta value
            beta = self.gibbs_processor(snp_count, iid_count, updated_betas, estimated_herit, variant_fraction, sm_dict)

            # Compute the effect size
            effect_size = beta / sm_dict[f"{self.ref_prefix}_{self.stds}"].flatten()

            self._write(sm_dict, beta, effect_size, chromosome)
            break

        return

    def compute_ld_scores(self, sm_dict, iid_count):
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
            for i, snp in enumerate(norm_snps):
                self.calculate_disequilibrium(i, snp, norm_snps, sm_dict, ld_scores, ld_dict, iid_count)

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

    def _estimate_heritability(self, sm_dict, average_ld, chromosome, iid_count, snp_count):
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

    def _infinitesimal_betas(self, sm_dict, estimate_herit, iid_count, snp_count):
        """
        Apply the infinitesimal shrink w LD (which requires LD information), from LDPred directly (mostly).

        """
        ld_window = self.ld_radius * 2

        updated_betas = np.empty(snp_count)
        for wi in range(0, snp_count, ld_window):
            window_snps = mc.snps_in_window(sm_dict[f"{self.ref_prefix}_{self.norm_snps}"], wi, snp_count, ld_window)

            dis = mc.shrink_r2_matrix(np.dot(window_snps, window_snps.T) / iid_count, iid_count)

            # numpy.eye is just an identity matrix
            a = np.array(((snp_count / estimate_herit) * np.eye(min(snp_count, (wi + (self.ld_radius * 2))) - wi))
                         + (iid_count / 1.0) * dis)

            # Update the betas for this window
            start_i = wi
            stop_i = min(snp_count, wi + ld_window)
            updated_betas[start_i: stop_i] = np.dot(np.linalg.pinv(a) * iid_count, sm_dict[self.beta][start_i: stop_i])
        return updated_betas

    def gibbs_processor(self, snp_count, iid_count, start_betas, est_herit, variant_fraction, sm_dict):
        """LDPred Gibbs Sampler"""
        # Set random seed to stabilize results
        np.random.seed(self.gibbs_random_seed)
        currant_betas = np.copy(start_betas)
        curr_post_means = np.zeros(snp_count)
        avg_betas = np.zeros(snp_count)
        iter_order = np.arange(snp_count)

        ld_dict = sm_dict[f"{self.ref_prefix}_{self.ld_dict}"]
        const_dict = self._const_dict_constructor(est_herit, variant_fraction, iid_count, snp_count)

        for k in range(self.gibbs_iter):
            # calculate the chromosome heritability with the current betas
            h2_est = max(0.00001, float(np.sum(currant_betas)))

            # Set alpha for the shrink
            alpha = self.set_alpha(est_herit, h2_est, iid_count)

            rand_ps = np.random.random(snp_count)
            rand_norms = stats.norm.rvs(0.0, 1, size=snp_count) * const_dict['rv_scalars']

            for i, snp_i in enumerate(iter_order):
                # Figure out what sample size and constants to use
                cd = self._get_constants(snp_i, const_dict)

                # Local (most recently updated) effect estimates
                local_betas = self.local_values(currant_betas, snp_i, snp_count)
                local_betas[min(self.ld_radius, snp_i)] = 0.0

                # Calculate the residual of beta hat from the dot product of the local LD matrix and local betas
                res_beta_hat_i = sm_dict[self.beta][snp_i] - np.dot(ld_dict[snp_i], local_betas)
                b2 = res_beta_hat_i ** 2

                # Calculate the posterior mean p
                postp = mc.posterior_mean(cd, b2, iid_count)

                curr_post_means[snp_i] = cd['hdmp_hdmpn'] * postp * res_beta_hat_i

                if rand_ps[i] < postp * alpha:
                    # Sample from the posterior Gaussian dist.
                    proposed_beta = rand_norms[snp_i] + cd['hdmp_hdmpn'] * res_beta_hat_i
                else:
                    # Sample 0
                    proposed_beta = 0.0

                # Update beta for this snp
                currant_betas[snp_i] = proposed_beta

            if k >= self.gibbs_burn_in:
                avg_betas += curr_post_means

        # Averaging over the posterior means instead of samples.
        return avg_betas / float(self.gibbs_iter - self.gibbs_burn_in)

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

    def _const_dict_constructor(self, chromosome_heritability, cp, iid_count, snp_count):
        """A bunch of constants where constructed in ldpred for the gibbs processor which is duplicated here"""
        causal_variants = snp_count * cp
        herit_by_cv = (chromosome_heritability / causal_variants)
        rv_scalars = np.zeros(snp_count)
        const_dict = {'Mp': causal_variants, 'hdmp': herit_by_cv}

        if iid_count is not None:
            # NOTE M = Number of snps and CP is the fraction of causal variants (CV)
            hdmpn = herit_by_cv + 1.0 / iid_count
            hdmp_hdmpn = (herit_by_cv / hdmpn)
            c_const = (cp / np.sqrt(hdmpn))
            d_const = (1.0 - cp) / (np.sqrt(1.0 / iid_count))

            const_dict['hdmpn'] = hdmpn
            const_dict['hdmp_hdmpn'] = hdmp_hdmpn
            const_dict['c_const'] = c_const
            const_dict['d_const'] = d_const
            rv_scalars[:] = self.gibbs_shrink * np.sqrt(hdmp_hdmpn * (1.0 / iid_count))

        else:
            raise NotImplementedError("Case control should have been caught by now but if not it isn't implemented")

        const_dict['rv_scalars'] = rv_scalars
        return const_dict

    def set_alpha(self, est_herit, h2_est, iid_count):
        """
        This allows a forced alpha shrink if estimates are way off compared to heritability estimates via gibbs tight
         which may improve MCMC convergence
        """
        if self.gibbs_tight:
            return min(1.0 - self.gibbs_zero_jump, 1.0 / h2_est, (est_herit + 1.0 / np.sqrt(iid_count)) / h2_est)
        else:
            return 1.0 - self.gibbs_zero_jump

    def _write(self, sm_dict, gibbs_effect_beta, gibbs_effect_size, chromosome):
        test_out = r"C:\Users\Samuel\Documents\Genetic_Examples\PolyTutOut\Working\TESTOUT.txt"

        bp_positions = mc.variant_array(self.bp_position.lower(), sm_dict[self.sm_variants])
        snp_ids = mc.variant_array(self.snp_id.lower(), sm_dict[self.sm_variants])
        nt1s = mc.variant_array("a1", sm_dict[self.sm_variants])
        nt2s = mc.variant_array("a2", sm_dict[self.sm_variants])

        with open(test_out, "w") as f:
            f.write('chrom    pos    sid    nt1    nt2    raw_beta      gibs_betas     ldpred_beta\n')
            for pos, sid, nt1, nt2, raw_beta, gibbs_beta, ldpred_beta in zip(bp_positions, snp_ids, nt1s, nt2s,
                                                                             sm_dict[self.log_odds], gibbs_effect_beta,
                                                                             gibbs_effect_size):

                f.write('%s    %d    %s    %s    %s    %0.4e    %0.4e   %0.4e\n' % (
                    chromosome, pos, sid, nt1, nt2, raw_beta, gibbs_beta, ldpred_beta))

    @staticmethod
    def _get_constants(snp_i, const_dict):
        """Managing differing specifications of constant"""
        if 'snp_dict' in const_dict:
            return const_dict['snp_dict'][snp_i]
        else:
            return const_dict
