from pyGeneticPipe.utils import error_codes as ec
from pyGeneticPipe.utils import misc as mc
from pyGeneticPipe.core.Input import Input
from csvObject import write_csv
from decimal import Decimal
from scipy import stats
import numpy as np
import time

import sys


class Gibbs(Input):
    def __init__(self, args):
        super().__init__(args)
        self.start_time = 0

    def construct_gibbs_weights(self, sm_dict, chromosome):

        # Check mandatory args are set
        self._assert_construct_gibbs_weights()
        print(f"Constructing Weights for Chromosome {chromosome}")

        # Update the betas via infinitesimal shrinkage using ld information to use as the start value for the variant
        # fraction method used within Gibbs
        inf_betas = self._infinitesimal_betas(sm_dict)
        print(inf_betas)
        print(f"Calculated LD and chromosome heritability for chromosome {chromosome} in "
              f"{round(time.time() - self.start_time, 2)} Seconds")

        for variant_fraction in self.gibbs_causal_fractions:
            self.start_time = time.time()

            # Run the LDPred gibbs processor to calculate a beta value
            beta = self.gibbs_processor(inf_betas, variant_fraction, sm_dict)

            print(beta)

            sum_sq_beta = np.sum(beta ** 2)

            print(sum_sq_beta)
        #     if sum_sq_beta > self.gm[f"{self.genome_key}_{self.herit}"]:
        #         print(f"Warning: Sum Squared beta is much large than estimated hertiability suggesting a lack of "
        #               f"convergence of Gibbs\n"
        #               f"{sum_sq_beta} > {self.gm[f'{self.genome_key}_{self.herit}']}")
        #
        #     # Compute the effect size then write to file
        #     print(sum_sq_beta)
        #     print(beta)
        #
            effect_size = beta / sm_dict[f"{self.ref_prefix}_{self.stds}"].flatten()
            print(effect_size)
            # self._write_weights(sm_dict, effect_size, chromosome, variant_fraction, beta)
            sys.exit()

        # Do the same for the infinitesimal model
        sm_dict[self.inf_dec] = inf_betas / sm_dict[f"{self.ref_prefix}_{self.stds}"].flatten()

    def _infinitesimal_betas(self, sm_dict):
        """
        Apply the infinitesimal shrink w LD (which requires LD information), from LDPred directly (mostly).

        """
        ld_window = self.ld_radius * 2
        snp_count = self.gm[self.count_snp]
        iid_count = self.gm[self.count_iid]

        updated_betas = np.empty(snp_count)
        for wi in range(0, snp_count, ld_window):
            start_i = wi
            stop_i = min(snp_count, wi + ld_window)
            current_window = stop_i - start_i

            # Load the snps in this window
            window_snps = mc.snps_in_window(sm_dict[f"{self.ref_prefix}_{self.norm_snps}"], wi, snp_count, ld_window)

            # Load the disequilibrium (D in LDPred)?
            dis = mc.shrink_r2_matrix(np.dot(window_snps, window_snps.T) / iid_count, iid_count)

            # numpy.eye is just an identity matrix, don't know what A standards for in LDPred
            a = np.array(((snp_count / self.gm[self.herit]) * np.eye(current_window)) + (self.sample_size / 1.0) * dis)

            # Update the betas
            updated_betas[start_i: stop_i] = np.dot(np.linalg.pinv(a) * self.sample_size,
                                                    sm_dict[self.beta][start_i: stop_i])

        return updated_betas

    def gibbs_processor(self, start_betas, variant_fraction, sm_dict):
        """LDPred Gibbs Sampler"""
        # Set random seed to stabilize results
        np.random.seed(self.gibbs_random_seed)
        currant_betas = np.copy(start_betas)
        curr_post_means = np.zeros(self.gm[self.count_snp])
        avg_betas = np.zeros(self.gm[self.count_snp])
        iter_order = np.arange(self.gm[self.count_snp])

        const_dict = self._const_dict_constructor(variant_fraction)

        for k in range(self.gibbs_iter):
            # calculate the chromosome heritability with the current betas
            h2_est = max(0.00001, float(np.sum(currant_betas ** 2)))

            # Set alpha for the shrink
            alpha = self._set_alpha(self.gm[self.herit], h2_est)

            rand_ps = np.random.random(self.gm[self.count_snp])
            rand_norms = stats.norm.rvs(0.0, 1, size=self.gm[self.count_snp]) * const_dict['rv_scalars']

            for i, snp_i in enumerate(iter_order):
                # Figure out what sample size and constants to use
                cd = self._get_constants(snp_i, const_dict)

                # Local (most recently updated) effect estimates
                local_betas = self.local_values(currant_betas, snp_i, self.gm[self.count_snp])
                local_betas[min(self.ld_radius, snp_i)] = 0.0

                # Calculate the residual of beta hat from the dot product of the local LD matrix and local betas
                res_beta_hat_i = sm_dict[self.beta][snp_i] - np.dot(sm_dict[self.ld_dict][snp_i], local_betas)
                b2 = res_beta_hat_i ** 2

                # Calculate the posterior mean p
                postp = mc.posterior_mean(cd, b2, self.sample_size)

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

    def _const_dict_constructor(self, variant_fraction):
        """A bunch of constants where constructed in ldpred for the gibbs processor which is duplicated here"""
        causal_variants = self.gm[self.count_snp] * variant_fraction
        herit_by_cv = (self.gm[self.herit] / causal_variants)
        rv_scalars = np.zeros(self.gm[self.count_snp])
        const_dict = {'Mp': causal_variants, 'hdmp': herit_by_cv}

        if self.sample_size is not None:
            # NOTE M = Number of snps and CP is the fraction of causal variants (CV)
            hdmpn = herit_by_cv + 1.0 / self.sample_size
            hdmp_hdmpn = (herit_by_cv / hdmpn)
            c_const = (variant_fraction / np.sqrt(hdmpn))
            d_const = (1.0 - variant_fraction) / (np.sqrt(1.0 / self.sample_size))

            const_dict['hdmpn'] = hdmpn
            const_dict['hdmp_hdmpn'] = hdmp_hdmpn
            const_dict['c_const'] = c_const
            const_dict['d_const'] = d_const
            rv_scalars[:] = self.gibbs_shrink * np.sqrt(hdmp_hdmpn * (1.0 / self.sample_size))

        else:
            raise NotImplementedError("Case control should have been caught by now but if not it isn't implemented")

        const_dict['rv_scalars'] = rv_scalars
        return const_dict

    def _set_alpha(self, est_herit, h2_est):
        """
        This allows a forced alpha shrink if estimates are way off compared to heritability estimates via gibbs tight
         which may improve MCMC convergence
        """
        if self.gibbs_tight:
            return min(1.0 - self.gibbs_zero_jump, 1.0 / h2_est, (est_herit + 1.0 / np.sqrt(self.sample_size)) / h2_est)
        else:
            return 1.0 - self.gibbs_zero_jump

    def _write_weights(self, sm_dict, effect_size, chromosome, name, gibbs_beta):
        """
        This will format all of our data into a csv file and store it in the working directory
        """
        # todo This might be generalisable to scores as well if placed within Input
        # Load name based on type
        if isinstance(name, (float, int, np.int8, np.float32)):
            file_name = f"{chromosome}_weights_p{Decimal(name):.2E}"
        else:
            file_name = f"{chromosome}_weights_{name}"

        if not isinstance(gibbs_beta, np.ndarray):
            gibbs_beta = ["NA" for _ in range(len(sm_dict[self.sm_variants]))]

        # Slice variant arrays into lists
        bp_positions = mc.variant_array(self.bp_position.lower(), sm_dict[self.sm_variants])
        snp_ids = mc.variant_array(self.snp_id.lower(), sm_dict[self.sm_variants])
        nt1s = mc.variant_array("a1", sm_dict[self.sm_variants])
        nt2s = mc.variant_array("a2", sm_dict[self.sm_variants])

        # Construct a list of lists, where sub lists represent the rows in the csv file
        rows = [[chromosome, pos, sid, nt1, nt2, beta, log_odds, ld_score, gibbs, ldpred_beta]
                for pos, sid, nt1, nt2, beta, log_odds, ld_score, gibbs, ldpred_beta in zip(
                bp_positions, snp_ids, nt1s, nt2s, sm_dict[self.beta], sm_dict[self.log_odds],
                sm_dict[self.ld_scores], gibbs_beta, effect_size)]

        # Write the file
        write_csv(self.working_dir, file_name, self.gibbs_headers, rows)
        print(f"Construct weights file for Chromosome {chromosome} for {file_name} in "
              f"{round(time.time() - self.start_time, 2)} Seconds")

    @staticmethod
    def _get_constants(snp_i, const_dict):
        """Managing differing specifications of constant"""
        if 'snp_dict' in const_dict:
            return const_dict['snp_dict'][snp_i]
        else:
            return const_dict

    def _assert_construct_gibbs_weights(self):
        """Check mandatory args"""
        assert self.working_dir, ec.missing_arg(self.operation, "Working_Directory")
        assert self.ld_radius, ec.missing_arg(self.operation, "LD_Radius")

        self.start_time = time.time()
