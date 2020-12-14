from pyGenicPipeline.utils import errors as ec
from pyGenicPipeline.utils import misc as mc
from pyGenicPipeline.core.Input import Input

from scipy import stats
import numpy as np
import time


class Gibbs(Input):
    def __init__(self, args):
        super().__init__(args)
        self.start_time = 0

    def construct_gibbs_weights(self, sm_dict, load_path, chromosome):
        """
        THis will construct weights for both the infinitesimal and variant fraction models, with the latter being run
        via a Gibbs processor. Both methods from LDPred
        """

        # Check mandatory args are set
        self._assert_construct_gibbs_weights()
        print(f"Constructing Weights for Chromosome {chromosome}")

        # Update the betas via infinitesimal shrinkage using ld information to use as the start value for the variant
        # fraction method used within Gibbs
        inf_betas, stds = self._infinitesimal_betas(sm_dict, load_path)
        print(f"Calculated infinitesimal for chromosome {chromosome} in {round(time.time() - self.start_time, 2)}"
              f" Seconds\n")

        # If the user wants to run the LD pred Gibbs then do so.
        if self.gibbs_run:
            self._gibbs_on_causal_fraction(chromosome, inf_betas, sm_dict)

        # Do the same for the infinitesimal model
        sm_dict[self.inf_dec] = inf_betas / stds.flatten()

    def _gibbs_on_causal_fraction(self, chromosome, inf_betas, sm_dict):
        """If gibbs is set to run, iterate through the causal fractions and calculate the gibbs beta"""
        print("WARNING - DEPRECIATED IN CURRENT VERSION")
        for variant_fraction in self.gibbs_causal_fractions:
            self.start_time = time.time()

            # Run the LDPred gibbs processor to calculate a beta value
            beta = self.gibbs_processor(inf_betas, variant_fraction, sm_dict)

            # Use the genome wide heritability to validate our sum_sq_betas, if we don't have enough causal variants at
            # this level we might end up with a lack of convergence
            sum_sq_beta = np.sum(beta ** 2)
            if sum_sq_beta > self.gm[f"{self.genome_key}_{self.herit}"]:
                ec.gibbs_convergence(variant_fraction, self.gm[self.count_snp], sum_sq_beta,
                                     self.gm[f'{self.genome_key}_{self.herit}'])

                # If True do not save this estimate and stop processing all future estimates, else save regardless
                if self.gibbs_breaker:
                    print("Gibbs_Breaker is turned on, so stopping")
                    break

            # Compute the effect size then write to file
            sm_dict[f"{self.gibbs}_{variant_fraction}"] = beta / sm_dict[self.stds].flatten()
            print(f"Construct weights file for Chromosome {chromosome} variant fraction of {variant_fraction} in "
                  f"{round(time.time() - self.start_time, 2)} Seconds\n")

    def _infinitesimal_betas(self, sm_dict, load_path):
        """
        Apply the infinitesimal shrink w LD (which requires LD information), from LDPred directly (mostly).

        """
        ld_window = self.ld_radius * 2
        snp_count = self.gm[self.count_snp]
        iid_count = self.gm[self.count_iid]

        updated_betas = np.empty(snp_count)
        stds = np.empty(snp_count)
        stds.shape = (snp_count, 1)

        for wi in range(0, snp_count, ld_window):
            print(f"Window {wi} / {snp_count}")
            start_i = wi
            stop_i = min(snp_count, wi + ld_window)
            current_window = stop_i - start_i

            # Load and normalise the snp dosage data for this window
            window_snp_names = self.variant_names(sm_dict)[start_i: min(snp_count, (start_i + ld_window))]
            window_snps, std = self.normalise_snps(self.gen_reference(load_path), window_snp_names, True)

            # Load the disequilibrium (D in LDPred)?
            dis = mc.shrink_r2_matrix(np.dot(window_snps, window_snps.T) / iid_count, iid_count)

            # numpy.eye is just an identity matrix, don't know what A standards for in LDPred
            a = np.array(((snp_count / self.gm[self.herit]) * np.eye(current_window)) + (self.sample_size / 1.0) * dis)

            # Update the betas
            updated_betas[start_i: stop_i] = np.dot(np.linalg.pinv(a) *
                                                    self.sample_size, sm_dict[self.beta][start_i: stop_i])
            stds[start_i: stop_i] = std

        return updated_betas, stds

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
