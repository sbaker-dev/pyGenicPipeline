import numpy as np
from scipy import stats
import sys
import time


def local_values(values, snp_index, radius, number_of_snps):
    """
    We want to construct a window of -r + r around each a given list of values where r is the radius. However, the first
    r and last N-r of the snps will not have r number of snps before or after them so we need to account for this by:

    Taking the maximum of (0, i-r) so that we never get a negative index
    Taking the minimum of (n_snps, (i + radius + 1)) to ensure we never get an index out of range

    :param values: A set of values to extract a local off
    :param snp_index: Index
    :param radius: radius
    :param number_of_snps: total number of snps

    :return: An array of shape snps of a maximum of 'radius' number of snps surrounding the current snp accessed via
        index.
    """
    return values[max(0, snp_index - radius): min(number_of_snps, (snp_index + radius + 1))]


def const_dict_constructor(chromosome_heritability, cp, n, n_snps):
    causal_variants = n_snps * cp
    hert_by_cv = (chromosome_heritability / causal_variants)
    rv_scalars = np.zeros(n_snps)
    const_dict = {'Mp': causal_variants, 'hdmp': hert_by_cv}

    if n is not None:
        # NOTE M = Number of snps and P is the fraction of causal variants (cv)
        hdmpn = hert_by_cv + 1.0 / n
        hdmp_hdmpn = (hert_by_cv / hdmpn)
        c_const = (cp / np.sqrt(hdmpn))
        d_const = (1.0 - cp) / (np.sqrt(1.0 / n))

        const_dict['hdmpn'] = hdmpn
        const_dict['hdmp_hdmpn'] = hdmp_hdmpn
        const_dict['c_const'] = c_const
        const_dict['d_const'] = d_const
        rv_scalars[:] = sampl_var_shrink_factor * np.sqrt(hdmp_hdmpn * (1.0 / n))

    else:
        raise NotImplementedError("Case control should have been caught by now but if not it isn't implemented")

    const_dict['rv_scalars'] = rv_scalars
    return const_dict


def set_alpha(chromosome_heritability, h2_est, n, tight_sampling, zero_jump_prob):
    if tight_sampling:
        # Force an alpha shrink if estimates are way off compared to heritability estimates.
        # (May improve MCMC convergence.)
        alpha = min(1.0 - zero_jump_prob, 1.0 / h2_est, (chromosome_heritability + 1.0 / np.sqrt(n)) / h2_est)
    else:
        alpha = 1.0 - zero_jump_prob
    return alpha

def get_constants(snp_i,const_dict):
    if 'snp_dict' in const_dict:
        return const_dict['snp_dict'][snp_i]
    else:
        return const_dict


def posterior_mean(cd, b2, n):
    d_const_b2_exp = cd['d_const'] * np.exp(-b2 * n / 2.0)
    numerator = cd['c_const'] * np.exp(-b2 / (2.0 * cd['hdmpn']))
    if not isinstance(d_const_b2_exp, complex):
        if not isinstance(numerator, complex) and (numerator != 0.0):
            postp = numerator / (numerator + d_const_b2_exp)
            assert type(postp) != complex, "Post mean not a real number"
            return postp
        else:
            return 0.0
    else:
        return 1.0


def gibs_processor(n_snps, sample_betas, n, start_betas, cp, chromosome_heritability, number_of_iterations, burn,
                   ld_dict, tight_sampling=False, zero_jump_prob = 0.01):

    # todo Allow seed to be set by user
    # Set random seed to stabilize results
    np.random.seed(42)
    currant_betas = np.copy(start_betas)
    curr_post_means = np.zeros(n_snps)
    avg_betas = np.zeros(n_snps)
    iter_order = np.arange(n_snps)

    # todo In a case_control environment then we can end up with differing N, but otherwise N is just provided by the
    #  user and is statistic

    const_dict = const_dict_constructor(chromosome_heritability, cp, n, n_snps)

    for k in range(number_of_iterations):
        # calculate the chromosome heritability with the current betas
        h2_est = max(0.00001, float(np.sum(currant_betas)))

        # Set alpha for the shrink
        alpha = set_alpha(chromosome_heritability, h2_est, n, tight_sampling, zero_jump_prob)

        rand_ps = np.random.random(n_snps)
        rand_norms = stats.norm.rvs(0.0, 1, size=n_snps) * const_dict['rv_scalars']

        for i, snp_i in enumerate(iter_order):
            # Figure out what sample size and constants to use
            cd = get_constants(snp_i, const_dict)

            # Local (most recently updated) effect estimates
            local_betas = local_values(currant_betas, snp_i, ld_radius, n_snps)
            local_betas[min(ld_radius, snp_i)] = 0.0

            # Calculate the residual of beta hat from the dot product of the local LD matrix and local betas
            res_beta_hat_i = sample_betas[snp_i] - np.dot(ld_dict[snp_i], local_betas)
            b2 = res_beta_hat_i ** 2

            # Calculate the posterior mean p
            postp = posterior_mean(cd, b2, n)

            curr_post_means[snp_i] = cd['hdmp_hdmpn'] * postp * res_beta_hat_i

            if rand_ps[i] < postp * alpha:
                # Sample from the posterior Gaussian dist.
                proposed_beta = rand_norms[snp_i] + cd['hdmp_hdmpn'] * res_beta_hat_i
            else:
                # Sample 0
                proposed_beta = 0.0

            # Update beta for this snp
            currant_betas[snp_i] = proposed_beta

        if k >= burn:
            avg_betas += curr_post_means

    # Averaging over the posterior means instead of samples.
    return avg_betas / float(number_of_iterations - burn)



seed_value = 42
pval_derived_betas = "g[beta]"
# The genome-wide heritability assumed by LDpred, which is then partitioned proportional to the number of SNPs on each
# chromosome which we calculate by default.
h2note = None
ns = 253288
# TODO: NOTE
ni = ns
# ns = "g[ns]"  # An array of n.....
p = [1, 0.3, "etc"]  # The fraction of causal variants used in the gibs sampler
ld_radius = radius = 183
verbose = True
num_inter = 100  # Number of iterations to run the gins sampler for
burn_in_iter = 10  # number of iterations to invalidate
# ld_dict = {}
start_beta = inf_reduced = 0
# boundarys = "genfile"
zero_jump_probability = 0.01  # who knows what the hell this is... its private so ???
sampl_var_shrink_factor = 1 # something to do with the bayes shrink
snp_lrld = None  # We don't need this, we filter out snps that don't meet this criteria
