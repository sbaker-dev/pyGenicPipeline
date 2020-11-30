import numpy as np
from scipy import stats
import sys


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


def posterior_mean(d_const_b2_exp, numerator):
    if np.isreal(d_const_b2_exp):
        if np.isreal(numerator):
            if numerator == 0.0:
                return 0.0
            else:
                postp = numerator / (numerator + d_const_b2_exp)
                assert np.isreal(
                    postp), 'The posterior mean is not a real number?  Possibly due to problems with summary stats, LD estimates, or parameter settings.'
                return postp
        else:
            return 0.0
    else:
        print("HERE")
        return 1.0


def gibs_processor(n_snps, sample_betas, n, start_betas, cp, chromosome_heritability, number_of_iterations, burn,
                   ld_dict, tight_sampling=False, zero_jump_prob = 0.01):

    # todo Allow seed to be set by user
    # Set random seed to stabilize results
    np.random.seed(42)
    currant_betas = np.copy(start_betas)
    post_means = np.zeros(n_snps)
    avg_betas = np.zeros(n_snps)

    # todo ?! What is this
    iter_order = np.arange(n_snps)


    # todo In a case_control environment then we can end up with differing N, but otherwise N is just provided by the
    #  user and is statistic

    # todo, do when initialised and then call?
    const_dict = const_dict_constructor(chromosome_heritability, cp, n, n_snps)

    for k in range(number_of_iterations):
        # calculate the chromosome heritability with the current betas
        h2_est = max(0.00001, float(np.sum(currant_betas)))

        # Set alpha for the shrink
        alpha = set_alpha(chromosome_heritability, h2_est, n, tight_sampling, zero_jump_prob)

        rand_ps = np.random.random(n_snps)
        rand_norms = stats.norm.rvs(0.0, 1, size=n_snps) * const_dict['rv_scalars']

        cbs = []
        for i, snp_i in enumerate(iter_order):
            start_i = max(0, snp_i - ld_radius)
            focal_i = min(ld_radius, snp_i)
            stop_i = min(n_snps, snp_i + ld_radius + 1)

            local_betas = currant_betas[start_i: stop_i].copy()
            local_betas[focal_i] = 0.0
            cbs.append(np.array(local_betas))

        # calculate posterior mean
        local_beta_array = np.array(cbs, dtype=object)
        ld_array = np.array([ld_dict[key] for key in ld_dict], dtype=object)
        dot_product = np.array([np.dot(ld, beta) for ld, beta in zip(ld_array, local_beta_array)])


        print(ld_array[1])
        print(local_beta_array[1])


        print("SAMPLE")
        print(sample_betas[0])
        print(sample_betas[1])

        print("DP")
        print(dot_product[0])
        print(dot_product[1])

        res_beta_hat = np.array(sample_betas - dot_product)
        b2 = np.array(res_beta_hat ** 2)
        print("RES")
        print(res_beta_hat[0])
        print(res_beta_hat[1])

        sys.exit()

        # deominator
        d_const_b2_exp = const_dict['d_const'] * np.exp(-b2 * n / 2.0)
        numerator = const_dict["c_const"] * np.exp(-b2 / (2.0 * const_dict['hdmpn']))

        # Calculate the post means and beta
        post_means_array = np.array([posterior_mean(d_c, n_c) for d_c, n_c in zip(d_const_b2_exp, numerator)])
        post_means = np.array(const_dict["hdmp_hdmpn"] * post_means_array * res_beta_hat)
        currant_betas = np.where(rand_ps < post_means_array * alpha, rand_norms + const_dict["hdmp_hdmpn"] *
                                 res_beta_hat, 0.0)

        # print(post_means[0])
        # print(post_means[1])
        # print(currant_betas[0])
        # print(currant_betas[1])

        print("")
        sys.exit()






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
num_inter = 60  # Number of iterations to run the gins sampler for
burn_in_iter = 5  # number of iterations to invalidate
# ld_dict = {}
start_beta = inf_reduced = 0
# boundarys = "genfile"
zero_jump_probability = 0.01  # who knows what the hell this is... its private so ???
sampl_var_shrink_factor = 1 # something to do with the bayes shrink
snp_lrld = None  # We don't need this, we filter out snps that don't meet this criteria
