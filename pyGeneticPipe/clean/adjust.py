from pyGeneticPipe.clean.gibs_processor import gibs_processor
import numpy as np
import h5py
import sys


def standardise_snps(g):
    """
    This loads in the snps and then creates a normalised snp.
    :return:
    """
    # todo Currently these three are stored seperatly but are only used so far for a normalised snp?
    raw_snps = g['raw_snps_ref'][...]
    snp_stds = g['snp_stds_ref'][...]
    snp_means = g['snp_means_ref'][...]
    n_snps, n_individuals = raw_snps.shape

    # Need to reformat the shape to construct snps
    snp_means.shape = (n_snps, 1)
    snp_stds.shape = (n_snps, 1)

    snps = np.array((raw_snps - snp_means) / snp_stds, dtype="float32")

    assert snps.shape == raw_snps.shape
    return snps, n_snps, n_individuals


def shrink_r2_matrix(distance_dp, n):
    """
    The values of R squared to be no less or greater than -1 and 1 respectively so the value calculated from the dot
    product is bounded to -1 and 1. Any value that is less than 1 divided by the number of individuals minus 1
    (to account for 0 indexing) is set to zero to reduce the complexity of the finalised r squared matrix

    :param distance_dp: The dot product calculated from the snps in ld and the transposed matrix of snps account for the
        number of individuals in the sample
    :param n: Then umber of individuals in the sample
    :return: Clipped and reduced r2 matrix
    """

    clipped_r2 = np.clip(distance_dp, -1, 1)
    clipped_r2[np.absolute(clipped_r2) < (1.0/(n-1))] = 0
    return clipped_r2


def _calculate_disequilibrium(snps, snp_i, snp, n_individuals, ld_dict, ld_scores):
    """
    This (i belie...the method was just called _calculate_D) calculates the disequilibrium of the current snps in ld

    :param snps: Normalised snps around the current normalised snp with a maximum length of (radius * 2) + 1
    :param snp_i: Current snp index
    :param snp: Current normalised snp
    :param n_individuals: Number of individuals in the sample
    :param ld_dict: dict for storing shrunk r2 matrix
    :param ld_scores: list for storing r2 matrix
    :return: Nothing, append to dict and list of ld_dict and ld_scores
    """

    distance_dp = np.dot(snp, snps.T) / n_individuals

    ld_dict[snp_i] = shrink_r2_matrix(distance_dp, n_individuals)

    r2s = distance_dp ** 2
    ld_scores[snp_i] = np.sum(r2s - ((1 - r2s) / (n_individuals - 2)), dtype='float32')


def snps_in_ld(snps, snp_index, radius, number_of_snps):
    """
    We want to construct a window of -r + r around each snp where r is the radius. However, the first r and last N-r of
    the snps will not have r number of snps before or after them so we need to account for this by:

    Taking the maximum of (0, i-r) so that we never get a negative index
    Taking the minimum of (n_snps, (i + radius + 1)) to ensure we never get an index out of range

    :param snps: Normalised snps
    :param snp_index: Index
    :param radius: radius
    :param number_of_snps: total number of snps

    :return: An array of shape snps of a maximum of 'radius' number of snps surrounding the current snp accessed via
        index.
    """
    return snps[max(0, snp_index - radius): min(number_of_snps, (snp_index + radius + 1))]


def snps_in_window(snps, window_start, number_of_snps, window_size):
    """
    When accessing a window of snps we don't look at the snps +/- radius from the current snp, but just iterate in
    chunks of r*2 through the snps as a window. However, the last iteration may be out of range so we take the minimum
    of the number of snps as a precaution to prevent that.

    :param snps: Normalised snps
    :param window_start: Start index from the window iteration
    :param number_of_snps: total number of snps
    :param window_size: The size, radius * 2, of the window
    :return:
    """

    return snps[window_start: min(number_of_snps, (window_start + window_size))]


def estimate_heritability(h2_calculated, g, chr_avg_ld_score, n, n_snps):
    """
    This will calculated the chromosome chi-squared lambda (maths from LDPred), and then take the maximum of 0.0001 or
    the  computed heritability of

    This chromosomes chi-sq lambda (or 1 if its less than 1 for reasons that are beyond me)
    ---------------------------------------------------------------------------------------
    Number of samples in the summary stats * (average ld score / number snps)

    :param g: temproy load directory. Will be accessed via sm_dict
    :param chr_avg_ld_score: avaerage of the score calcualted in the step before
    :param n: number of samples in summary stats
    :param n_snps: Number of snps that passed screening on this chromosome

    :return: estimated heritaiblity (h2 in ldpred)
    """
    if h2_calculated:
        return h2_calcualted
    else:
        betas = g['betas'][...]

        sum_beta_sq = np.sum(betas ** 2)

        char_chi_sq_lambda = np.mean((n * sum_beta_sq) / float(n_snps))

        return max(0.0001, (max(1.0, float(char_chi_sq_lambda)) - 1) / (n * (chr_avg_ld_score / n_snps)))


def _multiple_hertiaiblity(updated_betas, n, n_snps, ld_scores):
    # todo Do we need both as this is just conjecture at this point?
    sum_beta_sq = np.sum(updated_betas ** 2)

    char_chi_sq_lambda = np.mean((n * sum_beta_sq) / float(n_snps))

    asddd = max(0.0001, (max(1.0, float(char_chi_sq_lambda)) - 1) / (n * (np.mean(ld_scores) / n_snps)))
    print(asddd)


def compute_ld_scores(n_individuals, n_snps, radius, snps):
    ld_dict = {}
    ld_scores = np.ones(n_snps)
    # todo genetic map seems to be a thing, we can look into that later
    # If we don't have a genetic map
    # compute disequilibrium.

    for i, snp in enumerate(snps):
        _calculate_disequilibrium(snps_in_ld(snps, i, radius, n_snps), i, snp, n_individuals, ld_dict, ld_scores)

    # Temporally store these until we know what they do
    ret_dict = {}
    ret_dict["ld_dict"] = ld_dict
    ret_dict["ld_scores"] = ld_scores
    return ld_scores, ld_dict


def infinitesimal_betas(g, h2, n, n_individuals, n_snps, radius, snps):
    """
    Apply the infinitesimal shrink w LD (which requires LD information). LDPRED

    """
    ld_window_size = radius * 2
    beta_hats = g['betas'][...]
    updated_betas = np.empty(n_snps)
    for wi in range(0, n_snps, ld_window_size):
        window_snps = snps_in_window(snps, wi, n_snps, ld_window_size)
        start_i = wi
        stop_i = min(n_snps, wi + ld_window_size)

        D = shrink_r2_matrix(np.dot(window_snps, window_snps.T) / n_individuals, n_individuals)

        # numpy.eye is just an identity matrix
        A = np.array(((n_snps / h2) * np.eye(min(n_snps, (wi + (radius * 2))) - wi)) + (n / 1.0) * D)

        updated_betas[start_i: stop_i] = np.dot(np.linalg.pinv(A) * n, beta_hats[start_i: stop_i])
    return updated_betas


def load_lrld_dict(chromosome):
    # Load Price et al. AJHG 2008 long range LD table.
    long_dict = {chromosome_key: {} for chromosome_key in range(1, 24)}

    # todo set system args to laod data for long-range-ld and hm3 snps
    lrld = r"C:\Users\Samuel\PycharmProjects\External Libaries\ldpred\reference\long-range-ld-price-2008hg38.txt"

    with open(lrld, 'r') as f:
        for line in f:
            chromosome_line, start_pos, end_pos, hild = line.split()
            try:
                long_dict[int(chromosome_line)][hild] = {'start_pos': int(start_pos), 'end_pos': int(end_pos)}
            except ValueError:
                continue

    # todo use chromosome (int via cleaner)
    return long_dict[1]


def filter_long_range(g, chrom_str):
    """
    This will load the information from our long range ld dict and filter out any snps that are in long range LD.
    """
    d = load_lrld_dict(chrom_str)
    positions = g['positions'][...]

    if len(d) != 0:
        for key in d.keys():
            long_filter = np.where((d[key]["start_pos"] < positions) & (positions < d[key]["end_pos"]), True, False)
            # += np.sum(filtered long range)
            # filter call of attributes



def call_main(coord_file, radius, n, ps, filter_long_range_ld, h2_calculated):
    df = h5py.File(coord_file, 'r')
    cord_data_g = df['cord_data']
    tt = 0
    for chrom_str in cord_data_g:
        g = cord_data_g[chrom_str]

        # todo This is just a standardised mesure that we could have calculated in input
        # todo n_snps and n_individuals be accessed via a property esk.
        snps, n_snps, n_individuals = standardise_snps(g)

        # Calculate the ld scores and a dict containing information (don't know what it does currently)
        ld_scores, ld_dict = compute_ld_scores(n_individuals, n_snps, radius, snps)

        # Estimate the partition heritability of this chromosome
        h2 = estimate_heritability(h2_calculated, g, np.mean(ld_scores), n, n_snps)

        # Update the betas via infinitesimal shrinkage using ld information
        updated_betas = infinitesimal_betas(g, h2, n, n_individuals, n_snps, radius, snps)

        # # heritibailtiy on betas post infinitesimal shrink (for testing only)
        # _multiple_hertiaiblity(updated_betas, n, n_snps, ld_scores)

        # todo, we probably want to do this in the filter stage of cleaner before standardisation so that n_snps ==
        #  number of filtered snps
        # Filter long range LD if set
        if filter_long_range_ld:
            filter_long_range(g, chrom_str)

        print(updated_betas)
        print(h2)
        beta_hats = g['betas'][...]

        for cp in ps:

            gibs_processor(n_snps, beta_hats, n, updated_betas, cp)

        print("F")
        break



if __name__ == '__main__':
    cf = r"C:\Users\Samuel\Documents\Genetic_Examples\PolyTutOut\EUR.coord"
    fff = "EUR.ld"
    rr = 183
    ns = 253288
    pss = [1, 0.3, 0.1, 0.03, 0.01, 0.003, 0.001]
    h2_calcualted = None

    call_main(cf, rr, ns, pss, True, h2_calcualted)




