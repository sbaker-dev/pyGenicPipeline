import numpy as np
import h5py


def standardise_snps(g):
    """
    This loads in the snps and then creates a normalised snp.
    :return:
    """
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


def call_main(coord_file, radius):
    df = h5py.File(coord_file, 'r')
    cord_data_g = df['cord_data']

    tt = 0
    for chrom_str in cord_data_g:
        g = cord_data_g[chrom_str]

        # todo This is just a standardised mesure that we could have calculated in input
        snps, n_snps, n_individuals = standardise_snps(g)

        ld_dict = {}
        ld_scores = np.ones(n_snps)

        # todo genetic map seems to be a thing, we can look into that later
        # If we don't have a genetic map
        # compute disequilibrium.
        for i, snp in enumerate(snps):
            _calculate_disequilibrium(snps_in_ld(snps, i, radius, n_snps), i, snp, n_individuals, ld_dict, ld_scores)

        print("E")
        ld_window_size = radius * 2
        ref_ld_matrices = []
        for wi in range(0, n_snps, ld_window_size):
            wi_distance = snps_in_window(snps, wi, n_snps, ld_window_size)
            ref_ld_matrices.append(shrink_r2_matrix(np.dot(wi_distance, wi_distance.T) / n_individuals, n_individuals))

        ret_dict = {}
        ret_dict["ld_dict"] = ld_dict
        ret_dict["ld_scores"] = ld_scores
        ret_dict["ref_ld_matrices"] = ref_ld_matrices

        print(ret_dict)

        break






if __name__ == '__main__':
    cf = r"C:\Users\Samuel\Documents\Genetic_Examples\PolyTutOut\EUR.coord"
    fff = "EUR.ld"
    rr = 183

    call_main(cf, rr)




