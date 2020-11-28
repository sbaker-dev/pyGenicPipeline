from pathlib import Path
import h5py

import numpy as np


def standardise_snps(g):
    """
    This loads in the snps and then creates a normalised snp.
    :return:
    """
    raw_snps = g['raw_snps_ref'][...]
    snp_stds = g['snp_stds_ref'][...]
    snp_means = g['snp_means_ref'][...]
    n_snps = raw_snps.shape[0]

    # Need to reformat the shape to construct snps
    snp_means.shape = (n_snps, 1)
    snp_stds.shape = (n_snps, 1)

    snps = np.array((raw_snps - snp_means) / snp_stds, dtype="float32")

    assert snps.shape == raw_snps.shape
    return snps, n_snps


def call_main(coord_file):
    df = h5py.File(coord_file, 'r')
    cord_data_g = df['cord_data']

    for chrom_str in cord_data_g:
        g = cord_data_g[chrom_str]

        # todo This is just a standardised mesure that we could have calculated in input
        snps, n_snps = standardise_snps(g)

        break


if __name__ == '__main__':
    cf = r"C:\Users\Samuel\Documents\Genetic_Examples\PolyTutOut\EUR.coord"
    fff = "EUR.ld"
    rr = 183

    call_main(cf)




