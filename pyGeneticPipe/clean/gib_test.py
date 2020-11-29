import numpy as np
import h5py
import gzip
import pickle
from scipy import linalg


def _load_ld_info_(local_ld_dict_file, compressed=True):

    if compressed:
        f = gzip.open(local_ld_dict_file, 'r')
    else:
        f = open(local_ld_dict_file, 'r')
    ld_dict = pickle.load(f)
    f.close()
    return ld_dict


def get_chr_heritiablity(cord_data_g, ld_scores_dict, n, tot_betas_sq, tot_num_snps):
    herit_dict = {}
    for chrom_str in cord_data_g:
        print(chrom_str)
        g = cord_data_g[chrom_str]
        betas = g['betas'][...]
        n_snps = len(betas)
        tot_num_snps += n_snps

        sum_beta_sq = np.sum(betas ** 2)
        tot_betas_sq += sum_beta_sq

        # todo number of snps is just the accept snps in Cleaner, could add this along side with adjust ld within
        #  current per chormosome information
        chr_chi_sq_lamda = np.mean((n * sum_beta_sq) / float(n_snps))
        chr_avg_ld_score = ld_scores_dict['chrom_dict'][chrom_str]['avg_ld_score']
        chr_h2_ld_score_est = max(0.0001, (max(1.0, float(chr_chi_sq_lamda)) - 1) / (n * (chr_avg_ld_score / n_snps)))
        herit_dict[chrom_str] = {'n_snps': n_snps, 'h2': chr_h2_ld_score_est}

    # todo this bit requries the whole geneome.
    L = ld_scores_dict['avg_gw_ld_score']
    chi_square_lambda = np.mean(n * tot_betas_sq / float(tot_num_snps))
    gw_h2_ld_score_est = max(0.0001, (max(1, chi_square_lambda) - 1) / (n * (L / tot_num_snps)))
    herit_dict['gw_h2_ld_score_est'] = gw_h2_ld_score_est
    return herit_dict


def main_call(coord_file, n, radius):
    df = h5py.File(coord_file, 'r')
    cord_data_g = df['cord_data']
    ld_dict = _load_ld_info_(r"C:\Users\Samuel\PycharmProjects\External Libaries\Tests\EUR.ld_ldradius183.pkl.gz")

    ld_scores_dict = ld_dict['ld_scores_dict']
    chrom_ld_dict = ld_dict['chrom_ld_dict']
    chrom_ref_ld_mats = ld_dict['chrom_ref_ld_mats']
    mean_n = n

    cord_data_g = df['cord_data']

    tot_num_snps = 0
    tot_betas_sq = 0
    ld_window_size = radius * 2

    herit_dict = get_chr_heritiablity(cord_data_g, ld_scores_dict, n, tot_betas_sq, tot_num_snps)

    for chrom_str in cord_data_g:

        g = cord_data_g[chrom_str]
        beta_hats = g['betas'][...]
        h2 = herit_dict[chrom_str]['h2']

        m = len(beta_hats)
        updated_betas = np.empty(m)

        for i, wi in enumerate(range(0, m, ld_window_size)):
            start_i = wi
            stop_i = min(m, wi + ld_window_size)
            D = chrom_ref_ld_mats[chrom_str][i]
            print(start_i)
            print(stop_i)
            print(D.shape)
            A = ((m / h2) * np.eye(min(m, (wi + (radius * 2))) - wi) + (n / 1.0) * D)
            A_inv = linalg.pinv(A)

            print(A_inv[0])
            break

            # updated_betas[start_i: stop_i] = np.dot(A_inv * n, beta_hats[start_i: stop_i])  # Adjust the beta_hats

            print("")
        break






if __name__ == '__main__':
    data_file = r"C:\Users\Samuel\Documents\Genetic_Examples\PolyTutOut\EUR.coord"
    rr = 183
    oout_file = "EUR.weight"
    ps = [1, 0.3, 0.1, 0.03, 0.01, 0.003, 0.001]
    ns = 253288
    h2 = None
    use_gw_h2 = False
    sampl_var_shrink_factor = 1
    incl_long_range_ld = False
    num_iter = 100
    verbose = True
    zero_jump_prob = 0.01
    burn_in = 10

    main_call(data_file, ns, rr)
