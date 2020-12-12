from pyGenicPipeline.support.ArgMaker import ArgMaker
from pyGenicPipeline.utils import error_codes as ec
from pyGenicPipeline.utils import misc as mc
from pyGenicPipeline.core.Input import Input

from csvObject import CsvObject
from colorama import Fore
from pathlib import Path
import numpy as np


class LDHerit(Input):
    def __init__(self, args):
        super().__init__(args)

        self._genome_dict = {"Genome-wide Stats": "Values", "Lambda Inflation": 0.0, "Mean LD Score": 0.0,
                             "Genome-Wide Heritability": 0.0}

    def genome_wide_heritability(self):
        """
        Once we have cleaned our data sets, we can store the genome wide data along side with some individual
        chromosome information that does not take the form of lists, in a yaml config file.
        """

        config_dict = {}
        cumulative_ld, sum_sq_beta, total_snps = self._heritability_by_chromosome(config_dict)

        # Check that we successfully found our genome wide attributes
        assert all([total_snps, sum_sq_beta, total_snps]) > 0, ec.all_missing(
            "Total_Snps = Sum_Sq_Beta = Cumulative_LD", "genome_wide_heritability")

        # Calculate genome wide heritability
        average_gw_ld_score = cumulative_ld / float(total_snps)
        chi_square_lambda = np.mean(self.sample_size * sum_sq_beta / float(total_snps))
        gw_h2_ld_score_est = max(0.0001, (max(1.0, float(chi_square_lambda)) - 1.0) /
                                 (self.sample_size * (average_gw_ld_score / total_snps)))
        config_dict["Genome"] = {f"{self.genome_key}_{self.avg_ld}": average_gw_ld_score,
                                 f"{self.genome_key}_{self.herit}": gw_h2_ld_score_est,
                                 "Genome_Description": "Genome-wide Statistics"}

        # Log warnings if applicable then print information to terminal
        if gw_h2_ld_score_est > self.ld_radius / 10.0:
            print(Fore.RED + ec.ld_radius_to_large())
        if gw_h2_ld_score_est > 1:
            print(Fore.RED + ec.heritability_to_large())

        self._genome_dict["Lambda Inflation"] = round(float(chi_square_lambda), 7)
        self._genome_dict["Mean LD Score"] = round(float(average_gw_ld_score), 7)
        self._genome_dict["Genome-Wide Heritability"] = round(float(gw_h2_ld_score_est), 7)
        mc.error_dict_to_terminal(self._genome_dict)

        # Construct config file
        ArgMaker().write_yaml_config_dict(config_dict, self.working_dir, "genome_wide_config")

    def _heritability_by_chromosome(self, config_dict):
        """
        This will calculate the heritability for each chromosome cumulatively the values for the genome-wide
        calculation
        """
        cumulative_ld = sum_sq_beta = total_snps = 0
        for file in mc.directory_iterator(self.clean_directory):
            load_file = CsvObject(Path(self.clean_directory, file), self.cleaned_types, set_columns=True)

            # Isolate the generic information
            chromosome, n_snps, n_iid = self._chromosome_from_load(load_file)

            # Calculate the heritability and average LD at a chromosome level
            heritability, average_ld = self._chromosome_heritability(load_file, chromosome, n_snps)

            # Cumulate ld, snps, and sum square beta
            cumulative_ld += np.sum(load_file.column_data[self.c_ld_score])
            total_snps += n_snps
            sum_sq_beta += np.sum(np.array(load_file.column_data[self.c_beta]) ** 2)

            # Store Values for config file
            chromosome_values = {self.herit: heritability, self.count_snp: n_snps, self.count_iid: n_iid,
                                 self.avg_ld: average_ld, "Description": f"Chromosome {chromosome}"}
            config_dict[chromosome] = chromosome_values
        return cumulative_ld, sum_sq_beta, total_snps

    def _chromosome_from_load(self, load_file):
        """
        We will need to extract the number of individuals, which requires re-parsing in the genetic file. We choose this
        based on the chromosome of the load file.
        """
        # Isolate the chromosome from the first row of data
        chromosome = load_file.row_data[0][self.c_chromosome]
        load_path = str(self.select_file_on_chromosome(chromosome, self.gen_directory, self.gen_type))

        # Isolate the genetic load file according to summary, and use this to isolate number of individuals
        _, core = self.construct_validation(load_path)

        return chromosome, load_file.column_length, core.iid_count

    def compute_ld_scores(self, sm_dict, snp_count, iid_count, ld_dict=False):
        """
        This will calculate the ld scores and create a dict of the ld reference for each snp in our normalised list of
        snps
        """
        ld_scores = np.ones(snp_count)

        if ld_dict:
            ld_dict = {}
        else:
            ld_dict = None

        norm_snps = sm_dict[f"{self.ref_prefix}_{self.norm_snps}"]
        for i, snp in enumerate(norm_snps):
            self._calculate_disequilibrium(i, snp, norm_snps, ld_scores, iid_count, snp_count, ld_dict)

        sm_dict[self.ld_scores] = ld_scores
        sm_dict[self.ld_dict] = ld_dict

    def _calculate_disequilibrium(self, snp_index, current_snp, norm_snps, ld_scores, iid_count, snp_count,
                                  ld_dict=None):
        """
        This will calculate the disequilibrium of the snps in a given radius window
        """

        # Create a window of normalised snps around the current snp with a maximum length of (self.ld_radius * 2) + 1
        snps_in_ld = self.local_values(norm_snps, snp_index, snp_count)

        # Calculate the distance dot product
        distance_dp = np.dot(current_snp, snps_in_ld.T) / iid_count
        if isinstance(ld_dict, dict):
            ld_dict[snp_index] = mc.shrink_r2_matrix(distance_dp, iid_count)

        # calculate the ld score
        r2s = distance_dp ** 2
        ld_scores[snp_index] = np.sum(r2s - ((1 - r2s) / (iid_count - 2)), dtype="float32")

    def _chromosome_heritability(self, load_file, chromosome, snp_count):
        """
        This will calculated the chromosome chi-squared lambda (maths from LDPred), and then take the maximum of 0.0001
        or the computed heritability of

        This chromosomes 1 - chi-sq lambda (or 1 if its less than 1 for reasons that are beyond me)
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
            average_ld = np.mean(load_file.column_data[self.c_ld_score])
            sum_beta_sq = np.sum(np.array(load_file.column_data[self.c_beta]) ** 2)
            chi_sq_lambda = np.mean((self.sample_size * sum_beta_sq) / snp_count)

            herit = max(0.0001, (max(1.0, float(chi_sq_lambda)) - 1) / (self.sample_size * (average_ld / snp_count)))
            return herit, average_ld
