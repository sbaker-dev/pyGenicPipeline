from pyGenicPipeline.utils.ArgMaker import ArgMaker
from pyGenicPipeline.utils import errors as ec
from pyGenicPipeline.utils import misc as mc
from pyGenicPipeline.core.Input import Input

from miscSupports import directory_iterator, terminal_time, write_pickle, load_pickle
from csvObject import CsvObject
from colorama import Fore
from pathlib import Path
import numpy as np
import time


class LDHerit(Input):
    def __init__(self, args):
        super().__init__(args)

        self._genome_dict = {"Genome-wide Stats": "Values", "Lambda Inflation": 0.0, "Mean LD Score": 0.0,
                             "Genome-Wide Heritability": 0.0}

    def suggest_ld_radius(self):
        """Suggest the size of LD that the user should be using"""
        total_snps = sum([CsvObject(Path(self.filter_directory, file)).column_length
                          for file in directory_iterator(self.filter_directory)])
        print(f"Suggested LD Radius based on total snps found after filtering / 3000 is {total_snps / 3000}")
        return total_snps / 3000

    def pgs_calculate_ld(self):
        """
        This will calculate the chromosome specific LD information which will be used later to construct genome_wide LD
        information, whilst also saving the values of the normalised Snps and standard errors

        :return: Dict containing the LD dict, scores, matrix, the normalised snps and the standard errors of the raw
            snps
        :rtype: dict
        """
        # Re load the filtered snps into memory
        t0 = time.time()
        sm_dict = self.sm_dict_from_csv(self.filter_directory, f"Filtered_{self.target_chromosome}.csv")

        # Use the snps we found to filter down our reference panel, create the normalised snps, set the sid/iid counts
        ref = self.construct_reference_panel()
        normalised_snps, std = self.normalise_snps(ref, self.snp_names(sm_dict), True)
        sid_count, iid_count = normalised_snps.shape

        # Calculate the Ld scores and ld reference dict
        ld_dict, ld_scores = self._calculate_ld_scores(normalised_snps, sid_count, iid_count)

        # Write the information as a combined dict to pickle file
        write_ld = {"LD_Dict": ld_dict, "LD_Scores": ld_scores, "Norm_Snps": normalised_snps, "Snp_Std": std}
        write_pickle(self.ld_directory, f"LD{self.target_chromosome}", write_ld)

        print(f"Constructed LD for chromosome {self.target_chromosome} in {time.time() - t0} seconds")
        return write_ld

    def _calculate_ld_scores(self, normalised_snps, sid_count, iid_count):
        """This will calculate the Linkage Disequilibrium """

        # Setup arrays and dicts, then iterate though the snps
        ld_scores = np.ones(sid_count)
        ld_dict = {}
        for i, snp in enumerate(normalised_snps):
            # Isolate the local snps in LD
            snps_in_ld = self.local_values(normalised_snps, i, sid_count)

            # Calculate the distance dot product, store the clipped r2 value in the ld_dict
            distance_dp = np.dot(snp, snps_in_ld.T) / iid_count
            ld_dict[i] = mc.shrink_r2_matrix(distance_dp, iid_count)

            # calculate the ld score and save it to the array
            r2s = distance_dp ** 2
            ld_scores[i] = np.sum(r2s - ((1 - r2s) / (iid_count - 2)), dtype="float32")

        print(f"Constructed LD scores {terminal_time()}")
        return ld_dict, ld_scores

    def distribute_heritability_genome_wide(self):
        """If we can't calculate heritability, distribute it from a provided float"""
        total_snps = 0
        config_dict = {}
        for file in directory_iterator(self.summary_directory):
            print(file)
            load_file = CsvObject(Path(self.summary_directory, file), self.cleaned_types, set_columns=True)

            # Isolate the generic information
            n_snps, n_iid = self._chromosome_from_load(load_file)
            chromosome_values = {self.count_snp: n_snps, self.count_iid: n_iid,
                                 "Description": f"Chromosome {self.target_chromosome}"}
            config_dict[self.target_chromosome] = chromosome_values
            total_snps += n_snps

        print(f"Suggested LD_Radius based on {total_snps} / 3000 is {total_snps / 3000}")

        for key, value in config_dict.items():
            config_dict[key][self.herit] = self.herit_calculated * (config_dict[key][self.count_snp] / total_snps)
        config_dict["Genome"] = {f"{self.genome}_{self.herit}": self.herit_calculated}
        ArgMaker().write_yaml_group_dict(config_dict, self.working_dir, "genome_wide_config")

    def _chromosome_from_load(self, load_file):
        """
        We will need to extract the number of individuals, which requires re-parsing in the genetic file. We choose this
        based on the chromosome of the load file.
        """
        # Set the target from the first row and column of data
        self.target_chromosome = load_file.row_data[0][0]

        # Isolate the genetic load file according to summary, and use this to isolate number of individuals
        ref = self.construct_reference_panel()

        return load_file.column_length, ref.iid_count

    def calculate_genome_wide_heritability(self):
        """
        Once we have cleaned our data sets, we can store the genome wide data along side with some individual
        chromosome information that does not take the form of lists, in a yaml config file.
        """

        t0 = time.time()
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

        # Construct the write dict
        config_dict[self.genome] = {f"{self.avg_ld}": average_gw_ld_score, f"{self.herit}": gw_h2_ld_score_est,
                                    "Genome_Description": "Genome-wide Statistics"}

        # Log warnings if applicable then print information to terminal
        if gw_h2_ld_score_est > self.ld_radius / 10.0:
            print(Fore.RED + ec.ld_radius_to_large())
        if gw_h2_ld_score_est > 1:
            print(Fore.RED + ec.heritability_to_large())

        # Log information to the terminal
        self._genome_dict["Lambda Inflation"] = round(float(chi_square_lambda), 7)
        self._genome_dict["Mean LD Score"] = round(float(average_gw_ld_score), 7)
        self._genome_dict["Genome-Wide Heritability"] = round(float(gw_h2_ld_score_est), 7)
        mc.error_dict_to_terminal(self._genome_dict, "calculate_genome_wide_heritability", t0)

        # Construct config file
        ArgMaker().write_yaml_group_dict(config_dict, Path(self.working_dir, "PGS"), "genome_wide_config")

    def _heritability_by_chromosome(self, config_dict):
        """
        This will calculate the heritability for each chromosome cumulatively the values for the genome-wide
        calculation
        """
        cumulative_ld = sum_sq_beta = total_snps = 0
        for file in directory_iterator(self.ld_directory):

            # Infer the chromosome of this file via remove the 'LD' prefix
            chromosome = file[2:]
            print(f"Processing Chromosome {chromosome}")

            # Load the ld data from the ld_directory and use it to set sid/iid count
            data = load_pickle(self.ld_directory, file)
            sid_count, iid_count = data[self.norm_snps].shape

            # Isolate the betas from the Filtered to calculate the sum squared betas, also ld scores from ld data
            betas = self.sm_dict_from_csv(self.filter_directory, f"Filtered_{chromosome}.csv")[self.beta]
            sum_beta_sq = np.sum(np.array(betas) ** 2)
            ld_scores = data[self.ld_scores]

            # Calculate the heritability and average LD at a chromosome level
            heritability, average_ld = self._chromosome_heritability(chromosome, ld_scores, sum_beta_sq, sid_count)

            # Cumulate ld, snps, and sum square beta
            cumulative_ld += np.sum(ld_scores)
            total_snps += sid_count
            sum_sq_beta += np.sum(sum_beta_sq)

            # Store Values for config file
            chromosome_values = {self.herit: heritability, self.count_snp: sid_count, self.count_iid: iid_count,
                                 self.avg_ld: average_ld, "Description": f"Chromosome {chromosome}"}
            config_dict[chromosome] = chromosome_values
        return cumulative_ld, sum_sq_beta, total_snps

    def _chromosome_heritability(self, chromosome, ld_scores, sum_beta_sq, snp_count):
        """
        This will calculated the chromosome chi-squared lambda (maths from LDPred), and then take the maximum of 0.0001
        or the computed heritability of

        This chromosomes 1 - chi-sq lambda (or 1 if its less than 1 for reasons that are beyond me)
        ---------------------------------------------------------------------------------------
        Number of samples in the summary stats * (average ld score / number snps)

        :return: estimated heritability (h2 in ldpred)
        """
        if isinstance(self.herit_calculated, dict):
            try:
                return self.herit_calculated[chromosome]
            except KeyError:
                raise Exception("You have said you will provided pre-calculated heritability but failed to find it for"
                                f"chromosome {chromosome}")
        else:
            average_ld = np.mean(ld_scores)
            chi_sq_lambda = np.mean((self.sample_size * sum_beta_sq) / snp_count)

            herit = max(0.0001, (max(1.0, float(chi_sq_lambda)) - 1) / (self.sample_size * (average_ld / snp_count)))
            return herit, average_ld
