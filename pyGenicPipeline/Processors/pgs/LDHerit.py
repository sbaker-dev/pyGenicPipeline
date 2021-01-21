from pyGenicPipeline.utils.ArgMaker import ArgMaker
from pyGenicPipeline.utils import errors as ec
from pyGenicPipeline.utils import misc as mc
from pyGenicPipeline.core.Input import Input

from miscSupports import directory_iterator, terminal_time
from csvObject import CsvObject
from colorama import Fore
from pathlib import Path
import numpy as np
import pickle
import time


class LDHerit(Input):
    def __init__(self, args):
        super().__init__(args)

        self._genome_dict = {"Genome-wide Stats": "Values", "Lambda Inflation": 0.0, "Mean LD Score": 0.0,
                             "Genome-Wide Heritability": 0.0}

    def calculate_ld(self):
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

        # Calculate the Ld scores, ld reference dict, and the ld matrix
        ld_dict, ld_scores = self._calculate_ld_scores(normalised_snps, sid_count, iid_count)
        ld_matrix = self._calculate_ld_matrix(normalised_snps, sid_count, iid_count)

        # Write the information as a combined dict to pickle file
        write_ld = {"LD_Dict": ld_dict, "LD_Scores": ld_scores, "LD_Matrix": ld_matrix, "Norm_Snps": normalised_snps,
                    "Snp_Std": std}
        with open(f"{self.ld_directory}/LD{self.target_chromosome}", "wb") as f:
            pickle.dump(write_ld, f)

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

    def _calculate_ld_matrix(self, normalised_snps, sid_count, iid_count):
        """This will use a window of snps rather than snps in local LD to calculate the disequilibrium"""
        ld_matrix = {}
        for index, window in enumerate(range(0, sid_count, self.ld_radius * 2)):
            # Isolate the snps in the current window
            snp_window = self.window_values(normalised_snps, sid_count, window)

            # Calculate the disequilibrium, then append the shrunk r2 matrix to the ld_matrix
            disequilibrium_dp = np.dot(snp_window, snp_window.T) / iid_count
            ld_matrix[index] = mc.shrink_r2_matrix(disequilibrium_dp, iid_count)

        print(f"Constructed LD matrix {terminal_time()}")
        return ld_matrix

    def calculate_genome_wide_heritability(self):
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

    def distribute_heritability_genome_wide(self):
        """If we can't calculate heritability, distribute it from a provided float"""
        total_snps = 0
        config_dict = {}
        for file in directory_iterator(self.clean_directory):
            print(file)
            load_file = CsvObject(Path(self.clean_directory, file), self.cleaned_types[:-1], set_columns=True)

            # Isolate the generic information
            chromosome, n_snps, n_iid = self._chromosome_from_load(load_file)
            chromosome_values = {self.count_snp: n_snps, self.count_iid: n_iid,
                                 "Description": f"Chromosome {chromosome}"}
            config_dict[chromosome] = chromosome_values
            total_snps += n_snps

        print(f"Suggested LD_Radius based on {total_snps} / 3000 is {total_snps / 3000}")

        for key, value in config_dict.items():
            config_dict[key][self.herit] = self.herit_calculated * (config_dict[key][self.count_snp] / total_snps)

        config_dict["Genome"] = {f"{self.genome_key}_{self.herit}": self.herit_calculated}
        ArgMaker().write_yaml_config_dict(config_dict, self.working_dir, "genome_wide_config")

    def _heritability_by_chromosome(self, config_dict):
        """
        This will calculate the heritability for each chromosome cumulatively the values for the genome-wide
        calculation
        """
        cumulative_ld = sum_sq_beta = total_snps = 0
        for file in directory_iterator(self.clean_directory):
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


    def _chromosome_heritability(self, load_file, chromosome, snp_count):
        """
        This will calculated the chromosome chi-squared lambda (maths from LDPred), and then take the maximum of 0.0001
        or the computed heritability of

        This chromosomes 1 - chi-sq lambda (or 1 if its less than 1 for reasons that are beyond me)
        ---------------------------------------------------------------------------------------
        Number of samples in the summary stats * (average ld score / number snps)

        :return: estimated heritability (h2 in ldpred)
        """
        if self.herit_calculated:
            try:
                return self.herit_calculated[chromosome]
            except KeyError:
                raise Exception("You have said you will provided pre-calculated heritability but failed to find it for"
                                f"chromosome {chromosome}")
        else:
            average_ld = np.mean(load_file.column_data[self.c_ld_score])
            sum_beta_sq = np.sum(np.array(load_file.column_data[self.c_beta]) ** 2)
            chi_sq_lambda = np.mean((self.sample_size * sum_beta_sq) / snp_count)

            herit = max(0.0001, (max(1.0, float(chi_sq_lambda)) - 1) / (self.sample_size * (average_ld / snp_count)))
            return herit, average_ld
