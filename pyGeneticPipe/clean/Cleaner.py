from pyGeneticPipe.geneticParsers.plink.plinkObject import PlinkObject
from pyGeneticPipe.geneticParsers.bgen.bgenObject import BgenObject
from pyGeneticPipe.geneticParsers.supportObjects import Variant
from pyGeneticPipe.utils import error_codes as ec
from pyGeneticPipe.utils import misc as mc
from pyGeneticPipe.core.Input import Input
from pysnptools.distreader import Bgen
from pysnptools.snpreader import Bed
from pathlib import Path
from scipy import stats
import numpy as np
import pickle
import gzip
import re


class Cleaner(Input):
    def __init__(self, args):
        super().__init__(args)

        self._error_dict = {"Invalid_Snps": 0, "Chromosome": 0, "Position": 0, "Effect_Size": 0, "P_Value": 0,
                            "Standard_Errors": 0, "Duplicate_Position": 0}

    def clean_summary_statistics(self):

        # Check for input arguments
        self._assert_clean_summary_statistics()

        valid_chromosomes = self._validation_chromosomes()

        for chromosome in valid_chromosomes:
            print(chromosome)

            load_path = str(self._select_file(chromosome))
            validation, core, indexer = self._load_variants(load_path)

            with mc.open_setter(self.summary_file)(self.summary_file) as file:
                # Skip header row
                file.readline()

                # For each line in the GWAS Summary file
                for index, line in enumerate(file):
                    if index % 10000 == 0 and self.debug:
                        print(f"{index}")

                    # Decode the line and extract the snp_id
                    line = mc.decode_line(line, self.zipped)
                    snp_id = line[self.sm_snp_id]
                    if (snp_id in validation) and (snp_id in core):
                        self._validate_summary_line(line, self._set_variant(snp_id, indexer))
                    else:
                        self._error_dict["Invalid_Snps"] += 1

                    file.close()
                    break

            return

    def _select_file(self, chromosome):
        """
        For a given chromosome, get the respective file
        :param chromosome: Current chromosome to be loaded
        :return: Path to the current file as a Path from pathlib
        """
        for file in mc.directory_iterator(self.load_directory):
            if Path(self.load_directory, file).suffix == self.load_type:
                if int(re.sub(r'[\D]', "", Path(self.load_directory, file).stem)) == chromosome:
                    return Path(self.load_directory, file)

        raise Exception(f"Failed to find any relevant file for {chromosome} in {self.load_directory}")

    def _validation_chromosomes(self):
        """
        This will create a dataset of all the chromosomes that we have to work with our validation group in the
        h5py file

        :return: A list of valid chromosomes
        :rtype: list
        """

        valid_chromosomes = []
        for file in mc.directory_iterator(self.load_directory):
            if Path(self.load_directory, file).suffix == self.load_type:
                valid_chromosomes.append(int(re.sub(r'[\D]', "", Path(self.load_directory, file).stem)))
        valid_chromosomes.sort()
        return valid_chromosomes

    def _load_variants(self, load_path):
        """
        Load variants, for .bgen or plink files, as a set of snps that exist within the current chromosome. Uses the
        validation percentage to construct a validation group, and returns the set of snps for each group. If hap_map_3
        is enabled, it will strip out snps not in hap_map_3.

        We will also need a way to index out the variant information, so we set the indexer according to the load type

        :param load_path: Current Chromosome file
        :return: Set of the validation and core set
        """
        hap_map_3 = self._load_hap_map_3()

        #  Set validation and core sets of sids based on the load type
        if self.load_type == ".bed":
            validation_size = self._set_validation_sample_size(Bed(load_path, count_A1=True).iid_count)
            validation = Bed(load_path, count_A1=True)[:validation_size, :].sid
            core = Bed(load_path, count_A1=True)[validation_size:, :].sid
            indexer = [PlinkObject(load_path).bim_index(), PlinkObject(load_path).bim_object()]

        elif self.load_type == ".bgen":
            # Bgen files store [variant id, rsid], we just want the rsid hence the [1]; see https://bit.ly/2J0C1kC
            validation_size = self._set_validation_sample_size(Bgen(load_path).iid_count)
            validation = [snp.split(",")[1] for snp in Bgen(load_path)[:validation_size, :].sid]
            core = [snp.split(",")[1] for snp in Bgen(load_path)[validation_size:, :].sid]
            indexer = BgenObject(load_path)

        else:
            raise Exception("Unknown load type set")

        # If we only want the hap_map_3 snps then check each snp against the set of hap_map_3
        if hap_map_3:
            validation = [snp for snp in validation if snp in hap_map_3]
            core = [snp for snp in core if snp in hap_map_3]

        return set(validation), set(core), indexer

    def _load_hap_map_3(self):
        """
        Users may wish to limit valid snps to those found within HapMap3. If they do, they need to provide a path to the
        hapmap3 snp file which will be check that it exists, have the snps extracted and return. Otherwise set to none

        :return: The set of the valid HapMap3 snps or None
        :rtype: set | None
        """

        if self.hap_map_3_file:
            print("WARNING: THIS IS UNTESTED CODE FROM LDPRED")
            # Construct path as an object and check it exists

            # If the HapMap3 file exists, then extract the snp ids and return them
            f = gzip.open(self.hap_map_3_file, 'r')
            hm3_sids = pickle.load(f)
            f.close()
            return hm3_sids
        else:
            return None

    def _assert_clean_summary_statistics(self):
        """
        clean_summary_statistics requires

        The project file, for writing too
        That the cleaning has not already been undertaken
        That the summary file path exists
        That the load type for the genetic data exists
        That the load directory containing the chromosome split data exists
        That the Validation_Size has been set
        """
        # Check parameters, validate that validation has been run, and that clean summary has not.
        assert self.project_file, ec.missing_arg(self.operation, "Project_Name")
        assert self.h5_summary not in self.project_file.keys(), ec.appending_error(self.project_name, self.h5_summary)
        assert self.summary_file, ec.missing_arg(self.operation, "Summary_Path")
        assert self.load_type, ec.missing_arg(self.operation, "Load_Type")
        assert self.load_directory, ec.missing_arg(self.operation, "Load_Directory")
        assert self.validation_size, ec.missing_arg(self.operation, "Validation_Size")

    def _set_validation_sample_size(self, full_sample_size):
        """
        This will return the value of the validation in terms of individuals rather than a percentage that the user
        specified

        :param full_sample_size: An integer of the number of samples in the full sample
        :type full_sample_size: int

        :return: Integer of the number of samples required in the validation sample
        :rtype: int
        """
        return int(full_sample_size * self.validation_size)

    def _set_variant(self, variant_id, indexer):
        """
        Get the variant id from genetic information based on load type
        :param variant_id: Current snp name to extract
        :param indexer: Indexer to extract from

        :return: Variant
        :rtype: Variant
        """
        if self.load_type == ".bgen":
            return indexer.get_variant(variant_id)
        else:
            index_dict, indexer = indexer
            return indexer.get_variant(index_dict[variant_id], True)

    def _validate_summary_line(self, line, variant):

        # If we have chromosomes in our summary statistics check the chromosome of the snp against the validation
        if self.sm_chromosome and line[self.sm_chromosome] != variant.chromosome:
            self._error_dict["Chromosome"] += 1
            return None

        # If we have base pair position in our summary then validate the base pair
        if self.sm_bp_position and int(line[self.sm_bp_position]) != variant.bp_position:
            self._error_dict["Position"] += 1
            return None

        # Set beta unless it is not a finite number
        beta = float(line[self.sm_effect_size])
        if not np.isfinite(beta):
            self._error_dict["Effect_Size"] += 1
            return None

        # Set p value as long as it is not zero or a non finite number
        p_value = float(line[self.sm_p_value])
        if not np.isfinite(p_value) or p_value == 0:
            self._error_dict["P_Value"] += 1
            return None

        # calculate the beta and beta odds
        beta, beta_odds = self._calculate_beta(beta, line, p_value)
        if not beta or not beta_odds:
            return None

    def _calculate_beta(self, beta, line, p_value):
        """
        Calculate both the beta, and the beta odds depending on the effect_type and if the user wants to constructed a
        standardised z score or not

        :param beta: beta from summary stats
        :type beta: float

        :param line: The line of the current file
        :type line: list

        :param p_value: p value from summary stats
        :type p_value: float

        :return: If successful, return the beta and beta odds otherwise None and None
        """

        beta_odds = self._beta_by_type(beta)

        # If effect type is Best linear unbiased prediction (BLUP) return beta
        if self.effect_type == "BLUP":
            return beta, beta_odds

        # If we want to compute z scores, compute them as long as standard errors are valid
        elif self.z_scores:
            se = float(line[self.sm_standard_errors])
            if not np.isfinite(se) or se == 0:
                self._error_dict["Standard_Errors"] += 1
                return None, None
            else:
                return self._beta_from_se(beta, beta_odds, se), beta_odds

        # Otherwise compute beta from p values
        else:
            return np.sign(beta) * (stats.norm.ppf(p_value / 2.0) / np.sqrt(self.sample_size)), beta_odds

    def _beta_from_se(self, beta, beta_typed, se):
        """Calculate z score with standard error"""
        if self.effect_type == "OR":
            abs_beta = np.absolute(1 - beta) / se
        else:
            abs_beta = np.absolute(beta) / se
        return np.sign(beta_typed) * (abs_beta / np.sqrt(self.sample_size))

    def _beta_by_type(self, beta):
        """
        If we are working with ods rations we need to take the log of the read beta
        """
        if self.effect_type == "OR":
            return np.log(beta)
        else:
            return beta
