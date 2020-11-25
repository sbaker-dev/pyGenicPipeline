from pyGeneticPipe.geneticParsers.supportObjects import Variant, Nucleotide, SMVariant
from pyGeneticPipe.geneticParsers.plink.plinkObject import PlinkObject
from pyGeneticPipe.geneticParsers.bgen.bgenObject import BgenObject
from pyGeneticPipe.utils import error_codes as ec
from pyGeneticPipe.utils import misc as mc
from pyGeneticPipe.core.Input import Input
from pysnptools.distreader import Bgen
from pysnptools.snpreader import Bed
from collections import Counter
from operator import itemgetter
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
                            "Standard_Errors": 0, "Duplicate_Position": 0, "Ambiguous_SNP": 0, "Non_Matching": 0,
                            "Non_Allowed_Allele": 0}

    def clean_summary_statistics(self):

        # Note - this is basiclly becoming the -main- of prs, so will want to extract the chromosome bit so that it can
        # run in a multi-core manner

        # Check for input arguments
        self._assert_clean_summary_statistics()
        valid_chromosomes = self._validation_chromosomes()

        for chromosome in valid_chromosomes:
            print(f"Starting Chromosome: {chromosome}")

            # Load the validation and core samples, as well as the indexer
            load_path = str(self._select_file(chromosome))

            validation, core = self._construct_validation(load_path)

            # Clean the summary statistics
            sm_variants = self._clean_summary_stats(load_path, validation, core)
            print(f"\nCleaned summary statistics for chromosome: {chromosome}\n{self._error_dict}")

            # Isolate the ordered on bp_position common snps to extract the dosage information on core and validation
            self._isolate_dosage(validation, core, sm_variants, load_path)

            return

    def _clean_summary_stats(self, load_path, validation, core):
        """
        This will take the validation and core sample of snps, and check the snp against both sets. If the snp exists in
        the validation files, then it will go to cleaning the summary statistics for this chromosome line by line.
        """

        validation_snps, core_snps, indexer = self._load_variants(load_path, validation, core)

        sm_variants = []
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

                # If the snp exists in both the validation and core snp samples then clean this line, else skip.
                if (snp_id in validation_snps) and (snp_id in core_snps):
                    sm_variant = self._validate_summary_line(line, self._set_variant(snp_id, indexer))
                    if sm_variant:
                        sm_variants.append(sm_variant)

                else:
                    self._error_dict["Invalid_Snps"] += 1

        # Given we have only accepted snps that are within the validation / core, we should never have more snps in
        # summary than within the validation. If we do, something has gone critically wrong.
        assert len(sm_variants) <= len(validation_snps), ec.snp_overflow(len(sm_variants), len(validation_snps))
        assert len(sm_variants) <= len(core_snps), ec.snp_overflow(len(sm_variants), len(core_snps))
        return sm_variants

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

    def _construct_validation(self, load_path):
        """
        We need to construct a validation sample from the percentage the user provided and the iid_count, this then
        returns this slice of the sample from the start up to this percentage (Uses int so may be slightly above or
        below the percentage provided based on rounding / floating point errors), and then the rest of the sample of the
        core set

        :param load_path: Path to the relevant load file
        :return: The validation and core sample class holders
        """

        # todo Before spliting in to validation and core, allow a sample size modifier to remove people out (ie for ukb)
        # Set validation and core sets of sids based on the load type
        if self.load_type == ".bed":
            validation_size = self._set_validation_sample_size(Bed(load_path, count_A1=True).iid_count)
            validation = Bed(load_path, count_A1=True)[:validation_size, :]
            core = Bed(load_path, count_A1=True)[validation_size:, :]
            return validation, core

        elif self.load_type == ".bgen":
            # Bgen files store [variant id, rsid], we just want the rsid hence the [1]; see https://bit.ly/2J0C1kC
            validation_size = self._set_validation_sample_size(Bgen(load_path).iid_count)
            validation = Bgen(load_path)[:validation_size, :]
            core = Bgen(load_path)[validation_size:, :]
            return validation, core

        else:
            raise Exception("Unknown load type set")

    def _load_variants(self, load_path, validation, core):
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
            validation = validation.sid
            core = core.sid
            indexer = [PlinkObject(load_path).bim_index(), PlinkObject(load_path).bim_object()]

        elif self.load_type == ".bgen":
            # Bgen files store [variant id, rsid], we just want the rsid hence the [1]; see https://bit.ly/2J0C1kC
            validation = [snp.split(",")[1] for snp in validation.sid]
            core = [snp.split(",")[1] for snp in core.sid]
            indexer = [BgenObject(load_path).index_of_snps(),  BgenObject(load_path)]

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
        index_dict, indexer = indexer

        if self.load_type == ".bgen":
            return indexer.get_variant(index_dict[variant_id])
        else:
            return indexer.get_variant(index_dict[variant_id], True)

    def _validate_summary_line(self, line, variant):
        """
        This will validate a given line in the summary statistics for all the possible headers that it could contain
        against a validation

        """
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
            self._error_dict["Standard_Errors"] += 1
            return None

        # Construct a summary nucleotide to check against our genetic validation
        sm_nucleotide = Nucleotide(line[self.sm_effect_allele], line[self.sm_alt_allele])

        # Check to see if the nucleotide is within the ambiguous snp set
        if (variant.nucleotide() in self.ambiguous_snps) or (sm_nucleotide.to_tuple() in self.ambiguous_snps):
            self._error_dict["Ambiguous_SNP"] += 1
            return None

        # Check if the nucleotides are sane (Ie it is only a t c and g)
        if (variant.a1 not in self.allowed_alleles) or (variant.a2 not in self.allowed_alleles) or \
                (sm_nucleotide.a1 not in self.allowed_alleles) or (sm_nucleotide.a2 not in self.allowed_alleles):
            self._error_dict["Non_Allowed_Allele"] += 1
            return None

        # Check if nucleotides need to be flipped
        beta, beta_odds = self._flip_nucleotide(variant, sm_nucleotide, beta, beta_odds)
        if not beta or not beta_odds:
            self._error_dict["Non_Matching"] += 1
            return None

        # If we have the info, then add this
        info = -1
        if self.sm_info:
            info = float(line[self.sm_info])

        # If we have frequency, then extract it
        frequency = -1
        if self.frequencies:
            frequency = self._sum_stats_frequencies()

        # Return a object with all this information if all the checks pass
        return SMVariant(variant.chromosome, variant.variant_id, variant.bp_position, variant.a1, variant.a2, beta,
                         beta_odds, p_value, info, frequency)

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

    def _flip_nucleotide(self, variant, sm_nucleotide, beta, beta_odds):
        """
        Checks to see if the summary stats requires flipping by comparing it to the genetic variant from our sample and
        the flipped variant.

        If not flipping is required, returns beta and beta odds as is
        If flipping is required, returns flipped beta and beta odds
        If flipping is required but failed, returns None, None
        """

        # Construct the opposite strand for the genetic variant to see if the summary can match the flipped, vs implying
        # that flipping is required
        gen_flipped = Nucleotide(self.allele_flip[variant.a1], self.allele_flip[variant.a2])

        # If either condition is not meet then we need to try to flip the nucleotide, else return beta and beta odds
        if not (np.all(variant.nucleotide(True) == sm_nucleotide.to_list())) or \
                (np.all(gen_flipped.to_list() == sm_nucleotide.to_list())):

            flip_nts = (variant.a2 == sm_nucleotide.a1 and variant.a1 == sm_nucleotide.a2) or (
                    gen_flipped.a2 == sm_nucleotide.a1 and gen_flipped.a1 == sm_nucleotide.a2)

            # If flip successful then invert beta and beta_odds, else return none
            if flip_nts:
                return -beta, -beta_odds
            else:
                return None, None
        else:
            return beta, beta_odds

    def _common_ordered_snps(self, sm_variants):
        """
        Because we require a snp to be within the validation and core samples in order to be accepted, we know that the
        size of sm_variants will be less than or equal to the number in the genetic samples. As such the common snps,
        are the snps in the summary statistics after cleaning.

        Given we will be using these snps for dosage, and dosage extraction requires snps to be in a set format, here
        we construct a list of snps and order them by their base pair position in a format of:

        .bed: ["rs123", "rs124", ... "rsN"]
        .bgen: ["rs123,rs123", "rs124,rs124", ... "rsN,rsN"]

        :param sm_variants: The variants found in the summary statistics cleaning stage for this chromosome
        :return: snps order by base pair position in the format of the load type.
        """

        if self.load_type == ".bed":
            variant_by_bp = [[variant.variant_id, variant.bp_position] for variant in sm_variants]
        elif self.load_type == ".bgen":
            variant_by_bp = [[variant.bgen_variant_id(), variant.bp_position] for variant in sm_variants]
        else:
            raise Exception(f"Critical Error: Unknown load type {self.load_type} found in _common_snps")

        return [rs_id for rs_id, bp in sorted(variant_by_bp, key=itemgetter(1))]

    def _isolate_dosage(self, validation, core, sm_variants, load_path):

        print(validation)
        print(core)
        print(len(sm_variants))

        # todo Currently we are in the validation stage, so we still need load path but this should be removed when
        #  ready

        if self.load_type == ".bed":
            validation = Bed(load_path, count_A1=True)
            ordered_common = validation[:, validation.sid_to_index(self._common_ordered_snps(sm_variants))].read().val

        elif self.load_type == ".bgen":
            print("Bgen load type, so need to restructure return type ... will take longer!")
            validation = Bgen(load_path)
            ordered_common = validation[:, validation.sid_to_index(self._common_ordered_snps(sm_variants))].read().val
        else:
            raise Exception(f"Critical Error: Unknown load type {self.load_type} found in _isolate_dosage")

        # LDPred used a system of snps*ids but pysnptools uses ids*snps so we need to invert them
        dosage = self.flip_list(ordered_common)

        print(dosage[0])

        frequency = np.sum(dosage, 1, dtype='float32') / (2 * float(len(dosage[0])))

        print(frequency)

        # todo This is very slow, and looks like it has the potential to be very memeory intensive try access a single
        #  snp at a time and construct a dosage that way, and then we can just hold dosage data for a single snp across
        #  individuals at a time. May also need to better configure flip_lists to take arguments in terms of config

    def flip_list(self, list_of_lists):
        """
        This will take a list of lists, and then flip it. It requires all sub lists to be the same length.

        NOTE: This is very heavy in performance and could use a speed up
        """

        list_of_keys = Counter([len(sub_list) for sub_list in list_of_lists])
        sub_key_length = list(list_of_keys.keys())

        value_dict = {0: 2, 1: 1,  2: 0}

        if len(list_of_keys.keys()) == 1:
            if self.load_type == ".bgen":
                # Load the dosages for each individual, but it gets loaded as [1, 0, 0] OR [0, 1, 0] OR [0, 0, 1] unlike
                # the plink file so we have to reorder this for consistency.
                return [np.array([np.where(sub[i] == 1)[0][0] for sub in list_of_lists])
                        for i in range(sub_key_length[0])]

            elif self.load_type == ".bed":
                # returns the operssite of ldpred so we need to remap 0 - 2 and 2 - 0
                return [np.array([value_dict[int(sub[i])] for sub in list_of_lists], dtype=np.int8)
                        for i in range(sub_key_length[0])]
            else:
                raise Exception(f"Critical Error: Unknown load type {self.load_type} found in _isolate_dosage")

        else:
            raise Exception(f"Sub lists should be all of the same length yet found lengths {sub_key_length}")

    def _sum_stats_frequencies(self):
        raise NotImplementedError("Frequencies are not yet implemented")
