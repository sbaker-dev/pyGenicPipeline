from pyGeneticPipe.geneticParsers.supportObjects import Variant, Nucleotide
from pyGeneticPipe.geneticParsers.plink.plinkObject import PlinkObject
from pyGeneticPipe.geneticParsers.bgen.bgenObject import BgenObject
from pyGeneticPipe.utils import error_codes as ec
from pyGeneticPipe.utils import misc as mc
from pyGeneticPipe.core.Input import Input
from pysnptools.distreader import Bgen
from pysnptools.snpreader import Bed
from colorama import Fore
from pathlib import Path
from scipy import stats
import numpy as np
import pickle
import gzip
import re


class Cleaner(Input):
    def __init__(self, args):
        super().__init__(args)

        # Depending on how this ends up being constructed we might want to log this out as a .txt rather than just
        # print this to a terminal.
        self._error_dict = {"Removal case": "Count", "Invalid_Snps": 0, self.chromosome: 0, self.bp_position: 0,
                            self.effect_size: 0, self.p_value: 0, self.standard_errors: 0, "Ambiguous_SNP": 0,
                            "Non_Allowed_Allele": 0, "Non_Matching": 0, "Filtered_Frequency": 0,
                            "Filtered_MAF": 0, "Monomorphic": 0, "Accepted_Snps": 0}
        self._summary_last_position = 0

    def clean_summary_statistics(self):

        # Note - this is basically becoming the -main- of prs, so will want to extract the chromosome bit so that it can
        # run in a multi-core manner
        # Make sure output is on a per chromosome level in a project directory so that the next step can read in all the
        #  data when it is required to be genome wide

        # Check for input arguments
        self._assert_clean_summary_statistics()
        valid_chromosomes = self._validation_chromosomes()

        for chromosome in valid_chromosomes:
            print(f"Starting Chromosome: {chromosome}")

            # Load the validation and core samples, as well as the indexer
            load_path = str(self._select_file(chromosome))
            validation, core = self._construct_validation(load_path)

            # Clean the summary statistics
            sm_variants = self._clean_summary_stats(load_path, validation, core, chromosome)
            if not sm_variants:
                print(f"No variants found for {chromosome}")
                return

            # Filter the summary stats
            self._filter_snps(load_path, sm_variants)

            # Log to terminal what has been filtered / removed
            self._error_dict_to_terminal(chromosome)

    def _clean_summary_stats(self, load_path, validation, core, chromosome):
        """
        This will take the validation and core sample of snps, and check the snp against both sets. If the snp exists in
        the validation files, then it will go to cleaning the summary statistics for this chromosome line by line.
        """

        validation_snps, core_snps, indexer = self._load_variants(load_path, validation, core)

        sm_variants = []
        sm_line = []
        with mc.open_setter(self.summary_file)(self.summary_file) as file:
            self._seek_to_start(chromosome, file)

            # For each line in the GWAS Summary file
            for index, line_byte in enumerate(file):
                if index % 10000 == 0 and self.debug:
                    print(f"{index}")

                # Decode the line and extract the snp_id
                line = mc.decode_line(line_byte, self.zipped)
                snp_id = line[self.sm_snp_id]

                # If the snp exists in both the validation and core snp samples then clean this line, else skip.
                if (snp_id in validation_snps) and (snp_id in core_snps):
                    sm_variants.append(self._set_variant(snp_id, indexer))
                    sm_line.append(line)

                else:
                    # If the chromosomes exist in summary statistics we can terminate this for loop when we are no
                    # longer in the right zone and set tell to seek to this position for the next chromosome
                    if (self.sm_chromosome is not None) and (int(line[self.sm_chromosome]) > chromosome):
                        self._summary_last_position = file.tell() - len(line_byte)
                        file.close()
                        break
                    else:
                        self._error_dict["Invalid_Snps"] += 1

        sm_dict = {self.sm_lines: np.array(sm_line), self.sm_variants: np.array(sm_variants)}
        sm_dict = self._validate_summary_lines(sm_dict)
        if not sm_dict:
            return None

        # Construct the order from the base pair position
        order = np.argsort(np.array([variant.bp_position for variant in sm_dict[self.sm_variants]]))

        # Given we have only accepted snps that are within the validation / core, we should never have more snps in
        # summary than within the validation. If we do, something has gone critically wrong.
        assert len(order) <= len(validation_snps), ec.snp_overflow(len(order), len(validation_snps))
        assert len(order) <= len(core_snps), ec.snp_overflow(len(order), len(core_snps))
        if self._filter_array(sm_dict, order) == 1:
            return None
        else:
            return sm_dict

    def _validate_summary_lines(self, sm_dict):
        """This will load in each possible header, and clean our dict of values by filtering"""

        # If we have chromosomes in our summary statistics check the chromosome of the snps against the validation
        if self.sm_chromosome is not None:
            if self._validation_equality(self.sm_chromosome, self.chromosome, sm_dict) == 1:
                return None

        # If we have base pair position in our summary then validate the base pair
        if self.bp_position is not None:
            if self._validation_equality(self.sm_bp_position, self.bp_position, sm_dict, int) == 1:
                return None

        # Clean the summary stats effect sizes for calculation of beta later
        if self._validation_finite(sm_dict, self.sm_effect_size, self.effect_size) == 1:
            return None

        # Clean the P values
        if self._validation_finite(sm_dict, self.sm_p_value, self.p_value) == 1:
            return None

        # If we are using z scores we need to load and clean the standard errors column
        if self.z_scores:
            if self._validation_finite(sm_dict, self.sm_standard_errors, self.standard_errors) == 1:
                return None

        # Use the raw beta, standard errors, and p value if required to construct beta and beta_odds
        self._validation_betas(sm_dict)

        # Check that the nucleotides are sane and flip them if required
        if self._validate_nucleotides(sm_dict) == 1:
            return None

        # Calculate the frequencies and set info if it exists
        sm_dict[self.frequency] = np.array([self._sum_stats_frequencies(line) for line in sm_dict[self.sm_lines]])
        sm_dict[self.info] = self._validate_info(sm_dict[self.sm_lines])

        return sm_dict

    @staticmethod
    def _line_array(line_index, line_array, type_np=None):
        """Construct an array of a single line index form line_array considering the type"""
        if type_np:
            return np.array([line[line_index] for line in line_array], dtype=type_np)
        else:
            return np.array([line[line_index] for line in line_array])

    @staticmethod
    def _variant_array(variant_key, variant_array):
        """Construct an array of a single variant_array's item using getitem via variant key"""
        return np.array([variant[variant_key] for variant in variant_array])

    @staticmethod
    def _filter_array(dict_to_filter, array_filter):
        """
        Filter out anything that is no longer required If the length of the array becomes zero, pass an error code of
        1. Otherwise return 0.
        """
        for key, value in zip(dict_to_filter.keys(), dict_to_filter.values()):
            dict_to_filter[key] = value[array_filter]

        if np.array([len(value) for value in dict_to_filter.values()])[0] == 0:
            return 1
        else:
            return 0

    def _validate_info(self, sm_line):
        """Construct infos if they exist in the summary stats else return an array of length of summary dict"""
        if self.sm_info is not None:
            infos = self._line_array(self.sm_info, sm_line, float)
        else:
            infos = np.empty(len(sm_line))
            infos.fill(-1)
        return infos

    def _validation_equality(self, line_index, variant_key, summary_dict, line_type=None):
        """
        Not all summary statistics may have chromosome or bp indexes and in this case information can be returned from
        the genetic variant. However if the information does exist, then we cross check to make sure it is equal in the
        summary and genetic files. If it is not, we filter out this snp. This is the generalised method which can also
        be used for other equality

        :param line_index: The index for the summary line to construct an array from the
        :type line_index: int

        :param line_type: The type of the summary line to be return as, defaults to none which will return a string
        :type line_type: None | type

        :param variant_key: Key to access Variant via getitem and set error dict
        :type variant_key: str

        :param summary_dict: The summary dictionary to hold information so that we can filter it
        :type summary_dict: dict

        :return: Nothing, construct and use the filter on summary_dict then stop
        """
        # Construct an array of summary and genetic chromosomes
        summary_array = self._line_array(line_index, summary_dict[self.sm_lines], line_type)
        variant_array = self._variant_array(variant_key.lower(), summary_dict[self.sm_variants])

        # Filter of True if the variant and summary match, else False which will remove this snp
        obj_filter = summary_array == variant_array
        self._error_dict[variant_key] = len(obj_filter) - np.sum(obj_filter)
        return self._filter_array(summary_dict, obj_filter)

    def _validation_finite(self, summary_dict, line_index, summary_key):
        """
        Numeric columns need to screened for values being finite and not equal to zero. Unlike _validation_equality this
        method also appended the information to summary_dict as it has created new information not within the genetic
        variants rather than just screening pre-existing information

        :param line_index: The index for the summary line to construct an array from the
        :type line_index: int

        :param summary_key: A string key that is used for accessing this attribute
        :type summary_key: str

        :param summary_dict: The summary dictionary to hold information so that we can filter it
        :type summary_dict: dict

        :return: Nothing, construct the filter and then filter all attributes within the summary dict
        """
        # Construct an array for this numeric summary_key and add it to summary dict under the name of summary_key
        summary_dict[summary_key] = self._line_array(line_index, summary_dict[self.sm_lines], float)

        # Filter out anything that is not finite or is equal to zero
        obj_filter = np.array([True if np.isfinite(obj) and obj != 0 else False for obj in summary_dict[summary_key]])
        self._error_dict[summary_key] = len(obj_filter) - np.sum(obj_filter)
        return self._filter_array(summary_dict, obj_filter)

    def _validation_betas(self, sm_dict):
        """
        Calculate both the beta, and the beta odds depending on the effect_type and if the user wants to constructed a
        standardised z score or not

        :param sm_dict: Summary dict of values to read from and write too
        :type sm_dict: dict

        :return: Nothing, append values to dicts when constructed
        """

        # Construct the log_odds based on the effect_type
        sm_dict[self.log_odds] = self._beta_by_type(sm_dict)

        # If effect type is Best linear unbiased prediction (BLUP) return effect size column as beta
        if self.effect_type == "BLUP":
            sm_dict[self.beta] = sm_dict[self.effect_size].copy()

        # If we want to compute z scores, compute them as long as standard errors are valid
        elif self.z_scores:
            # The betas need to be altered for z scores
            if self.effect_type == "OR":
                abs_beta = np.array([(np.absolute(1 - beta) / se)
                                     for beta, se in zip(sm_dict[self.effect_size], sm_dict[self.standard_errors])])
            else:
                abs_beta = np.array([(np.absolute(beta) / se)
                                     for beta, se in zip(sm_dict[self.effect_size], sm_dict[self.standard_errors])])

            sm_dict[self.beta] = np.array([np.sign(beta_t) * (ab / np.sqrt(self.sample_size))
                                           for beta_t, ab in zip(sm_dict[self.log_odds], abs_beta)])

        # Otherwise compute beta from p values
        else:
            # probability density function
            pdf = stats.norm.ppf(sm_dict[self.p_value] / 2.0)
            sm_dict[self.beta] = np.array([np.sign(beta) * (pdf / np.sqrt(self.sample_size))
                                           for beta, pdf in zip(sm_dict[self.effect_size], pdf)])

    def _validate_nucleotides(self, sm_dict):
        """
        This wil validate the nucleotides against ambiguous snps, invalid snps, and flip the snps if possible whilst
        removing them if flipping fails.

        :param sm_dict: The summary dictionary to hold information so that we can filter it
        :type sm_dict: dict

        :return: Nothing, filter the arrays if required and flip betas if the allele is flipped
        """

        # Ambiguous
        # Construct the summary nucleotide
        effected_allele = self._line_array(self.sm_effect_allele, sm_dict[self.sm_lines])
        alt_allele = self._line_array(self.sm_alt_allele, sm_dict[self.sm_lines])
        sm_dict[self.nucleotide] = np.array([Nucleotide(e, a) for e, a in zip(effected_allele, alt_allele)])

        # Filter out any snps where the summery or variant Nucleotide is ambiguous
        filter_ambiguous = [False if (sm_nuc.to_tuple() in self.ambiguous_snps) or
                                     ((var_nuc.a1, var_nuc.a2) in self.ambiguous_snps)
                            else True
                            for sm_nuc, var_nuc in zip(sm_dict[self.nucleotide], sm_dict[self.sm_variants])]
        self._error_dict["Ambiguous_SNP"] = len(filter_ambiguous) - np.sum(filter_ambiguous)
        if self._filter_array(sm_dict, filter_ambiguous) == 1:
            return 1

        # Sainity Check
        # Filter out any snps that do not pass a sanity check (Only a t c and g)
        allowed_filter = [False if (sm_nuc.a1 not in self.allowed_alleles) or
                                   (sm_nuc.a2 not in self.allowed_alleles) or
                                   (var_nuc.a1 not in self.allowed_alleles) or
                                   (var_nuc.a2 not in self.allowed_alleles)
                          else True
                          for sm_nuc, var_nuc in zip(sm_dict[self.nucleotide], sm_dict[self.sm_variants])]
        self._error_dict["Non_Allowed_Allele"] = len(allowed_filter) - np.sum(allowed_filter)
        if self._filter_array(sm_dict, allowed_filter) == 1:
            return 1

        # Determine Flipping
        # Construct a flip status of 1, 0, -1 for No flipping, failed flipping, and flipped successfully which we can
        # multiple our betas by
        sm_dict["Flip"] = np.array([self._flip_nucleotide(var_nuc, sm_nuc) for var_nuc, sm_nuc in
                                    zip(sm_dict[self.sm_variants], sm_dict[self.nucleotide])])
        filter_flipped = np.array([False if flipped == 0 else True for flipped in sm_dict["Flip"]])
        self._error_dict["Non_Matching"] = len(filter_flipped) - np.sum(filter_flipped)
        if self._filter_array(sm_dict, filter_flipped) == 1:
            return 1

        # Now we have filtered away any errors, multiple the dicts beta and log_odds elements by 1 or -1 based on no
        # flipping or requiring flipping
        sm_dict[self.beta] = sm_dict[self.beta] * sm_dict["Flip"]
        sm_dict[self.log_odds] = sm_dict[self.log_odds] * sm_dict["Flip"]

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
            indexer = [BgenObject(load_path).index_of_snps(), BgenObject(load_path)]

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
            # If the HapMap3 file exists, then extract the snp ids and return them
            f = gzip.open(self.hap_map_3_file, 'r')
            hm3_sids = pickle.load(f)
            f.close()
            return set(hm3_sids)
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

    def _beta_by_type(self, sm_dict):
        """
        If we are working with Odds ratios we need to take the log of the read beta
        """
        if self.effect_type == "OR":
            return np.array([np.log(beta) for beta in sm_dict[self.effect_size]])
        else:
            return sm_dict[self.effect_size].copy()

    def _flip_nucleotide(self, variant, sm_nucleotide):
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
                return -1
            else:
                return 0
        else:
            return 1

    def _isolate_raw_snps(self, gen_file, sm_dict):
        """
        This will isolate the raw snps for a given bed or bgen file

        :param gen_file: Genetic file you wish to load from
        :param sm_dict: dict of clean information
        :return: raw snps
        """
        ordered_common = gen_file[:, gen_file.sid_to_index(self._extract_variant_name(sm_dict))].read().val

        # bed returns 2, 1, 0 rather than 0, 1, 2 although it says its 0, 1, 2; so this inverts it
        if self.load_type == ".bed":
            return np.array([abs(snp - 2) for snp in ordered_common.T])

        # We have a [1, 0, 0], [0, 1, 0], [0, 0, 1] array return for 0, 1, 2 respectively. So if we multiple the arrays
        # by their index position and then sum them we get [0, 1, 2]
        elif self.load_type == ".bgen":
            return sum(np.array([snp * i for i, snp in enumerate(ordered_common.T)]))

        else:
            raise Exception(f"Critical Error: Unknown load type {self.load_type} found in _isolate_dosage")

    def _extract_variant_name(self, sm_dict):
        """
        Different file types have different naming standards.

        .bed: ["rs123", "rs124", ... "rsN"]
        .bgen: ["rs123,rs123", "rs124,rs124", ... "rsN,rsN"]

        This will standardise the names to be a list of type equivalent to bed
        :param sm_dict: dict of clean information
        :return: list of snp names
        """
        if self.load_type == ".bed":
            return [variant.variant_id for variant in sm_dict[self.sm_variants]]
        elif self.load_type == ".bgen":
            print("Bgen load type, so need to restructure return type ... will take a bit longer longer!")
            return [variant.bgen_variant_id() for variant in sm_dict[self.sm_variants]]
        else:
            raise Exception(f"Critical Error: Unknown load type {self.load_type} found in _isolate_dosage")

    def _filter_snps(self, load_path, sm_dict):
        # testing validation for comparision TEMP
        validation = Bgen(load_path)

        # Construct the genetic raw snps and genetric freqs
        sm_dict["Validation_Raw_Snps"] = self._isolate_raw_snps(validation, sm_dict)
        sm_dict["Freqs"] = np.sum(sm_dict["Validation_Raw_Snps"], 1, dtype='float32') / (2 * float(validation.iid_count))

        # If the frequencies in the summary stats are not just a list of -1 errors then screen genetic snp frequencies
        if (self.freq_discrepancy < 1) and (np.sum(sm_dict[self.frequency] == -1) != len(sm_dict[self.frequency])):
            freq_filter = np.logical_or(np.absolute(sm_dict[self.frequency] - sm_dict["Freqs"]) < self.freq_discrepancy,
                                        np.absolute(sm_dict[self.frequency] + (sm_dict["Freqs"] - 1)) < self.freq_discrepancy)
            freq_filter = np.logical_or(freq_filter, sm_dict[self.frequency] <= 0)
            self._error_dict["Filtered_Frequency"] = len(freq_filter) - np.sum(freq_filter)
            if self._filter_array(sm_dict, freq_filter) == 0:
                return None

        if self.maf_min > 0:
            maf_filter = (sm_dict["Freqs"] > self.maf_min) * (sm_dict["Freqs"] < (1 - self.maf_min))
            self._error_dict["Filtered_MAF"] = len(maf_filter) - np.sum(maf_filter)
            if self._filter_array(sm_dict, maf_filter) == 0:
                return None

        # Do the same for std
        stds = np.std(sm_dict["Validation_Raw_Snps"], 1, dtype='float32')
        monomorphic_filter = stds > 0
        self._error_dict["Monomorphic"] = len(monomorphic_filter) - np.sum(monomorphic_filter)

        if self._filter_array(sm_dict, monomorphic_filter) == 0:
            return None
        else:
            self._error_dict["Accepted_Snps"] = len(sm_dict[self.sm_variants])

    def _seek_to_start(self, chromosome, file):
        """
        Seek to the start position of the summary file for this chromosome if chromosomes where in the summary file
        otherwise it will just skip the header
        """
        if chromosome == 1:
            start_line = file.readline()
            self._summary_last_position = len(start_line)
        else:
            file.seek(self._summary_last_position)

    def _sum_stats_frequencies(self, line):
        """
        This will check to see if the user is using Psychiatric Genomics Consortium Summary Stats and if so use the
        case/control frequency and N to construct the frequency. Otherwise, it will attempt to use the MAF column and
        if that is also not set it will return -1
        """

        # If we have Psychiatric Genomics Consortium Summary stats use these
        if (self.sm_case_freq is not None) and (self.sm_control_freq is not None):

            # If the frequency's are indexed but are NA then return -1
            if line[self.sm_case_freq] in (".", "NA") or line[self.sm_control_freq] in (".", "NA"):
                return -1

            # If the n for the case and control are known, then calculate a more accurate frequency
            elif (self.sm_case_n is not None) and (self.sm_control_freq is not None):
                if line[self.sm_case_n] in (".", "NA") or line[self.sm_control_n] in (",", "NA"):
                    return -1
                else:
                    case_n = float(line[self.case_n])
                    control_n = float(line[self.control_n])
                    tot_n = case_n + control_n
                    a_scalar = case_n / float(tot_n)
                    u_scalar = control_n / float(tot_n)
                    return (float(line[self.sm_case_freq]) * a_scalar) + (float(line[self.sm_control_freq]) * u_scalar)

            # If n is not know we just divide by 2
            else:
                return (float(line[self.sm_case_freq]) + float(line[self.sm_control_freq])) / 2.0

        # If we have MAF values in the summary stats then we can use those as the frequency
        elif self.sm_minor_allele_freq is not None:
            if line[self.sm_minor_allele_freq] not in (",", "NA"):
                return float(line[self.sm_minor_allele_freq])
            else:
                return -1

        # If we have no frequency information just return -1
        else:
            return -1

    def _error_dict_to_terminal(self, chromosome):
        """Print the error dict for this chromosome then reset the initialised to default 0"""
        print(f"\nCleaned summary statistics for chromosome: {chromosome}")
        for index, (k, v) in enumerate(zip(self._error_dict.keys(), self._error_dict.values())):
            if index == 0:
                print(Fore.LIGHTCYAN_EX + "{:<25} {}".format(k, v))
                print(Fore.LIGHTCYAN_EX + "----------------------------------------")
            else:
                print("{:<25} {}".format(k, v))

        # Reset values to 0
        for k, v in zip(self._error_dict.keys(), self._error_dict.values()):
            if isinstance(v, int):
                self._error_dict[k] = 0
