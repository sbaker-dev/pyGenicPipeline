from pyGeneticPipe.geneticParsers.supportObjects import Variant, Nucleotide, SMVariant
from pyGeneticPipe.geneticParsers.plink.plinkObject import PlinkObject
from pyGeneticPipe.geneticParsers.bgen.bgenObject import BgenObject
from pyGeneticPipe.utils import error_codes as ec
from pyGeneticPipe.utils import misc as mc
from pyGeneticPipe.core.Input import Input
from pysnptools.distreader import Bgen
from pysnptools.snpreader import Bed
from operator import itemgetter
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

        self._error_dict = {"Removal case": "Count", "Invalid_Snps": 0, "Chromosome": 0, "Position": 0,
                            "Effect_Size": 0, "P_Value": 0, "Standard_Errors": 0, "Duplicate_Position": 0,
                            "Ambiguous_SNP": 0, "Non_Matching": 0, "Non_Allowed_Allele": 0, "Filtered_Frequency": 0,
                            "Filtered_MAF": 0, "Monomorphic": 0, "Accepted_Snps": 0}
        self._summary_last_position = 0

    def clean_summary_statistics(self):

        # Note - this is basically becoming the -main- of prs, so will want to extract the chromosome bit so that it can
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
            sm_variants = self._clean_summary_stats(load_path, validation, core, chromosome)

            # Filter the summary stats
            self._filter_snps(load_path, sm_variants)
            return

    def _clean_summary_stats(self, load_path, validation, core, chromosome):
        """
        This will take the validation and core sample of snps, and check the snp against both sets. If the snp exists in
        the validation files, then it will go to cleaning the summary statistics for this chromosome line by line.
        """

        validation_snps, core_snps, indexer = self._load_variants(load_path, validation, core)

        sm_variants = []
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
                    sm_variant = self._validate_summary_line(line, self._set_variant(snp_id, indexer))
                    if sm_variant:
                        sm_variants.append(sm_variant)

                else:
                    # If the chromosomes exist in summary statistics we can terminate this for loop when we are no
                    # longer in the right zone and set tell to seek to this position for the next chromosome
                    if (self.sm_chromosome is not None) and (int(line[self.sm_chromosome]) > chromosome):
                        self._summary_last_position = file.tell() - len(line_byte)
                        file.close()
                        break
                    else:
                        self._error_dict["Invalid_Snps"] += 1

        # Given we have only accepted snps that are within the validation / core, we should never have more snps in
        # summary than within the validation. If we do, something has gone critically wrong.
        assert len(sm_variants) <= len(validation_snps), ec.snp_overflow(len(sm_variants), len(validation_snps))
        assert len(sm_variants) <= len(core_snps), ec.snp_overflow(len(sm_variants), len(core_snps))

        # We then need to order the snps on the base pair position
        variant_by_bp = [[variant, variant.bp_position] for variant in sm_variants]
        return np.array([variant for variant, bp in sorted(variant_by_bp, key=itemgetter(1))])

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
        if (self.sm_chromosome is not None) and line[self.sm_chromosome] != variant.chromosome:
            self._error_dict["Chromosome"] += 1
            return None

        # If we have base pair position in our summary then validate the base pair
        if (self.sm_bp_position is not None) and int(line[self.sm_bp_position]) != variant.bp_position:
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
        if self.sm_info is not None:
            info = float(line[self.sm_info])

        # If we have frequency, then extract it
        frequency = self._sum_stats_frequencies(line)

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

    def _isolate_raw_snps(self, gen_file, sm_variants):
        """
        This will isolate the raw snps for a given bed or bgen file

        :param gen_file: Genetic file you wish to load from
        :param sm_variants: summary statistics variants
        :return: raw snps
        """
        ordered_common = gen_file[:, gen_file.sid_to_index(self._extract_variant_name(sm_variants))].read().val

        # bed returns 2, 1, 0 rather than 0, 1, 2 although it says its 0, 1, 2; so this inverts it
        if self.load_type == ".bed":
            return np.array([abs(snp - 2) for snp in ordered_common.T])

        # We have a [1, 0, 0], [0, 1, 0], [0, 0, 1] array return for 0, 1, 2 respectively. So if we multiple the arrays
        # by their index position and then sum them we get [0, 1, 2]
        elif self.load_type == ".bgen":
            return sum(np.array([snp * i for i, snp in enumerate(ordered_common.T)]))

        else:
            raise Exception(f"Critical Error: Unknown load type {self.load_type} found in _isolate_dosage")

    def _extract_variant_name(self, sm_variants):
        """
        Different file types have different naming standards.

        .bed: ["rs123", "rs124", ... "rsN"]
        .bgen: ["rs123,rs123", "rs124,rs124", ... "rsN,rsN"]

        This will standardise the names to be a list of type equivalent to bed
        :param sm_variants: list of SMVariant
        :return: list of snp names
        """
        if self.load_type == ".bed":
            return [variant.variant_id for variant in sm_variants]
        elif self.load_type == ".bgen":
            print("Bgen load type, so need to restructure return type ... will take a bit longer longer!")
            return [variant.bgen_variant_id() for variant in sm_variants]
        else:
            raise Exception(f"Critical Error: Unknown load type {self.load_type} found in _isolate_dosage")

    def _filter_snps(self, load_path, sm_variants):
        # testing validation for comparision TEMP
        validation = Bed(load_path, count_A1=True)

        # for each core / validation we need to do _isoalte dosage
        validation_raw_snps = self._isolate_raw_snps(validation, sm_variants)
        # print(np.mean(validation_raw_snps, axis=1))

        freqs = np.sum(validation_raw_snps, 1, dtype='float32') / (2 * float(validation.iid_count))

        print("PRE FILTER")
        print(len(sm_variants), len(validation_raw_snps), len(freqs))

        # If the frequencies in the summary stats are not just a list of -1 errors then screen genetic snp frequencies
        summary_freqs = np.array([variant.frequency for variant in sm_variants])
        if (self.freq_discrepancy < 1) and (np.sum(summary_freqs == -1) != len(summary_freqs)):
            freq_filter = np.logical_or(np.absolute(summary_freqs - freqs) < self.freq_discrepancy,
                                        np.absolute(summary_freqs + (freqs - 1)) < self.freq_discrepancy)
            freq_filter = np.logical_or(freq_filter, summary_freqs <= 0)
            filtered_freqs = len(freq_filter) - np.sum(freq_filter)
            assert filtered_freqs >= 0
            self._error_dict["Filtered_Frequency"] += filtered_freqs

            sm_variants = sm_variants[freq_filter]
            validation_raw_snps = validation_raw_snps[freq_filter]
            freqs = freqs[freq_filter]

        if self.maf_min > 0:
            maf_filter = (freqs > self.maf_min) * (freqs < (1 - self.maf_min))
            maf_filtered_snps = len(maf_filter) - np.sum(maf_filter)
            assert maf_filtered_snps >= 0
            self._error_dict["Filtered_MAF"] += maf_filtered_snps

            # Filter out anything that fails maf
            sm_variants = sm_variants[maf_filter]
            validation_raw_snps = validation_raw_snps[maf_filter]
            freqs = freqs[maf_filter]

        # Do the same for std
        stds = np.std(validation_raw_snps, 1, dtype='float32')
        monomorphic_filter = stds > 0
        monomorphic_filtered_snps = len(monomorphic_filter) - np.sum(monomorphic_filter)
        assert monomorphic_filtered_snps >= 0
        self._error_dict["Monomorphic"] += monomorphic_filtered_snps
        sm_variants = sm_variants[monomorphic_filter]
        validation_raw_snps = validation_raw_snps[monomorphic_filter]
        freqs = freqs[monomorphic_filter]

        self._error_dict["Accepted_Snps"] = len(sm_variants)

        print("POST MONOMORPHIC FILTER")
        print(len(sm_variants), len(validation_raw_snps), len(freqs))

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
