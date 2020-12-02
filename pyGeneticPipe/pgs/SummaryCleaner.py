from pyGeneticPipe.geneticParsers.variantObjects import Variant, Nucleotide
from pyGeneticPipe.geneticParsers.bgenObject import BgenObject
from pyGeneticPipe.geneticParsers.bimObject import BimObject
from pyGeneticPipe.utils import error_codes as ec
from pyGeneticPipe.utils import misc as mc
from pyGeneticPipe.core.Input import Input
from scipy import stats
import numpy as np
import time


class SummaryCleaner(Input):
    def __init__(self, args):
        super().__init__(args)

        # Depending on how this ends up being constructed we might want to log this out as a .txt rather than just
        # print this to a terminal.
        self._sum_error_dict = {"Removal case": "Count", "Invalid_Snps": 0, f"Miss Matching {self.chromosome}": 0,
                                f"Miss Matching {self.bp_position}": 0, f"Non Finite {self.effect_size}": 0,
                                f"Non Finite {self.p_value}": 0, f"Non Finite {self.standard_errors}": 0,
                                "Ambiguous_SNP": 0, "Non_Allowed_Allele": 0, "Non_Matching": 0}
        self._summary_last_position = 0

    def clean_summary_statistics(self, chromosome, load_path, validation, core):
        """
        This will take the summary statistics and access the validatable snps, found by cross referencing the genetic
        validation and core samples, and clean them of possible errors. It then returns a ordered on base pair position
        dictionary of information required for constructing poly-genetic scores
        """
        # Check for input arguments
        t0 = self._assert_clean_summary_statistics()
        print(f"Starting Chromosome: {chromosome}")

        # Clean the summary lines to only include validatable snps from our genetic samples that exit in this chromosome
        core_snps, sm_line, sm_variants, validation_snps = self._valid_snps_lines_and_variants(
            chromosome, core, load_path, validation)

        # Construct the summary dict with our summary lines and Variants objects of our valid snps
        sm_dict = {self.sm_lines: np.array(sm_line), self.sm_variants: np.array(sm_variants)}

        # Clean the summary lines of valid snps for potential errors, if we ever wipe all our samples return None
        sm_dict = self._validate_summary_lines(sm_dict)
        if not sm_dict:
            return None

        # Construct the order from the base pair position
        order = np.argsort(mc.variant_array(self.bp_position.lower(), sm_dict[self.sm_variants]))

        # Given we have only accepted snps that are within the validation / core, we should never have more snps in
        # summary than within the validation. If we do, something has gone critically wrong.
        assert len(order) <= len(validation_snps), ec.snp_overflow(len(order), len(validation_snps))
        assert len(order) <= len(core_snps), ec.snp_overflow(len(order), len(core_snps))

        # In this case we can order the array using filter array as well, and we return this ordered dict
        mc.filter_array(sm_dict, order)

        # Log to terminal what has been filtered / removed
        t1 = mc.error_dict_to_terminal(self._sum_error_dict)
        print(f"Cleaned summary stats for Chromosome {chromosome} in {round(t1 - t0, 2)} Seconds\n")
        return sm_dict

    def _valid_snps_lines_and_variants(self, chromosome, core, load_path, validation):
        """
        We will load our variants from our validation and core samples and use those to check if the snp found in the
        summary line is within our validation and core sample sets of snps. If this is the case, then we will add the
        line to sm_line as well as a Variant object of the current snp valid snp to sm_variants
        """
        validation_snps, core_snps, indexer = self._load_variants(load_path, validation, core)

        sm_variants = []
        sm_line = []
        with mc.open_setter(self.summary_file)(self.summary_file) as file:
            self._seek_to_start(chromosome, file)

            # For each line in the GWAS Summary file
            for index, line_byte in enumerate(file):

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
                        self._sum_error_dict["Invalid_Snps"] += 1

        print("Extracted snps from summary file\n")
        return core_snps, sm_line, sm_variants, validation_snps

    def _validate_summary_lines(self, sm_dict):
        """This will load in each possible header, and clean our dict of values by filtering"""

        # If we have chromosomes in our summary statistics check the chromosome of the snps against the validation
        if self.sm_chromosome is not None:
            if not self._validation_equality(self.sm_chromosome, self.chromosome, sm_dict):
                return None

        # If we have base pair position in our summary then validate the base pair
        if self.bp_position is not None:
            if not self._validation_equality(self.sm_bp_position, self.bp_position, sm_dict, int):
                return None

        # Clean the summary stats effect sizes for calculation of beta later
        if not self._validation_finite(sm_dict, self.sm_effect_size, self.effect_size):
            return None

        # Clean the P values
        if not self._validation_finite(sm_dict, self.sm_p_value, self.p_value):
            return None

        # If we are using z scores we need to load and clean the standard errors column
        if self.z_scores:
            if not self._validation_finite(sm_dict, self.sm_standard_errors, self.standard_errors):
                return None

        # Use the raw beta, standard errors, and p value if required to construct beta and beta_odds
        self._validation_betas(sm_dict)

        # Check that the nucleotides are sane and flip them if required
        sm_dict = self._validate_nucleotides(sm_dict)
        if not sm_dict:
            return None
        else:
            # Calculate the frequencies and set info if it exists
            sm_dict[self.freq] = np.array([self._sum_stats_frequencies(line) for line in sm_dict[self.sm_lines]])
            sm_dict[self.info] = self._validate_info(sm_dict[self.sm_lines])
            return sm_dict

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
        summary_array = mc.line_array(line_index, summary_dict[self.sm_lines], line_type)
        variant_array = mc.variant_array(variant_key.lower(), summary_dict[self.sm_variants])

        # Filter of True if the variant and summary match, else False which will remove this snp
        obj_filter = summary_array == variant_array

        self._sum_error_dict[f"Miss Matching {variant_key}"] = len(obj_filter) - np.sum(obj_filter)
        return mc.filter_array(summary_dict, obj_filter)

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
        summary_dict[summary_key] = mc.line_array(line_index, summary_dict[self.sm_lines], float)

        # Filter out anything that is not finite or is equal to zero
        obj_filter = np.array([True if np.isfinite(obj) and obj != 0 else False for obj in summary_dict[summary_key]])
        self._sum_error_dict[f"Non Finite {summary_key}"] = len(obj_filter) - np.sum(obj_filter)
        return mc.filter_array(summary_dict, obj_filter)

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
        effected_allele = mc.line_array(self.sm_effect_allele, sm_dict[self.sm_lines])
        alt_allele = mc.line_array(self.sm_alt_allele, sm_dict[self.sm_lines])
        sm_dict[self.nucleotide] = np.array([Nucleotide(e, a) for e, a in zip(effected_allele, alt_allele)])

        # Filter out any snps where the summery or variant Nucleotide is ambiguous
        filter_ambiguous = [False if (sm_nuc.to_tuple() in self.ambiguous_snps) or
                                     ((var_nuc.a1, var_nuc.a2) in self.ambiguous_snps)
                            else True
                            for sm_nuc, var_nuc in zip(sm_dict[self.nucleotide], sm_dict[self.sm_variants])]
        self._sum_error_dict["Ambiguous_SNP"] = len(filter_ambiguous) - np.sum(filter_ambiguous)
        if not mc.filter_array(sm_dict, filter_ambiguous):
            return None

        # Sanity Check
        # Filter out any snps that do not pass a sanity check (Only a t c and g)
        allowed_filter = [False if (sm_nuc.a1 not in self.allowed_alleles) or
                                   (sm_nuc.a2 not in self.allowed_alleles) or
                                   (var_nuc.a1 not in self.allowed_alleles) or
                                   (var_nuc.a2 not in self.allowed_alleles)
                          else True
                          for sm_nuc, var_nuc in zip(sm_dict[self.nucleotide], sm_dict[self.sm_variants])]
        self._sum_error_dict["Non_Allowed_Allele"] = len(allowed_filter) - np.sum(allowed_filter)
        if not mc.filter_array(sm_dict, allowed_filter):
            return None

        # Determine Flipping
        # Construct a flip status of 1, 0, -1 for No flipping, failed flipping, and flipped successfully which we can
        # multiple our betas by
        sm_dict["Flip"] = np.array([self._flip_nucleotide(var_nuc, sm_nuc) for var_nuc, sm_nuc in
                                    zip(sm_dict[self.sm_variants], sm_dict[self.nucleotide])])
        filter_flipped = np.array([False if flipped == 0 else True for flipped in sm_dict["Flip"]])
        self._sum_error_dict["Non_Matching"] = len(filter_flipped) - np.sum(filter_flipped)
        if not mc.filter_array(sm_dict, filter_flipped):
            return None

        # Now we have filtered away any errors, multiple the dicts beta and log_odds elements by 1 or -1 based on no
        # flipping or requiring flipping
        sm_dict[self.beta] = sm_dict[self.beta] * sm_dict["Flip"]
        sm_dict[self.log_odds] = sm_dict[self.log_odds] * sm_dict["Flip"]
        return sm_dict

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

    def _validate_info(self, sm_line):
        """Construct infos if they exist in the summary stats else return an array of length of summary dict"""
        if self.sm_info is not None:
            infos = mc.line_array(self.sm_info, sm_line, float)
        else:
            infos = np.empty(len(sm_line))
            infos.fill(-1)
        return infos

    def _load_variants(self, load_path, validation, core):
        """
        Load variants, for .bgen or plink files, as a set of snps that exist within the current chromosome. Uses the
        validation percentage to construct a validation group, and returns the set of snps for each group. If hap_map_3
        is enabled, it will strip out snps not in hap_map_3.

        We will also need a way to index out the variant information, so we set the indexer according to the load type

        :param load_path: Current Chromosome file
        :return: Set of the validation and core set
        """

        #  Set validation and core sets of sids based on the load type
        if self.load_type == ".bed":
            validation = validation.sid
            core = core.sid
            indexer = [BimObject(load_path).construct_index(), BimObject(load_path)]

        elif self.load_type == ".bgen":
            # Bgen files store [variant id, rsid], we just want the rsid hence the [1]; see https://bit.ly/2J0C1kC
            validation = [snp.split(",")[1] for snp in validation.sid]
            core = [snp.split(",")[1] for snp in core.sid]
            indexer = [BgenObject(load_path).index_of_snps(), BgenObject(load_path)]

        else:
            raise Exception("Unknown load type set")

        # If we only want the hap_map_3 snps then check each snp against the set of hap_map_3
        if self.hap_map_3:
            hap_map_3_snps = self.load_hap_map_3()
            validation = [snp for snp in validation if snp in hap_map_3_snps]
            core = [snp for snp in core if snp in hap_map_3_snps]

        return set(validation), set(core), indexer

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
        assert self.summary_file, ec.missing_arg(self.operation, "Summary_Path")
        assert self.load_type, ec.missing_arg(self.operation, "Load_Type")
        assert self.load_directory, ec.missing_arg(self.operation, "Load_Directory")
        assert self.validation_size, ec.missing_arg(self.operation, "Validation_Size")

        return time.time()

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

    def _seek_to_start(self, chromosome, file):
        """
        If we are running without multi-core positioning then we can log the differing seek values to jump to the next
        position when calling the next chromosome. Otherwise it will just skip the header
        """
        if self.multi_core_splitter:
            start_line = file.readline()
            self._summary_last_position = len(start_line)
        else:
            if chromosome == 1:
                start_line = file.readline()
                self._summary_last_position = len(start_line)
            else:
                file.seek(self._summary_last_position)
