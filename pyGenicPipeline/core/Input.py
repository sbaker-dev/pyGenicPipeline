from pyGenicPipeline.utils import errors as ec
from pyGenicPipeline.utils import misc as mc
from .Loaders import *

from miscSupports import load_yaml, directory_iterator, flatten
from pysnptools.distreader import Bgen
from pysnptools.snpreader import Bed
from csvObject import CsvObject
from collections import Counter
from pyGenicParser import *
from pathlib import Path
import numpy as np
import pickle
import gzip
import re


class Input(CommonGenetic, SummaryLoader, ArgsParser):
    def __init__(self, args):
        super().__init__(args)

        # General operational parameters
        self.operation = self._set_current_job(self.args["Operation"])

        # Todo this is bad, needs to be removed at some point
        self._config = self.config["Other"]

        # Set filter information
        self.make_sub_directory("PGS", "Filter_Cleaned")
        self.filter_directory = Path(self.working_dir, "PGS", "Filter_Cleaned")
        self.maf_min = self._config["Min_Maf"]
        self.freq_discrepancy = self._config["Max_Freq_Discrepancy"]
        self.filter_index = self.args["Filter_Index"]
        self._filter_iter_size = self.args["Filter_Range"]
        self.clean_headers, self._clean_dict = self._set_cleaned_headers()

        # Gibbs information
        self.make_sub_directory("PGS", "Weights")
        self.weights_directory = Path(self.working_dir, "PGS", "Weights")
        self.gm = self._set_genome()
        self.ld_radius = self.args["LD_Radius"]
        self.herit_calculated = self.args["Heritability_Calculated"]
        self.gibbs_causal_fractions = self._set_causal_fractions()
        # todo set these up externally
        self.gibbs_run = False
        self.gibbs_iter = 100
        self.gibbs_burn_in = 10
        self.gibbs_shrink = 1
        self.gibbs_zero_jump = 0.01
        self.gibbs_random_seed = 42
        self.gibbs_tight = None
        self.gibbs_headers, self._gibbs_header_dict = self._set_gibbs_headers()
        self.gibbs_breaker = True

        # Score information
        self.make_sub_directory("PGS", "Chromosome_Scores")
        self.scores_directory = Path(self.working_dir, "PGS", "Chromosome_Scores")
        self.phenotype_file = mc.validate_path(self.args["Phenotype"])
        self.covariates_file = mc.validate_path(self.args["Covariates"])

        if (self.sm_case_freq is not None) or (self.sm_control_n is not None):
            raise NotImplementedError("Psychiatric Genomics Consortium Summary stats are untested and unfinished!")

    @staticmethod
    def _set_current_job(operation_dict):
        """
        Set the current job to be processed

        :param operation_dict:
            May be None if the user is using the main method as an object
            May be a string if the user is using a job submission form
            May be a dict if the user is using the method for development or natively in python
        :type operation_dict: None | str | dict

        :return:
            If None, Returns None
            If str, Returns operation_dict
            If Dict, Asserts that only 1 job is selected and then returns the job name as a string

        :raises TypeError, AssertionError:
            TypeError if job is not a None, str, or dict.
            AssertionError if the job dict contains more than a single job
        """
        if not operation_dict:
            return None
        elif isinstance(operation_dict, str):
            return operation_dict
        elif isinstance(operation_dict, dict):
            job = [job_name for job_name, run in zip(operation_dict.keys(), operation_dict.values()) if run]
            assert len(job) == 1, ec.job_violation(job)
            return job[0]
        else:
            raise TypeError(ec.job_type(type(operation_dict)))

    def _load_local_data(self, access_key):
        """
        This will set a dataset path that has been embedded into the package as a non yaml file sourced from LDPred if
        the current arg is set to be true

        :param access_key: The key to access the data file
        :type access_key: str

        :return: Path to the relevant file if request is not equal to None, else None
        :rtype: Path | None
        """

        if self.args[access_key]:
            package_root = Path(__file__).parent.parent

            if access_key == "HapMap3":
                access_path = Path(package_root, "Data", "hm3_sids.txt.gz")
                assert access_path, ec.path_invalid(access_path, "_load_local_data")
                return access_path
            elif access_key == "Filter_Long_Range_LD":
                access_path = Path(package_root, "Data", "long-range-ld-price-2008hg38.txt")
                assert access_path, ec.path_invalid(access_path, "_load_local_data")
                return access_path
            else:
                raise Exception(f"Unknown Key provided to _load_local_data: {access_key}")

        else:
            return None

    def validation_chromosomes(self):
        """
        This will create a dataset of all the chromosomes that we have to work with

        :return: A list of valid chromosomes
        :rtype: list
        """

        valid_chromosomes = []
        for file in directory_iterator(self.gen_directory):
            if Path(self.gen_directory, file).suffix == self.gen_type:
                valid_chromosomes.append(int(re.sub(r'[\D]', "", Path(self.gen_directory, file).stem)))
        valid_chromosomes.sort()
        return valid_chromosomes

    def validate_chromosomes(self, chromosome_set):
        """
        It is possible for non valid chromosomes, this will validate for numeric or known maps from str chromosomes to
        numeric, for example X: 23,  and flag and error if it fails to find a map.
        """
        ok_chromosomes = []
        for chromosome in chromosome_set:
            try:
                ok_chromosomes.append(int(chromosome))
            except ValueError:
                try:
                    ok_chromosomes.append(self._chromosome_map[chromosome])
                except KeyError:
                    raise Exception(f"Found chromosome {chromosome} which could not be mapped!")

        return np.unique(ok_chromosomes)

    def select_file_on_chromosome(self):
        """
        For a given chromosome, get the respective file from the genetic directory

        :return: Path to the current file as a Path from pathlib
        """
        for file in directory_iterator(self.gen_directory):
            if Path(self.gen_directory, file).suffix == self.gen_type:
                try:
                    if int(re.sub(r'[\D]', "", Path(self.gen_directory, file).stem)) == self.target_chromosome:
                        return str(Path(self.gen_directory, file).absolute())
                except (ValueError, TypeError):
                    continue

        raise Exception(f"Failed to find any relevant file for {self.target_chromosome} in {self.gen_directory}")

    @staticmethod
    def _set_validation_size(validation_size):
        """
        If Validation_Size is set, validate it is between 0 and 1 and then return it, otherwise None.

        :param validation_size: The size of the validation group
        :type validation_size: None | float

        :return: A None or a float
        :rtype: None | float
        """

        if not validation_size:
            return None
        else:
            assert 0 <= float(validation_size) <= 1, ec.validation_size_invalid(validation_size)
            return float(validation_size)

    def gen_reference(self, load_path):
        """Get the pysnptools reference via the load type"""
        if self.gen_type == ".bed":
            return Bed(load_path, count_A1=True)
        elif self.gen_type == ".bgen":
            if self._snp_tools:
                return Bgen(load_path)
            else:
                return BgenObject(load_path)
        else:
            raise Exception("Unknown load type set")

    def _set_validation_sample_size(self, full_sample_size):
        """
        This will return the value of the validation in terms of individuals rather than a percentage that the user
        specified

        :param full_sample_size: An integer of the number of samples in the full sample
        :type full_sample_size: int

        :return: Integer of the number of samples required in the validation sample
        :rtype: int
        """
        return int((full_sample_size * self.population_percent) * self.validation_size)

    def construct_validation(self, load_path):
        """
        We need to construct a validation sample from the percentage the user provided and the iid_count, this then
        returns this slice of the sample from the start up to this percentage (Uses int so may be slightly above or
        below the percentage provided based on rounding / floating point errors), and then the rest of the sample of the
        core set

        :param load_path: Path to the relevant load file
        :return: The validation and ref sample class holders
        """

        # todo Before splitting in to validation and core, allow a sample size modifier to remove people (ie for ukb)
        # Set validation and core sets of sids based on the load type
        validation_size = self._set_validation_sample_size(self.gen_reference(load_path).iid_count)
        validation = self.gen_reference(load_path)[:validation_size, :]
        ref = self.gen_reference(load_path)[validation_size:, :]
        return validation, ref

    def load_variants(self, load_path, validation, ref):
        """
        Load variants, for .bgen or plink files, as a set of snps that exist within the current chromosome. Uses the
        validation percentage to construct a validation group, and returns the set of snps for each group. If hap_map_3
        is enabled, it will strip out snps not in hap_map_3.

        We will also need a way to index out the variant information, so we set the indexer according to the load type

        :return: Set of the validation and core set of snps, as well as an indexer to extract information from them
        """

        #  Set validation and core sets of sids based on the load type
        if self.gen_type == ".bed":
            validation = validation.sid
            ref = ref.sid
            indexer = PlinkObject(load_path, True)

        elif self.gen_type == ".bgen":
            indexer = BgenObject(load_path)
            if self._snp_tools:
                print("Loading bgen with PySnpTools\n")
                # Bgen files store [variant id, rs_id], we just want the rs_id hence the [1]; see https://bit.ly/2J0C1kC
                validation = [snp.split(",")[1] for snp in validation.sid]
                ref = [snp.split(",")[1] for snp in ref.sid]

            else:
                print("Loading bgen with custom pybgen via pyGenicParser\n")
                validation = validation.sid_array()
                ref = ref.sid_array()

        else:
            raise Exception("Unknown load type set")

        # If we only want the hap_map_3 snps then check each snp against the set of hap_map_3
        if self.hap_map_3:
            hap_map_3_snps = self.load_hap_map_3()
            validation = [snp for snp in validation if snp in hap_map_3_snps]
            ref = [snp for snp in ref if snp in hap_map_3_snps]

        # Check for duplicates that may be loaded later in the pipeline if we don't filter them out and will otherwise
        # not be detected due to returning a set
        v_count = len(validation)
        validation = [snp for (snp, count) in Counter(validation).items() if count == 1]

        r_count = len(ref)
        ref = [snp for (snp, count) in Counter(ref).items() if count == 1]

        # Count total duplicates
        duplicates = np.sum([(v_count - len(validation)) + (r_count - len(ref))])
        return set(flatten([validation, ref])), indexer, duplicates

    def variant_names(self, sm_dict):
        """Variant names differ in pysnptools bgen, so account for this and just return rs_id's"""
        return mc.variant_array(self.snp_id.lower(), sm_dict[self.sm_variants])

    def chunked_snp_names(self, sm_dict, chunk_return=False):
        """
        Even a couple of 10's of thousands of snps will lead to memory issues especially if there are large numbers of
        individuals in the data set. This will load the variant names that have been cleaned, then chunk them by a
        dimension calculated from the number of variant names by the filter_iter_size set by the user.

        :param sm_dict: Dict of cleaned summary statistics
        :type sm_dict: dict

        :param chunk_return: If True returns the chunk size as well as the variants, otherwise just the variants
        :type chunk_return: bool

        :return:  np.array_split numpy arrays of snp names, or this combined with the chunk size
        """

        # Extract variant names
        variant_names = self.variant_names(sm_dict)

        # Calculate the number of chunks required, then return the variant names split on chunk size
        chunks = int(np.ceil(len(variant_names) / self._filter_iter_size))
        if chunk_return:
            return np.array_split(variant_names, chunks), chunks
        else:
            return np.array_split(variant_names, chunks)

    def isolate_raw_snps(self, gen_file, variant_names):
        """
        This will isolate the raw snps for a given bed or bgen file

        :param gen_file: Genetic file you wish to load from
        :param variant_names: The snp names to isolate the dosage for
        :return: raw snps
        """

        # bed returns 2, 1, 0 rather than 0, 1, 2 although it says its 0, 1, 2; so this inverts it
        if self.gen_type == ".bed":
            ordered_common = gen_file[:, gen_file.sid_to_index(variant_names)].read(dtype=np.int8).val
            raw_snps = np.array([abs(snp - 2) for snp in ordered_common.T], dtype=np.int8)

        # We have a [1, 0, 0], [0, 1, 0], [0, 0, 1] array return for 0, 1, 2 respectively. So if we multiple the arrays
        # by their index position and then sum them we get [0, 1, 2]
        elif self.gen_type == ".bgen":
            if self._snp_tools:
                # Re-construct the variant_id-rs_id
                v_dict = {snp[1]: snp[0] for snp in [snp.split(",") for snp in gen_file.sid]}
                variant_names = [f"{v_dict[rs_id]},{rs_id}" for rs_id in variant_names]

                ordered_common = gen_file[:, gen_file.sid_to_index(variant_names)].read(dtype=np.int8).val
                raw_snps = sum(np.array([snp * i for i, snp in enumerate(ordered_common.T)], dtype=np.int8))
            else:
                raw_snps = gen_file.dosage_from_sid(variant_names)

        else:
            raise Exception(f"Critical Error: Unknown load type {self.gen_type} found in _isolate_dosage")

        assert len(raw_snps) == len(variant_names), "Failed to filter out duplicates"
        return raw_snps

    def normalise_snps(self, gen_file, variant_names, std_return=False):
        """For gibbs we use normalised snps, this process will use the information we have filtered to construct it"""

        raw_snps = self.isolate_raw_snps(gen_file, variant_names)

        # Get the number of snps and individuals in the filtered dict
        n_snps, n_individuals = raw_snps.shape

        # Need to reformat the shape to construct the normalised snps
        raw_means = np.mean(raw_snps, 1, dtype='float32')
        raw_means.shape = (n_snps, 1)
        raw_stds = np.std(raw_snps, 1, dtype='float32')
        raw_stds.shape = (n_snps, 1)

        # Use this information to construct a normalised snps
        normalised_snps = np.array((raw_snps - raw_means) / raw_stds, dtype="float32")
        assert normalised_snps.shape == raw_snps.shape

        if std_return:
            return normalised_snps, raw_stds
        else:
            return normalised_snps, None

    def genetic_phenotypes(self, gen_file, load_path):
        """
        Load the full genetic data for this chromosome and isolate any information that can be isolated from it. In this
        case, .bed load types can access more than bgen due to the ability to extract sex from the .fam file.
        """

        ph_dict = {}
        # For plink files, load the fam file then extract the fid, iid and sex information
        if self.gen_type == ".bed":
            ph_dict[self.fam] = np.array(PlinkObject(load_path).get_family_identifiers())
            ph_dict[self.fid] = mc.variant_array(self.fid.lower(), ph_dict[self.fam])
            ph_dict[self.iid] = mc.variant_array(self.iid.lower(), ph_dict[self.fam])
            ph_dict.pop(self.fam, None)

        # Bgen doesn't have a fam equivalent, so just load the fid and iid
        elif self.gen_type == ".bgen":
            # todo update to allow for sex and missing if we have loaded .sample
            if self._snp_tools:
                ids = gen_file.iid
            else:
                # todo - This won't work. Ids only returns iid not iid + fid
                ids = gen_file.iid_array()

            ph_dict[self.fid] = np.array([fid for fid, iid in ids])
            ph_dict[self.iid] = np.array([iid for fid, iid in ids])

        else:
            raise Exception("Unknown load type set")

        return ph_dict

    def load_hap_map_3(self):
        """
        Users may wish to limit valid snps to those found within HapMap3. If they do, we access them via the local file
        and return them as a set
        """

        # If the HapMap3 file exists, then extract the snp ids as a set and return them
        f = gzip.open(self.hap_map_3, 'r')
        hm3_sids = pickle.load(f)
        f.close()
        return set(hm3_sids)

    def load_lr_ld_dict(self):
        """
        This will read in the long rang ld dict from Price et al. AJHG 2008 long range LD tables taken from LDPred and
        then filter out the keys relevant to this chromosome.
        """
        long_dict = {chromosome_key: {} for chromosome_key in range(1, 24)}
        with open(str(self.lr_ld_path.absolute()), 'r') as f:
            for line in f:
                chromosome_line, start_pos, end_pos, hild = line.split()
                try:
                    long_dict[int(chromosome_line)][hild] = {'start_pos': int(start_pos), 'end_pos': int(end_pos)}
                except ValueError:
                    continue
        return long_dict





    def _set_causal_fractions(self):
        """
        If the user has provided a set of causal fractions of variants to use for the gibbs sampler then use those, else
        use the default that LDPred used.
        """
        if self.args["Causal_Fractions"]:
            return self.args["Causal_Fractions"]
        else:
            return [1, 0.3, 0.1, 0.03, 0.01, 0.003, 0.001]

    def sm_dict_from_csv(self, load_path, different_types=None):
        """
        Load a saved cleaned file from summary statistics for a given chromosome found at the load_path provided, and
        use this to construct the sm_dict we pass between methods
        """
        header_length = len(CsvObject(load_path).headers)
        if header_length == len(self.cleaned_types):
            load_file = CsvObject(load_path, self.cleaned_types, set_columns=True)
        elif header_length == len(self.cleaned_types) - 1:
            load_file = CsvObject(load_path, self.cleaned_types[:-1], set_columns=True)
        elif different_types:
            print("No sm dict construction - returning headers and column values")
            load_file = CsvObject(load_path, different_types, set_columns=True)
            return {header: np.array(columns) for header, columns in zip(load_file.headers, load_file.column_data)}
        else:
            print("No Know specification - returning untyped headers and column values")
            load_file = CsvObject(load_path, set_columns=True)
            return {header: columns for header, columns in zip(load_file.headers, load_file.column_data)}

        chromo = load_file.column_data[self.c_chromosome]
        bp_pos = load_file.column_data[self._clean_dict[self.bp_position]]
        snp_id = load_file.column_data[self._clean_dict[self.snp_id]]
        effect = load_file.column_data[self._clean_dict[self.effect_allele]]
        alt = load_file.column_data[self._clean_dict[self.alt_allele]]
        log = load_file.column_data[self._clean_dict[self.log_odds]]
        beta = load_file.column_data[self.c_beta]
        freq = load_file.column_data[self._clean_dict[self.freq]]

        sm_variants = [Variant(ch, bp, sn, ef, al) for ch, bp, sn, ef, al in zip(chromo, bp_pos, snp_id, effect, alt)]
        return {self.sm_variants: np.array(sm_variants), self.log_odds: np.array(log), self.beta: np.array(beta),
                self.freq: np.array(freq)}

    @property
    def _chromosome_map(self):
        """
        By default we assume chromosome are numeric, but can catch X and turn it into 23. If someone wants to specific
        a new map, for example if they have a y chromosome and want to map it to a number, they can do so with a custom
        dict. Otherwise we just return the dict that was embedded in ldpred
        """
        if self.args["Chromosome_Map"]:
            return self.args["Chromosome_Map"]
        else:
            return {"X": 23, "chr_x": 23, "chrom_x": 23}

    def _set_cleaned_headers(self):
        """Construct headers to be used for writing and reading cleaned files"""
        cleaned_headers = [self.chromosome, self.bp_position, self.snp_id, self.effect_allele, self.alt_allele,
                           self.log_odds, self.beta, self.freq, self.ld_scores]

        cleaned_dict = {header: i for i, header in enumerate(cleaned_headers)}
        return cleaned_headers, cleaned_dict

    def _set_gibbs_headers(self):
        """Construct the headers that will be used in the writing of weights"""

        gibbs_headers = [self.chromosome, self.bp_position, self.snp_id, self.effect_allele, self.alt_allele,
                         self.beta, self.log_odds, self.ld_scores, self.gibbs_beta, self.effect_size]

        gibbs_dict = {header: i for i, header in enumerate(gibbs_headers)}

        return gibbs_headers, gibbs_dict

    def local_values(self, values, snp_index, number_of_snps):
        """
        We want to construct a window of -r + r around each a given list of values where r is the radius. However, the
        first r and last N-r of the snps will not have r number of snps before or after them so we need to account for
        this by:

        Taking the maximum of (0, i-r) so that we never get a negative index
        Taking the minimum of (n_snps, (i + radius + 1)) to ensure we never get an index out of range

        :param values: A set of values to extract a local off
        :param snp_index: Index
        :param number_of_snps: total number of snps

        :return: An array of shape snps of a maximum of 'radius' number of snps surrounding the current snp accessed via
            index.
        """
        return values[max(0, snp_index - self.ld_radius): min(number_of_snps, (snp_index + self.ld_radius + 1))]

    def _set_genome(self):
        """Load the genome file if it has been produced, else return None"""
        genome_path = Path(self.working_dir, "genome_wide_config.yaml")
        if genome_path.exists():
            return load_yaml(genome_path)
        else:
            return None

    @staticmethod
    def check_sm_dict(sm_dict):
        """Validate that sm dict still exists"""
        if not sm_dict:
            raise Exception("Failed to find any snps!\n\n")

    # todo: Often, this will not be 9 but 8 long. We need to make it 8, then correct when its 9
    @property
    def cleaned_types(self):
        """The types of each column in the cleaned data"""
        return [int, int, str, str, str, float, float, float, float]


    @property
    def f_std(self):
        return f"{self.filter_key}_{self.stds}"

    @property
    def f_freq(self):
        return f"{self.filter_key}_{self.freq}"


    @property
    def c_chromosome(self):
        """Chromosome header index in Cleaned Data file"""
        return self._clean_dict[self.chromosome]

    @property
    def c_ld_score(self):
        """LD_Score header index in Cleaned Data file"""
        return self._clean_dict[self.ld_scores]

    @property
    def c_beta(self):
        """Beta header index in Cleaned Data file"""
        return self._clean_dict[self.beta]
