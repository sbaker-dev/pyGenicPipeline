from pyGeneticPipe.utils import error_codes as ec
from pyGeneticPipe.utils import misc as mc
from pysnptools.distreader import Bgen
from pysnptools.snpreader import Bed
from pathlib import Path
import numpy as np
import pickle
import gzip
import h5py
import re


class Input:
    def __init__(self, args):
        # General operational parameters
        self.args = self._set_args(args)
        self._config = mc.load_yaml(Path(Path(__file__).parent, "Keys.yaml"))
        self.debug = self.args["Debug"]
        self.working_dir = self._validate_path(self.args["Working_Directory"], False)
        self.operation = self._set_current_job(self.args["Operation"])
        self.multi_core_splitter = self.args["Multi_Core_Splitter"]

        # The project file for this project
        self.project_name = self.args['Project_Name']
        self.project_file = self._create_project_file()
        self.summary_file = self._validate_path(self.args["Summary_Path"])
        self.load_directory = self._validate_path(self.args["Load_Directory"])
        self.hap_map_3 = self._load_local_data("HapMap3")
        self.lr_ld_path = self._load_local_data("Filter_Long_Range_LD")
        self.load_type = self.args["Load_Type"]
        self.validation_size = self._set_validation_size(self.args["Validation_Size"])

        # Set summary statistics information if required
        self.zipped, self.sample_size = self._set_summary_stats()
        self._summary_headers = self._set_summary_headers()
        self.effect_type = self._set_effect_type(self.args["Summary_Effect_Type"])
        self.z_scores = self._set_z_scores(self.args["Z_Scores"])
        self.ambiguous_snps, self.allowed_alleles, self.allele_flip = self._configure_alleles()
        self.maf_min = self._config["min_maf"]
        self.freq_discrepancy = self._config["max_freq_discrepancy"]

        if (self.sm_case_freq is not None) or (self.sm_control_n is not None):
            raise NotImplementedError("Psychiatric Genomics Consortium Summary stats are untested and unfinished!")

    @staticmethod
    def _set_args(args):
        """
        Args may be set as a dict, or as a yaml file with its path passed as the args. If the later then we need to load
        in the dict of values
        """
        if isinstance(args, dict):
            return args
        else:
            yaml_path = Path(args)
            assert (yaml_path.exists() and yaml_path.suffix == ".yaml"), ec.path_invalid(yaml_path, "_set_args")

            return mc.load_yaml(yaml_path)

    @staticmethod
    def _set_current_job(operation_dict):
        """
        Set the current job from a dict of possible jobs or a string of the job
        :param operation_dict: If A dict, each key is a method_call to be done via getattr and the value is a True or
            False bool. Only one job should be true for each process. If a string, then just the string of method call
        :return: the current job name to be processed via getattr
        """
        if not operation_dict:
            return None
        elif isinstance(operation_dict, str):
            return operation_dict
        else:
            job = [job_name for job_name, run in zip(operation_dict.keys(), operation_dict.values()) if run]
            assert len(job) == 1, ec.job_violation(job)
            return job[0]

    @staticmethod
    def _validate_path(path, allow_none=True):
        """
        We have multiple types of files and directories, some may be allow to be None as they will not be required
        whilst others like the working directory will always be required. This method is a generalisation of individual
        setters.

        :param path: Path to a directory or file
        :type path: str

        :param allow_none: Defaults to True, if true if a path is set to none it will just return None. If False, an
            assertion will be run to validate that it is not none. In both cases, should the file not be None, then the
            path is validated via Path.exists()
        :type allow_none: Bool

        :return: Path to the current file or directory if None return is not allowed, otherwise the Path return is
            optional and the return may be none.
        """
        if allow_none and not path:
            return None
        else:
            assert path and Path(path).exists(), ec.path_invalid(path, "_validate_path")
            return Path(path)

    def _load_local_data(self, access_key):
        """
        This will load a dataset that has been embedded into the package as a non yaml file sourced from LDPred
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

    def _configure_alleles(self):
        """
        Yaml storage of tuples/ sets didn't work so this configures a list of lists into a set of tuples for ambiguous,
        sets a set of allow alleles, and also returns the dict of allele_flip
        """
        ambiguous = set(tuple(ambiguous) for ambiguous in self._config["ambiguous_snps"])
        return ambiguous, set(self._config["allowed_alleles"]), self._config["allele_flip"]

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

    def _set_summary_stats(self):
        """
        If we are reading in the summary statistics file then validate its path, construct a valid set of snps in a set
        and map the chromosome and bp_position to each valid snp in a dict. Validate the type of zipped structure and
        the sample size of the study

        :return: summary_path, snp_pos_map, valid_snp, gz_status, sample_size
        """

        # Stop if not required
        if not self.summary_file:
            return None, None

        # Check the sample size from this study has been provided
        sample_size = self.args["Summary_Sample_Size"]
        assert (sample_size is not None) and (sample_size > 0), ec.sample_size()

        # Determine if the summary file is g-zipped
        gz_status = (self.summary_file.suffix == ".gz")
        return gz_status, sample_size

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

    def select_file_on_chromosome(self, chromosome):
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

    def _check_header(self, sum_header, headers, summary_headers):
        """
        We need to standardise our headers, and locate where the current header is in our summary file in terms of a
        base zero index.

        :param sum_header: Standardised header
        :param headers: summary statistics headers
        :return: None if not found else the index of the header in our file for this standardised header
        :rtype: None | int
        """
        header_indexes = [i for i, h in enumerate(headers) if h in summary_headers[sum_header]]

        assert len(header_indexes) < 2, ec.ambiguous_header(sum_header, headers, summary_headers[sum_header])
        if len(header_indexes) == 0:
            assert sum_header not in self._config["mandatory_headers"], ec.mandatory_header(
                sum_header, headers, summary_headers[sum_header])
            return None
        else:
            return header_indexes[0]

    def _set_summary_headers(self):
        """
        We may have users using custom headers, or they may be using a format we already have covered

        Note
        -----
        In GUI we basically want to call the header check to make sure we align columns correctly. If not they can set
        it themselves

        :return: The headers to validate
        """
        if not self.summary_file:
            return None

        custom_headers = self.args["Custom_Summary_Header"]
        if custom_headers:
            # Recast so that the values are in a list so they can be checked by the same method as defaults
            header_sets = {key: [v] for key, v in zip(custom_headers.keys(), custom_headers.values())}
        else:
            # Based on known summary statistics from LDPred sum_stats_parsers.py
            header_sets = self._config["header_keys"]

        with mc.open_setter(self.summary_file)(self.summary_file) as file:

            # Determine if we have custom headers or not via _loaded_sum_headers
            raw_headers = file.readline()
            headers = {header: self._check_header(header, mc.decode_line(raw_headers, self.zipped), header_sets)
                       for header in header_sets}

            file.close()
            return headers

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

    def construct_validation(self, load_path):
        """
        We need to construct a validation sample from the percentage the user provided and the iid_count, this then
        returns this slice of the sample from the start up to this percentage (Uses int so may be slightly above or
        below the percentage provided based on rounding / floating point errors), and then the rest of the sample of the
        core set

        :param load_path: Path to the relevant load file
        :return: The validation and core sample class holders
        """

        # todo Before splitting in to validation and core, allow a sample size modifier to remove people (ie for ukb)
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

    def _set_effect_type(self, effect_type):
        """
        Set the effect type of the betas for GWAS summary stats if set
        """
        if effect_type:
            assert effect_type in self._config["effect_types"], ec.invalid_effect_type(
                self._config["effect_types"], effect_type)
            return effect_type
        else:
            return None

    def _set_z_scores(self, set_z_scores):
        """
        If the user wants to compute z scores, then standard_errors most be set but otherwise it isn't a mandatory
        header.

        :param set_z_scores: A bool of if z scores should be calculated or not
        :type set_z_scores: bool

        :return: True if assertion of standard errors is also True, if set_z_scores is None then return None
        :rtype: None | bool
        """
        if set_z_scores:
            assert self.sm_standard_errors is not None, ec.z_scores_with_standard_errors
            return True
        else:
            return None

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

    def _create_project_file(self):
        """
        Setup the hdf5 file for this project

        :rtype: h5py.File
        """
        if not self.project_name:
            return None

        # Construct the path and validate it
        project_file = Path(self.working_dir, self.project_name)

        # Creates file if it doesn't exist, or overrides it, ie for fast debugging iterations, if set.
        if (project_file.exists() and self.args["Override"]) or not project_file.exists():
            return h5py.File(Path(self.working_dir, self.project_name), "w")

        # If we don't want to override it, we load it in appending mode
        elif project_file.exists() and not self.args["Override"]:
            return h5py.File(Path(self.working_dir, self.project_name), "a")

        # Should be impossible to reach here
        else:
            raise Exception(f"CRITICAL ERROR: Cannot load file - unknown reason")

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

    @property
    def chromosome(self):
        """Key used for accessing Chromosome headers, groups or other attributes"""
        return "Chromosome"

    @property
    def sm_chromosome(self):
        """Summary stat chromosome header index in GWAS summary file"""
        return self._summary_headers[self.chromosome]

    @property
    def snp_id(self):
        """Key used for accessing SNP_ID headers, groups or other attributes"""
        return "SNP_ID"

    @property
    def sm_snp_id(self):
        """Snp/variant id header index in GWAS summary file"""
        return self._summary_headers[self.snp_id]

    @property
    def effect_allele(self):
        """Key used for accessing Effect_Allele headers, groups or other attributes"""
        return "Effect_Allele"

    @property
    def sm_effect_allele(self):
        """Effect allele header index in GWAS summary file"""
        return self._summary_headers[self.effect_allele]

    @property
    def alt_allele(self):
        """Key used for accessing Alt_Allele headers, groups or other attributes"""
        return "Alt_Allele"

    @property
    def sm_alt_allele(self):
        """Alt allele header index in GWAS summary file"""
        return self._summary_headers[self.alt_allele]

    @property
    def bp_position(self):
        """Key used for accessing bp_position headers, groups or other attributes"""
        return "BP_Position"

    @property
    def sm_bp_position(self):
        """Base pair position header index in GWAS summary file"""
        return self._summary_headers[self.bp_position]

    @property
    def p_value(self):
        """Key used for accessing P_Value headers, groups or other attributes"""
        return "P_Value"

    @property
    def sm_p_value(self):
        """P value position header index in GWAS summary file"""
        return self._summary_headers[self.p_value]

    @property
    def effect_size(self):
        """Key used for accessing Effect_size headers, groups or other attributes"""
        return "Effect_size"

    @property
    def sm_effect_size(self):
        """Effect size header index in GWAS summary file"""
        return self._summary_headers[self.effect_size]

    @property
    def minor_allele_freq(self):
        """Key used for accessing Minor_Allele_Freq headers, groups or other attributes"""
        return "Minor_Allele_Freq"

    @property
    def sm_minor_allele_freq(self):
        """Minor allele Frequency header index in GWAS summary file"""
        return self._summary_headers[self.minor_allele_freq]

    @property
    def info(self):
        """Key used for accessing Minor_Allele_Freq headers, groups or other attributes"""
        return "Info"

    @property
    def sm_info(self):
        """Info score header index in GWAS summary file"""
        return self._summary_headers[self.info]

    @property
    def standard_errors(self):
        """Key used for accessing Minor_Allele_Freq headers, groups or other attributes"""
        return "Standard_Errors"

    @property
    def sm_standard_errors(self):
        """Standard error header index in GWAS summary file"""
        return self._summary_headers[self.standard_errors]

    @property
    def case_freq(self):
        """Key used for accessing Case_Freq headers, groups or other attributes"""
        return "Case_Freq"

    @property
    def sm_case_freq(self):
        """Case_Freq header index in GWAS summary file"""
        return self._summary_headers[self.case_freq]

    @property
    def case_n(self):
        """Key used for accessing Case_N headers, groups or other attributes"""
        return "Case_N"

    @property
    def sm_case_n(self):
        """Case_N header index in GWAS summary file"""
        return self._summary_headers[self.case_freq]

    @property
    def control_freq(self):
        """Key used for accessing Control_Freq headers, groups or other attributes"""
        return "Control_Freq"

    @property
    def sm_control_freq(self):
        """Control_Freq header index in GWAS summary file"""
        return self._summary_headers[self.case_freq]

    @property
    def control_n(self):
        """Key used for accessing Control_N headers, groups or other attributes"""
        return "Control_N"

    @property
    def sm_control_n(self):
        """Control_N header index in GWAS summary file"""
        return self._summary_headers[self.case_freq]

    @property
    def sm_lines(self):
        """Key for accessing lines in the summary stats file"""
        return "SM_Lines"

    @property
    def sm_variants(self):
        """Key for accessing variants associated with the snps found in the summary stats file"""
        return "SM_Variants"

    @property
    def beta(self):
        """Key used for accessing Beta headers, groups or other attributes"""
        return "Beta"

    @property
    def log_odds(self):
        """Key used for accessing Log_Odds headers, groups or other attributes"""
        return "Log_Odds"

    @property
    def nucleotide(self):
        """Key used for accessing Nucleotide headers, groups or other attributes"""
        return "Nucleotide"

    @property
    def frequency(self):
        """Key used for accessing Frequency headers, groups or other attributes"""
        return "Frequency"
