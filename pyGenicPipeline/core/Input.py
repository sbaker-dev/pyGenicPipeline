from pyGenicPipeline.utils import errors as ec
from pyGenicPipeline.utils import misc as mc
from .Loaders import *

from miscSupports import load_yaml
from csvObject import CsvObject
from pyGenicParser import *
from pathlib import Path
import numpy as np


class Input(CommonGenetic, SummaryLoader, ArgsParser):
    def __init__(self, args):
        super().__init__(args)

        # General operational parameters
        self.operation = self._set_current_job(self.args["Operation"])

        # Todo this is bad, needs to be removed at some point
        self._config = self.config["Other"]

        # Set filter information
        self.make_sub_directory("PGS", "Filtered")
        self.filter_directory = Path(self.working_dir, "PGS", "Filtered")
        self.maf_min = self._config["Min_Maf"]
        self.freq_discrepancy = self._config["Max_Freq_Discrepancy"]
        self.filter_index = self.args["Filter_Index"]
        self._filter_iter_size = self.args["Filter_Range"]

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
        self.make_sub_directory("PGS", "Scores")
        self.scores_directory = Path(self.working_dir, "PGS", "Scores")
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

    def _set_causal_fractions(self):
        """
        If the user has provided a set of causal fractions of variants to use for the gibbs sampler then use those, else
        use the default that LDPred used.
        """
        if self.args["Causal_Fractions"]:
            return self.args["Causal_Fractions"]
        else:
            return [1, 0.3, 0.1, 0.03, 0.01, 0.003, 0.001]

    def sm_dict_from_csv(self, directory, name):
        """
        Load a saved cleaned file from summary statistics for a given chromosome found at the load_path provided, and
        use this to construct the sm_dict we pass between methods
        """
        load_file = CsvObject(Path(directory, name), self.cleaned_types, set_columns=True)

        chromo = load_file.column_data[self.summary_dict[self.chromosome]]
        bp_pos = load_file.column_data[self.summary_dict[self.bp_position]]
        snp_id = load_file.column_data[self.summary_dict[self.snp_id]]
        effect = load_file.column_data[self.summary_dict[self.effect_allele]]
        alt = load_file.column_data[self.summary_dict[self.alt_allele]]
        log = load_file.column_data[self.summary_dict[self.log_odds]]
        beta = load_file.column_data[self.summary_dict[self.beta]]
        freq = load_file.column_data[self.summary_dict[self.freq]]

        sm_variants = [Variant(ch, bp, sn, ef, al) for ch, bp, sn, ef, al in zip(chromo, bp_pos, snp_id, effect, alt)]
        return {self.sm_variants: np.array(sm_variants), self.log_odds: np.array(log), self.beta: np.array(beta),
                self.freq: np.array(freq)}

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

    @property
    def f_std(self):
        return f"{self.filter_key}_{self.stds}"

    @property
    def f_freq(self):
        return f"{self.filter_key}_{self.freq}"

    @property
    def c_ld_score(self):
        """LD_Score header index in Cleaned Data file"""
        return self.summary_dict[self.ld_scores]

    @property
    def c_beta(self):
        """Beta header index in Cleaned Data file"""
        return self.summary_dict[self.beta]
