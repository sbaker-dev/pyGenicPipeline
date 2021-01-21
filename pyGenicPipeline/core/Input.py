from pyGenicPipeline.utils import errors as ec
from pyGenicPipeline.utils import misc as mc
from .Loaders import *

from miscSupports import directory_iterator
from csvObject import CsvObject
from pyGenicParser import *
from pathlib import Path
import numpy as np
import re


class Input(SummaryLoader, FilterLoader, LDLoader, CommonGenetic, ArgsParser):
    def __init__(self, args):
        super().__init__(args)

        # General operational parameters
        self.operation = self._set_current_job(self.args["Operation"])

        # Score information
        self.make_sub_directory("PGS", "Scores")
        self.scores_directory = Path(self.working_dir, "PGS", "Scores")
        self.phenotype_file = mc.validate_path(self.args["Phenotype"])
        self.covariates_file = mc.validate_path(self.args["Covariates"])



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

    def validation_chromosomes(self):
        """
        This will create a dataset of all the chromosomes that we have to work with

        :Note: Mostly used externally to aid multi-core processing

        :return: A list of valid chromosomes
        :rtype: list
        """

        valid_chromosomes = []
        for file in directory_iterator(self.gen_directory):
            if Path(self.gen_directory, file).suffix == self.gen_type:
                valid_chromosomes.append(int(re.sub(r'[\D]', "", Path(self.gen_directory, file).stem)))
        valid_chromosomes.sort()
        return valid_chromosomes

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
        variant_names = self.snp_names(sm_dict)

        # Calculate the number of chunks required, then return the variant names split on chunk size
        chunks = int(np.ceil(len(variant_names) / self.filter_iter_size))
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

    @property
    def f_std(self):
        return f"{self.filter_key}_{self.stds}"

    @property
    def f_freq(self):
        return f"{self.filter_key}_{self.freq}"
