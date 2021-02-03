from pyGenicPipeline.utils.misc import validate_path, set_current_job
from .Loaders import *

from miscSupports import directory_iterator
from csvObject import CsvObject
from pyGenicParser import *
from pathlib import Path
import numpy as np
import re


class Input(SummaryLoader, FilterLoader, LDLoader, GibbLoader, CommonGenetic, ArgsParser):
    def __init__(self, args):
        """
        This Class Loads all other Loaders that represent specific or common attributes/methods that are used across
        processors.

        Methods or attributes set within this class are dependent on multiple attributes from its inherited sub class or
        are specific to Input

        :param args: Args Passed by the users
        """
        super().__init__(args)

        # General operational parameters
        self.operation = set_current_job(self.args["Operation"])
        self.iter_size = self.args["Iter_Size"]

        # Score information
        self.make_sub_directory("PGS", "Scores")
        self.scores_directory = Path(self.working_dir, "PGS", "Scores")
        # todo expose this but have a default
        self.score_type = "INF_Weights"
        self.phenotype_file = validate_path(self.args["Phenotype"])
        self.covariates_file = validate_path(self.args["Covariates"])

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
        individuals in the data set. This chunk snp_names by a dimension calculated from the number of variant names by
        the snp_iteration set by the user.

        :param sm_dict: Dict of cleaned summary statistics or a list of snp names
        :type sm_dict: dict | list

        :param chunk_return: If True returns the chunk size as well as the variants, otherwise just the variants
        :type chunk_return: bool

        :return:  np.array_split numpy arrays of snp names, or this combined with the chunk size
        """

        # Determine the snp names base on the type provided to sm_dict
        if isinstance(sm_dict, dict):
            # Extract variant names
            variant_names = self.snp_names(sm_dict)
        elif isinstance(sm_dict, list):
            variant_names = sm_dict
        else:
            raise TypeError(f"Chunk snps expected sm_dict or a list of snp names but found {type(sm_dict)}")

        # Calculate the number of chunks required, then return the variant names split on chunk size
        chunks = int(np.ceil(len(variant_names) / self.iter_size))
        if chunk_return:
            return np.array_split(variant_names, chunks), chunks
        else:
            return np.array_split(variant_names, chunks)

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
