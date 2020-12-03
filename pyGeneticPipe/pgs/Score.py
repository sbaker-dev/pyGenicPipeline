from pyGeneticPipe.utils import error_codes as ec
from pyGeneticPipe.utils import misc as mc
from pyGeneticPipe.core.Input import Input
import numpy as np


class Score(Input):
    def __init__(self, args):
        super().__init__(args)

    def construct_pgs(self):

        phenotype_map = self._construct_phenotype_dict()
        return

    def _construct_phenotype_dict(self):
        # Isolate a file within our sample
        valid_chromosomes = self.validation_chromosomes()
        load_path = str(self.select_file_on_chromosome(valid_chromosomes[0]))

        # Construct the validation and reference sample
        # todo: we want to construct a score for everyone, but based on reference validation sample of individuals?
        # todo: so here, we would want a separate method to load all of the samples?
        validation, core = self.construct_validation(load_path)

        # Extract the FIDs and IIDs from the sample
        if self.covariates:
            raise NotImplementedError("Covariates not yet supported")

        return core.iid

    def _assert_construct_pgs(self):
        assert self.ld_radius, ec.missing_arg(self.operation, "LD_Radius")
