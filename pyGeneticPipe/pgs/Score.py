from pyGeneticPipe.core.Input import Input


class Score(Input):
    def __init__(self, args):
        super().__init__(args)

    def construct_pgs(self):
        phenotype_map = self._construct_phenotype_dict()
        print(phenotype_map)
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
