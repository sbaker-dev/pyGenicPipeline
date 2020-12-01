from pyGeneticPipe.support.ShellMaker import ShellMaker
from pyGeneticPipe.utils.misc import terminal_time
from pyGeneticPipe.pgs.SummaryCleaner import SummaryCleaner
from pyGeneticPipe.pgs.FilterSnps import FilterSnps
from pyGeneticPipe.core.Input import Input
from colorama import init


class Main(ShellMaker, SummaryCleaner, FilterSnps, Input):
    def __init__(self, args):
        """
        This Class inherits all other classes that can be used, and then execute the job via getattr
        :param args: Json args that have been set via GUI or manually
        """
        init(autoreset=True)
        super().__init__(args)
        print(f"Starting {self.operation}: {terminal_time()}")
        getattr(Main, self.operation)(self)

    def _pgs_construct_weights(self):
        """
        This will clean and coordinate the summary statistics and merge it with our genetic data via Cleaners
        clean_summary_stats.

        Then this information will be passed into the LDPred's gibbs estimator to construct a weight that can then be
        validated for convergence based on genome wide heritability and be used to construct the pgs in
        _pgs_construct_scores
        """
        if self.multi_core_splitter:
            self.clean_summary_statistics(self.multi_core_splitter)

        else:
            valid_chromosomes = self._validation_chromosomes()
            for chromosome in valid_chromosomes:
                self.clean_summary_statistics(chromosome)


