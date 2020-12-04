from pyGeneticPipe.pgs.SummaryCleaner import SummaryCleaner
from pyGeneticPipe.support.ShellMaker import ShellMaker
from pyGeneticPipe.pgs.FilterSnps import FilterSnps
from pyGeneticPipe.utils.misc import terminal_time
from pyGeneticPipe.pgs.LDHerit import LDHerit
from pyGeneticPipe.core.Input import Input
from csvObject import write_csv
from colorama import init, Fore
import time


class Main(ShellMaker, SummaryCleaner, FilterSnps, LDHerit, Input):
    def __init__(self, args):
        """
        This Class inherits all other classes that can be used, and then execute the job via getattr
        :param args: Json args that have been set via GUI or manually
        """
        init(autoreset=True)
        super().__init__(args)
        if self.operation is not None:
            print(f"Starting {self.operation}: {terminal_time()}")
            getattr(Main, self.operation)(self)

    def pgs_clean_and_coordinate(self):
        """
        This will clean and coordinate the summary statistics and merge it with our genetic data via Cleaners
        clean_summary_stats.

        Then this information will be passed into the LDPred's gibbs estimator to construct a weight that can then be
        validated for convergence based on genome wide heritability and be used to construct the pgs in
        _pgs_construct_scores
        """
        if self.multi_core_splitter:
            self.chromosome_clean_and_coordinate(self.multi_core_splitter)

        else:
            valid_chromosomes = self.validation_chromosomes()
            for chromosome in valid_chromosomes:
                self.chromosome_clean_and_coordinate(chromosome)

    def pgs_scores(self):
        """
        This will construct pgs scores based on the cleaned data produced from pgs_clean_and_coordinate and the genome
        type file create from genome_wide_heritability.
        """
        if self.multi_core_splitter:
            self.pgs_chromosome_scores(self.multi_core_splitter)
        else:
            valid_chromosomes = self.validation_chromosomes()
            for chromosome in valid_chromosomes:
                self.pgs_chromosome_scores(chromosome)

    def chromosome_clean_and_coordinate(self, chromosome):
        """This takes the value of a current chromosome and constructs the weights described in pgs_construct_weights"""
        # Load the validation and core samples, as well as the indexer
        start_time = time.time()
        load_path = str(self.select_file_on_chromosome(chromosome, self.gen_directory, self.gen_type))
        validation, core = self.construct_validation(load_path)

        # Then we need to take these samples to construct valid snps, these snps are extract for this chromosome from
        # our summary stats, and then cleaned for possible errors.
        sm_dict = self.clean_summary_statistics(chromosome, load_path, validation, core)

        # Filter our genetic types for snps, such as those that have undesirable frequencies.
        sm_dict = self.filter_snps(self.val_prefix, validation, sm_dict, chromosome)
        sm_dict = self.filter_snps(self.ref_prefix, core, sm_dict, chromosome)

        # Compute the chromosome specific ld scores and heritability
        self.compute_ld_scores(sm_dict, len(sm_dict[self.sm_variants]), core.iid_count)

        # Construct rows to right out
        rows_out = []
        for v, log_odds, beta, freq, ld in zip(sm_dict[self.sm_variants], sm_dict[self.log_odds], sm_dict[self.beta],
                                               sm_dict[self.freq], sm_dict[self.ld_scores]):
            rows_out.append(v.items() + [log_odds, beta, freq, ld])

        write_csv(self.clean_directory, f"Cleaned_{chromosome}", self.clean_headers, rows_out)
        print(Fore.LIGHTCYAN_EX + f"Finished {self.operation} for chromosome {chromosome} at {terminal_time()}.\n"
                                  f"Total time spent was {round(time.time() - start_time, 2)} Seconds\n")

    def pgs_chromosome_scores(self, chromosome):

        load_path = self.select_file_on_chromosome(chromosome, self.clean_directory, ".csv")

        # Construct the dict of values we need for this run from our cleaned data
        sm_dict = self.sm_dict_from_csv(load_path)

        # Load the genetic core
        load_path = self.select_file_on_chromosome(chromosome, self.gen_directory, self.gen_type)
        _, core = self.construct_validation(load_path)

        sm_dict = self.filter_snps(self.ref_prefix, core, sm_dict, chromosome)
        print(sm_dict.keys())


        # print(sm_dict)
