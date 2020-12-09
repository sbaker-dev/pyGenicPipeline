from pyGenicPipeline.pgs.SummaryCleaner import SummaryCleaner
from pyGenicPipeline.support.ShellMaker import ShellMaker
from pyGenicPipeline.pgs.FilterSnps import FilterSnps
from pyGenicPipeline.utils.misc import terminal_time
from pyGenicPipeline.pgs.LDHerit import LDHerit
from pyGenicPipeline.core.Input import Input
from pyGenicPipeline.pgs.Gibbs import Gibbs
from pyGenicPipeline.pgs.Score import Score
from csvObject import write_csv
from colorama import init, Fore
import numpy as np
import time


class Main(ShellMaker, SummaryCleaner, FilterSnps, LDHerit, Gibbs, Score, Input):
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
        validation, ref = self.construct_validation(load_path)

        # Then we need to take these samples to construct valid snps, these snps are extract for this chromosome from
        # our summary stats, and then cleaned for possible errors.
        sm_dict = self.clean_summary_statistics(chromosome, load_path, validation, ref)

        # Filter our genetic types for snps, such as those that have undesirable frequencies.
        sm_dict = self.filter_snps(self.val_prefix, validation, sm_dict, chromosome)
        sm_dict = self.filter_snps(self.ref_prefix, ref, sm_dict, chromosome)

        # Compute the chromosome specific ld scores and heritability
        self.compute_ld_scores(sm_dict, len(sm_dict[self.sm_variants]), ref.iid_count)

        # Construct rows to right out
        rows_out = []
        for v, log_odds, beta, freq, ld in zip(sm_dict[self.sm_variants], sm_dict[self.log_odds], sm_dict[self.beta],
                                               sm_dict[self.freq], sm_dict[self.ld_scores]):
            rows_out.append(v.items() + [log_odds, beta, freq, ld])

        write_csv(self.clean_directory, f"Cleaned_{chromosome}", self.clean_headers, rows_out)
        print(Fore.LIGHTCYAN_EX + f"Finished {self.operation} for chromosome {chromosome} at {terminal_time()}.\n"
                                  f"Total time spent was {round(time.time() - start_time, 2)} Seconds\n")

    def pgs_chromosome_scores(self, chromosome):

        # Assert we have the genome file form genome_wide_heritability, set dict to of this chromosome and genome
        assert self.gm, "missing g"
        self.gm = {**self.gm[chromosome], **self.gm[self.genome_key]}

        # Construct the dict of values we need for this run from our cleaned data
        load_path = str(self.select_file_on_chromosome(chromosome, self.clean_directory, ".csv"))
        sm_dict = self.sm_dict_from_csv(load_path)

        # Load the genetic reference as was in the first stage
        load_path = str(self.select_file_on_chromosome(chromosome, self.gen_directory, self.gen_type))
        _, ref = self.construct_validation(load_path)

        # Create the normalised snps
        sm_dict = self.filter_snps(self.ref_prefix, ref, sm_dict, chromosome)

        # Compute the ld scores and dict
        self.compute_ld_scores(sm_dict, self.gm[self.count_snp], self.gm[self.count_iid], ld_dict=True)

        # Construct the weighted betas for each snp for use in construction of scores
        self.construct_gibbs_weights(sm_dict, chromosome)

        # Construct the Poly-genetic Scores
        self.construct_chromosome_pgs(sm_dict, load_path, chromosome)

    def debug_cross_check(self):
        """
        This is designed to allow us to cross check with LDPred, following the data from:
        https://choishingwan.github.io/PRS-Tutorial/
        """
        # Load the validation and core samples, as well as the indexer
        chromosome = 1
        self.gm = {**self.gm[chromosome], **self.gm[self.genome_key]}
        load_path = str(self.select_file_on_chromosome(chromosome, self.gen_directory, self.gen_type))
        ref = self.gen_reference(load_path)
        validation = self.gen_reference(load_path)

        # Then we need to take these samples to construct valid snps, these snps are extract for this chromosome from
        # our summary stats, and then cleaned for possible errors.
        sm_dict = self.clean_summary_statistics(chromosome, load_path, validation, ref)

        # Filter our genetic types for snps, such as those that have undesirable frequencies.
        # sm_dict = self.filter_snps(self.val_prefix, validation, sm_dict, chromosome)
        sm_dict = self.filter_snps(self.ref_prefix, ref, sm_dict, chromosome)

        # Compute the chromosome specific ld scores and heritability
        self.compute_ld_scores(sm_dict, len(sm_dict[self.sm_variants]), ref.iid_count, ld_dict=True)

        # Mirror test environment gm
        self.gm[self.count_snp] = 5693
        self.gm[self.count_iid] = 483
        self.gm[self.herit] = 0.04553305821357676
        self.gibbs_causal_fractions = [1]

        # This is the start of scores
        # Construct the Weight
        self.construct_gibbs_weights(sm_dict, chromosome)

        # 0.11157818410953191 ish
        print(sm_dict[self.inf_dec])
        print(np.sum(sm_dict[self.inf_dec]))

        # 0.21582699762327068 ish
        print(sm_dict[f"{self.gibbs}_{self.gibbs_causal_fractions[0]}"])
        print(np.sum(sm_dict[f"{self.gibbs}_{self.gibbs_causal_fractions[0]}"]))

        # Construct the Poly-genetic Scores
        self.construct_chromosome_pgs(sm_dict, load_path, chromosome)
