from abc import ABC

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
from pathlib import Path
import numpy as np
import time


class Main(ShellMaker, SummaryCleaner, FilterSnps, LDHerit, Gibbs, Score, Input, ABC):
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

    def pgs_clean_and_filter(self):
        """
        This will clean and coordinate the summary statistics and merge it with our genetic data via Cleaners
        clean_summary_stats.

        Then this information will be passed into the LDPred's gibbs estimator to construct a weight that can then be
        validated for convergence based on genome wide heritability and be used to construct the pgs in
        _pgs_construct_scores
        """
        if self.multi_core_splitter:
            self.chromosome_clean_and_filter(self.multi_core_splitter)

        else:
            valid_chromosomes = self.validation_chromosomes()
            for chromosome in valid_chromosomes:
                self.chromosome_clean_and_filter(chromosome)

    def pgs_clean_summary_stats(self):
        """Filtering may take a lot of computing, so we can split these operations. This is the Cleaner"""
        if self.multi_core_splitter:
            self.chromosome_clean_summary_stats(self.multi_core_splitter)
        else:
            valid_chromosomes = self.validation_chromosomes()
            for chromosome in valid_chromosomes:
                self.chromosome_clean_summary_stats(chromosome)

    def pgs_filter_cleaned(self):
        """Filtering may take a lot of computing, so we can split these operations. This is the Filter"""
        if self.multi_core_splitter:
            self.chromosome_filter(self.multi_core_splitter)
        else:
            valid_chromosomes = self.validation_chromosomes()
            for chromosome in valid_chromosomes:
                self.chromosome_filter(chromosome)

    def pgs_ld_scores(self):
        """Compute LD scores for Genome-wide regression"""
        if self.multi_core_splitter:
            self.chromosome_ld_score(self.multi_core_splitter)
        else:
            valid_chromosomes = self.validation_chromosomes()
            for chromosome in valid_chromosomes:
                self.chromosome_ld_score(chromosome)

    def pgs_weights(self):
        if self.multi_core_splitter:
            self.chromosome_weights(self.multi_core_splitter)
        else:
            valid_chromosomes = self.validation_chromosomes()
            for chromosome in valid_chromosomes:
                self.chromosome_weights(chromosome)

    def pgs_weights_and_scores(self):
        """
        This will construct pgs scores based on the cleaned data produced from pgs_clean_and_coordinate and the genome
        type file create from genome_wide_heritability.
        """
        if self.multi_core_splitter:
            self.chromosome_pgs_weights_and_scores(self.multi_core_splitter)
        else:
            valid_chromosomes = self.validation_chromosomes()
            for chromosome in valid_chromosomes:
                self.chromosome_pgs_weights_and_scores(chromosome)

    def chromosome_clean_and_filter(self, chromosome):
        """This takes the value of a current chromosome and constructs the weights described in pgs_construct_weights"""
        # Load the validation and core samples, as well as the indexer
        start_time = time.time()
        load_path = str(self.select_file_on_chromosome(chromosome, self.gen_directory, self.gen_type))
        validation, ref = self.construct_validation(load_path)

        # Then we need to take these samples to construct valid snps, these snps are extract for this chromosome from
        # our summary stats, and then cleaned for possible errors.
        sm_dict = self.clean_summary_statistics(chromosome, load_path, validation, ref)
        self.check_sm_dict(sm_dict)

        # Filter our genetic types for snps, such as those that have undesirable frequencies.
        self.filter_snps(sm_dict, ref, chromosome, validation)
        print(f"Found {len(sm_dict[self.sm_variants])} Snps after filtering")

        # Compute the chromosome specific ld scores and heritability
        self.compute_ld_scores(sm_dict, ref, len(sm_dict[self.sm_variants]), ref.iid_count, load_path)
        self._write_full_cleaned(sm_dict, chromosome, start_time)

    def chromosome_pgs_weights_and_scores(self, chromosome):

        # Assert we have the genome file form genome_wide_heritability, set dict to of this chromosome and genome
        assert self.gm, "missing g"
        self.gm = {**self.gm[chromosome], **self.gm[self.genome_key]}

        # Construct the dict of values we need for this run from our cleaned data
        load_path = str(self.select_file_on_chromosome(chromosome, self.clean_directory, ".csv"))
        sm_dict = self.sm_dict_from_csv(load_path)

        # Load the genetic reference as was in the first stage
        load_path = str(self.select_file_on_chromosome(chromosome, self.gen_directory, self.gen_type))
        _, ref = self.construct_validation(load_path)

        # Compute the ld scores and dict
        if self.gibbs:
            self.compute_ld_scores(sm_dict, ref, self.gm[self.count_snp], self.gm[self.count_iid], load_path,
                                   ld_dict=True)

        # Construct the weighted betas for each snp for use in construction of scores
        self.construct_gibbs_weights(sm_dict, load_path, chromosome)

        # Construct the Poly-genetic Scores
        self.construct_chromosome_pgs(sm_dict, load_path, chromosome)

    def chromosome_clean_summary_stats(self, chromosome):
        """Just the summary stats cleaneer from clean and filter"""
        # Load the validation and core samples, as well as the indexer
        start_time = time.time()
        load_path = str(self.select_file_on_chromosome(chromosome, self.gen_directory, self.gen_type))
        validation, ref = self.construct_validation(load_path)

        # Then we need to take these samples to construct valid snps, these snps are extract for this chromosome from
        # our summary stats, and then cleaned for possible errors.
        sm_dict = self.clean_summary_statistics(chromosome, load_path, validation, ref)
        self.check_sm_dict(sm_dict)
        print(f"Found {len(sm_dict[self.sm_variants])}, will create {len(self.chunked_snp_names(sm_dict)[0])} Chunks")

        rows_out = []
        for v, log_odds, beta, freq, in zip(sm_dict[self.sm_variants], sm_dict[self.log_odds], sm_dict[self.beta],
                                            sm_dict[self.freq]):
            rows_out.append(v.items() + [log_odds, beta, freq])

        write_csv(self.clean_directory, f"Cleaned_{chromosome}", self.clean_headers[:-1], rows_out)
        print(Fore.LIGHTCYAN_EX + f"Finished {self.operation} for chromosome {chromosome} at {terminal_time()}.\n"
                                  f"Total time spent was {round(time.time() - start_time, 2)} Seconds\n")

    def chromosome_filter(self, chromosome):
        # Construct the dict of values from our saved csv
        load_path = str(self.select_file_on_chromosome(chromosome, self.clean_directory, ".csv"))
        sm_dict = self.sm_dict_from_csv(load_path)

        # Chunk the data that was saved
        snp_list, freqs, bp_positions = self.chunked_snp_names(sm_dict)

        # Filter based on the current filter_index
        filter_dict = {self.snp_id: snp_list[self.filter_index], self.freq: freqs[self.filter_index],
                       self.bp_position: bp_positions[self.filter_index],
                       self.filter_key: np.full(len(snp_list[self.filter_index]), True)}

        # Now load the genetic files
        load_path = str(self.select_file_on_chromosome(chromosome, self.gen_directory, self.gen_type))
        validation, ref = self.construct_validation(load_path)

        # Filter on ref and validation
        filter_dict = self.filter_snp_chunk(ref, filter_dict, chromosome)
        filter_dict = self.filter_snp_chunk(validation, filter_dict, chromosome)

        # Save the filter array
        np.save(Path(self.filter_directory, f"Filtered_{self.filter_index}_{chromosome}"), filter_dict[self.filter_key])

    def chromosome_ld_score(self, chromosome):
        # Construct the dict of values we need for this run from our cleaned data
        start_time = time.time()
        load_path = str(self.select_file_on_chromosome(chromosome, self.clean_directory, ".csv"))
        sm_dict = self.sm_dict_from_csv(load_path)

        # Load the genetic reference as was in the first stage
        load_path = str(self.select_file_on_chromosome(chromosome, self.gen_directory, self.gen_type))
        _, ref = self.construct_validation(load_path)

        # Compute the ld scores and dict
        self.compute_ld_scores(sm_dict, ref, len(sm_dict[self.sm_variants]), ref.iid_count, load_path)
        self._write_full_cleaned(sm_dict, chromosome, start_time)

    def chromosome_weights(self, chromosome):
        assert self.gm, "missing g"
        self.gm = {**self.gm[chromosome], **self.gm[self.genome_key]}

        # Construct the dict of values we need for this run from our cleaned data
        load_path = str(self.select_file_on_chromosome(chromosome, self.clean_directory, ".csv"))
        sm_dict = self.sm_dict_from_csv(load_path)

        # Load the genetic reference as was in the first stage
        load_path = str(self.select_file_on_chromosome(chromosome, self.gen_directory, self.gen_type))
        _, ref = self.construct_validation(load_path)

        # Compute the ld scores and dict
        if self.gibbs_run:
            self.compute_ld_scores(sm_dict, ref, self.gm[self.count_snp], self.gm[self.count_iid], load_path,
                                   ld_dict=True)

        # Construct the weighted betas for each snp for use in construction of scores
        self.construct_gibbs_weights(sm_dict, load_path, chromosome)

        # Save the keys that have weights to a file
        columns = {key: v for key, v in sm_dict.items() if self.gibbs in key or key == self.inf_dec}
        headers = list(columns.keys())
        values = np.array([v for v in columns.values()]).T.flatten()

        write_csv(self.weights_directory, f"Weights_{chromosome}", headers, values)
        print(f"Constructed weights file for {chromosome}")

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
        self.check_sm_dict(sm_dict)

        # Filter our genetic types for snps, such as those that have undesirable frequencies.
        # sm_dict = self.filter_snps(self.val_prefix, validation, sm_dict, chromosome)
        self.filter_snps(sm_dict, ref, chromosome)

        # Compute the chromosome specific ld scores and heritability
        self.compute_ld_scores(sm_dict, ref, len(sm_dict[self.sm_variants]), ref.iid_count, load_path, ld_dict=True)

        # Mirror test environment gm
        self.gm[self.count_snp] = 5693
        self.gm[self.count_iid] = 483
        self.gm[self.herit] = 0.04553305821357676
        self.gibbs_causal_fractions = [1]
        # self.gibbs_run = True

        # This is the start of scores
        # Construct the Weight
        self.construct_gibbs_weights(sm_dict, load_path, chromosome)

        print(sm_dict[self.inf_dec])
        print(np.sum(sm_dict[self.inf_dec]))
        print(f"0.11157818410953191 is the target (ish)\n")

        if self.gibbs_run:
            print(sm_dict[f"{self.gibbs}_{self.gibbs_causal_fractions[0]}"])
            print(np.sum(sm_dict[f"{self.gibbs}_{self.gibbs_causal_fractions[0]}"]))
            print(f"0.21582699762327068 is that target (ish)\n")

        # Construct the Poly-genetic Scores
        self.construct_chromosome_pgs(sm_dict, load_path, chromosome)

    def _write_full_cleaned(self, sm_dict, chromosome, start_time):
        # Construct rows to right out
        rows_out = []
        for v, log_odds, beta, freq, ld in zip(sm_dict[self.sm_variants], sm_dict[self.log_odds], sm_dict[self.beta],
                                               sm_dict[self.freq], sm_dict[self.ld_scores]):
            rows_out.append(v.items() + [log_odds, beta, freq, ld])

        write_csv(self.clean_directory, f"Cleaned_{chromosome}", self.clean_headers, rows_out)
        print(Fore.LIGHTCYAN_EX + f"Finished {self.operation} for chromosome {chromosome} at {terminal_time()}.\n"
                                  f"Total time spent was {round(time.time() - start_time, 2)} Seconds\n")
