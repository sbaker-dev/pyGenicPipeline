from abc import ABC

from pyGenicPipeline.Processors import *
from ..Input import Input

from miscSupports import terminal_time
from csvObject import write_csv
from colorama import Fore
import time


class PGS(SummaryCleaner, Input, ABC):
    def __init__(self, args):
        super().__init__(args)

    def pgs_clean_summary_stats(self):
        """Filtering may take a lot of computing, so we can split these operations. This is the Cleaner"""
        # Load the validation and core samples, as well as the indexer
        start_time = time.time()

        load_path = str(self.select_file_on_chromosome(self.gen_directory, self.gen_type))
        validation, ref = self.construct_validation(load_path)

        # Then we need to take these samples to construct valid snps, these snps are extract for this chromosome from
        # our summary stats, and then cleaned for possible errors.
        sm_dict = self.clean_summary_statistics(load_path, validation, ref)
        self.check_sm_dict(sm_dict)
        print(f"Found {len(sm_dict[self.sm_variants])}, will create {len(self.chunked_snp_names(sm_dict)[0])} Chunks")

        rows_out = []
        for v, log_odds, beta, freq, in zip(sm_dict[self.sm_variants], sm_dict[self.log_odds], sm_dict[self.beta],
                                            sm_dict[self.freq]):
            rows_out.append(v.items() + [log_odds, beta, freq])

        write_csv(self.clean_directory, f"Cleaned_{self.target_chromosome}", self.clean_headers[:-1], rows_out)
        print(Fore.LIGHTCYAN_EX + f"Finished {self.operation} for chromosome {self.target_chromosome} at "
                                  f"{terminal_time()}.\nTotal time spent was {round(time.time() - start_time, 2)}"
                                  f" Seconds\n")
