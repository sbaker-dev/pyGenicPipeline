from pyGenicPipeline.pgs import *
from ..Input import Input


class PGS(SummaryCleaner, FilterSnps, LDHerit, Gibbs, Score, Input):
    def __init__(self, args):
        super().__init__(args)

    def pgs_clean_summary_stats(self):
        """Filtering may take a lot of computing, so we can split these operations. This is the Cleaner"""
        pass
