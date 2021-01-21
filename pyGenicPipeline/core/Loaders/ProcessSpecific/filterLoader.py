from pyGenicPipeline.core.Loaders.argsParser import ArgsParser

from pathlib import Path


class FilterLoader(ArgsParser):
    def __init__(self, args):
        super().__init__(args)

        # Set the filter directory
        self.make_sub_directory("PGS", "Filtered")
        self.filter_directory = Path(self.working_dir, "PGS", "Filtered")

        # Load the iteration size, to chunk our list to be memory concise
        self.filter_iter_size = self.args["Filter_Iter_Size"]

        # Load the Filtering parameters
        self.maf_min = self.args["Filter_Min_Maf"]
        self.freq_discrepancy = self.args["Filter_Max_Freq_Discrepancy"]

    @property
    def f_std(self):
        """Filtering uses freqs from summary, this key is the distinction for conistancy within the filter dict"""

        return f"{self.filter_key}_{self.stds}"

    @property
    def f_freq(self):
        """Filtering uses freqs from summary, this key is the distinction between that and freqs from raw_snps"""
        return f"{self.filter_key}_{self.freq}"