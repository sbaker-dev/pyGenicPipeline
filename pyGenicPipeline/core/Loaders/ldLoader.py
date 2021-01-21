from .argsParser import ArgsParser

from miscSupports import load_yaml
from pathlib import Path


class LDLoader(ArgsParser):
    def __init__(self, args):
        super().__init__(args)

        # Construct Sub Directory
        self.make_sub_directory("PGS", "LDRef")
        self.ld_directory = Path(self.working_dir, "PGS", "LDRef")

        self.ld_radius = self.args["LD_Radius"]

        # Load genomic heritability calculations if set and the pre-calculated input if that was set
        self.gm = self._set_genome()
        self.herit_calculated = self.args["Heritability_Calculated"]

    def local_values(self, norm_snps, snp_index, number_of_snps):
        """
        We want to construct a window of -r + r around each a given list of normalised snps where r is the radius.
        However, the first r and last N-r of the snps will not have r number of snps before or after them so we need to
        account for this by:

        Taking the maximum of (0, i-r) so that we never get a negative index
        Taking the minimum of (n_snps, (i + radius + 1)) to ensure we never get an index out of range

        :param norm_snps: Normalised filtered snps
        :type norm_snps: np.ndarray

        :param snp_index: Index of the current snp
        :type snp_index: int

        :param number_of_snps: Snp count of the Normalised filtered snps
        :type number_of_snps: int

        :return: An array of shape snps of a maximum of 'radius' number of snps surrounding the current snp accessed via
            index.
        """
        return norm_snps[max(0, snp_index - self.ld_radius): min(number_of_snps, (snp_index + self.ld_radius + 1))]

    def window_values(self, norm_snps, number_of_snps, window):
        """
        This chunks the snps into windows

        Where sid count is large, this will most commonly return an array of (ld_radius * 2, iid_count) bar the last
        entry where their may not be (ld_radius * 2) snps left, so the window will be shorter

        :param norm_snps: Normalised filtered snps
        :type norm_snps: np.ndarray

        :param number_of_snps: Snp count of the Normalised filtered snps
        :type number_of_snps: int

        :param window: The current window starting value from a range of (0, snp_count, ld_radius * 2)
        :type window: int

        :return: A slice of the normalised snps
        :rtype: np.ndarray
        """
        return norm_snps[window: min(number_of_snps, window + self.ld_radius * 2)]

    def _set_genome(self):
        """Load the genome file if it has been produced, else return None"""
        genome_path = Path(self.working_dir, "genome_wide_config.yaml")
        if genome_path.exists():
            return load_yaml(genome_path)
        else:
            return None
