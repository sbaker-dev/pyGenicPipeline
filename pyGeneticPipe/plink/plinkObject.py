"""
This is a modified version of pandas plink found at https://github.com/limix/pandas-plink to act more like plinkio
available at https://github.com/mfranberg/libplinkio
"""
from pathlib import Path
from pyGeneticPipe.utils import error_codes as ec
from pyGeneticPipe.plink.rowObjects import BimObject


class PlinkObject:
    def __init__(self, plink_path):
        self.base_path = plink_path
        self._bed_path, self._bim_path, self._fam_path = self.validate_paths()

    def bim_object(self):
        """
        Turns the bim file into a set of objects to be called.

        :return: A BimObject for each line in the bim file
        """
        with open(self._bim_path) as f:
            return [BimObject(line) for line in f]

    def validate_paths(self):
        """
        Users may submit a path to a specific file within plink, such as a .bed/.bim/.fam or they just provide the root
        name. This method validates and returns the paths.

        :return: The path to the bed, bim, and fam file in that order
        """
        # Construct path as an object
        ld_path = Path(self.base_path)

        # Check file home directory can be reached
        assert ld_path.parent.exists(), ec.path_invalid(ld_path.parent, "_set_ld_ref")

        # If a file has a plink suffix take the stem of the name otherwise just take the name
        if ld_path.suffix == (".bed" or ".bim" or ".fam"):
            bed = Path(f"{str(ld_path.parent)}/{ld_path.stem}.bed")
            bim = Path(f"{str(ld_path.parent)}/{ld_path.stem}.bim")
            fam = Path(f"{str(ld_path.parent)}/{ld_path.stem}.fam")
        else:
            bed = Path(f"{str(ld_path.parent)}/{ld_path.name}.bed")
            bim = Path(f"{str(ld_path.parent)}/{ld_path.name}.bim")
            fam = Path(f"{str(ld_path.parent)}/{ld_path.name}.fam")

        # Check the files exists then return with mode of plink, no bgen object and a bed, bim and fam file.
        assert bed.exists(), ec.path_invalid(bed, "_set_ld_ref")
        assert bim.exists(), ec.path_invalid(bim, "_set_ld_ref")
        assert fam.exists(), ec.path_invalid(fam, "_set_ld_ref")
        return bed, bim, fam
