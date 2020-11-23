"""
This is a modified version of pandas plink found at https://github.com/limix/pandas-plink to act more like plinkio
available at https://github.com/mfranberg/libplinkio
"""
from pyGeneticPipe.geneticParsers.plink.bimObject import BimObject
from pyGeneticPipe.geneticParsers.plink.bedObject import BedObject
from pyGeneticPipe.utils import error_codes as ec
from pathlib import Path


class PlinkObject:
    def __init__(self, plink_path):
        self.base_path = plink_path
        self.bed_path, self.bim_path, self.fam_path = self.validate_paths()

    def bed_object(self):
        """Bed Object"""
        print("Starting bed")
        variant_number, sample_number = self.get_dimensions()
        return BedObject(self.bed_path, variant_number, sample_number)

    def bim_object(self):
        """Bim Object"""
        return BimObject(self.bim_path)

    def bim_index(self):
        """Index of Bim files"""
        return BimObject(self.bim_path).construct_index()

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
        if (ld_path.suffix == ".bed") or (ld_path.suffix == ".bim") or (ld_path == ".fam"):
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

    def get_dimensions(self):
        """
        When working with Bed files we need to know the number of variants, the length of the bim file, and the number
        of individuals, the length of the fam file.

        :return: The length of the variants and number of samples
        """
        # Get the number of variants in the bim file
        with open(self.bim_path, "r") as f:
            variant_length = len(list(f))
        f.close()

        # Get the number of samples in the fam file
        with open(self.fam_path, "r") as f:
            number_of_samples = len(list(f))
        f.close()
        return variant_length, number_of_samples
