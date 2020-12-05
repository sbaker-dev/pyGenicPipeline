from pyGeneticPipe.geneticParsers.variantObjects import BimVariant
from pyGeneticPipe.utils import error_codes as ec
from pathlib import Path


class PlinkObject:
    def __init__(self, genetic_path):
        self.bed_file, self.bim_file, self.fam_file = self.validate_paths(genetic_path)
        self.bim_file = open(self.bim_file, "r")

    def construct_bim_index(self):
        """
        Bim files need to be index via seek, so we can extract a given snp loci without having to store all of this of
        then in memory
        """
        indexer = {}
        cumulative_seek = 0
        for line in self.bim_file:
            chromosome, variant_id, morgan_pos, bp_position, a1, a2 = line.split()
            indexer[variant_id] = cumulative_seek
            cumulative_seek += len(line)

        self.bim_file.close()
        return indexer

    def get_variant(self, seek, as_variant=False):
        """
        Extract a given variants loci based on the seek index from construct_index
        :param seek: The seek index
        :param as_variant: If you want it as a standardised across parameter variant, or a Bim Variant with morgan pos
        :return: The line
        """
        self.bim_file.seek(seek)
        chromosome, variant_id, morgan_pos, bp_position, a1, a2 = self.bim_file.readline().split()

        if as_variant:
            return BimVariant(chromosome, variant_id, morgan_pos, bp_position, a1, a2).to_variant()
        else:
            return BimVariant(chromosome, variant_id, morgan_pos, bp_position, a1, a2)

    @staticmethod
    def validate_paths(genetic_path):
        """
        Users may submit a path to a specific file within plink, such as a .bed/.bim/.fam or they just provide the root
        name. This method validates and returns the paths.

        :return: The path to the bed, bim, and fam file in that order
        """
        # Construct path as an object
        ld_path = Path(genetic_path)

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
