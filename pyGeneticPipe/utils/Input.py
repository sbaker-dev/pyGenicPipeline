from pyGeneticPipe.utils import error_codes as ec
from pathlib import Path


class Input:
    def __init__(self, args):
        self.debug = args["Debug"]
        self.ld_ref_mode, self.bgen, self.bed, self.bim, self.fam = self._set_ld_ref(args["LD_Reference_Genotype"])

    @staticmethod
    def _set_ld_ref(ref_path):
        # Construct path as an object
        ld_path = Path(ref_path)

        # Check file directory can be reached
        assert ld_path.parent.exists(), ec.path_invalid(ld_path.parent, "_set_ld_ref")

        # Check the mode we are running in, and return args of mode, bgenObject and the 3 plink files accordingly
        if ld_path.suffix == ".bgen":
            # return "bgen", bgenObject(ld_path), None, None, None
            raise NotImplementedError("Reading bgen files not yet implemented")

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
        return "plink", None, bed, bim, fam
