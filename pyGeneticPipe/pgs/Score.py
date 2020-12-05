from pyGeneticPipe.geneticParsers.plinkObject import PlinkObject
from pyGeneticPipe.utils import error_codes as ec
from pyGeneticPipe.utils import misc as mc
from pyGeneticPipe.core.Input import Input
import numpy as np


class Score(Input):
    def __init__(self, args):
        super().__init__(args)

    def construct_pgs(self, sm_dict, chromosome):

        phenotype_map = self._construct_phenotype_dict(chromosome)
        return

    def _construct_phenotype_dict(self, chromosome):
        # Load the genetic embedded information
        core, sc_dict = self._genetic_phenotypes(chromosome)

        print(sc_dict.keys())

        # Extract the FIDs and IIDs from the sample
        if self.covariates:
            # todo check for sex /
            raise NotImplementedError("Covariates not yet supported")

    def _genetic_phenotypes(self, chromosome):
        """
        Load the full genetic data for this chromosome and isolate any information that can be isolated from it. In this
        case, .bed load types can access more than bgen due to the ability to extract sex from the .fam file.
        """
        load_path = str(self.select_file_on_chromosome(chromosome, self.gen_directory, self.gen_type))

        # load The full sample, and use it to parse in any identifies that exist internally
        core = self.gen_reference(load_path)

        sc_dict = {}
        # For plink files, load the fam file then extract the fid, iid and sex information
        if self.gen_type == ".bed":
            sc_dict[self.fam] = np.array(PlinkObject(load_path).get_family_identifiers())
            sc_dict[self.fid] = mc.variant_array(self.fid.lower(), sc_dict[self.fam])
            sc_dict[self.iid] = mc.variant_array(self.iid.lower(), sc_dict[self.fam])
            sc_dict[self.sex] = mc.variant_array(self.sex.lower(), sc_dict[self.fam])

        # Bgen doesn't have a fam equivalent, so just load the fid and iid
        elif self.gen_type == ".bgen":
            ids = core.iid
            sc_dict[self.fid] = [fid for fid, iid in ids]
            sc_dict[self.iid] = [iid for fid, iid in ids]
            sc_dict[self.sex] = [-1 for _ in range(len(ids))]

        else:
            raise Exception("Unknown load type set")

        return core, sc_dict

    def _assert_construct_pgs(self):
        assert self.ld_radius, ec.missing_arg(self.operation, "LD_Radius")
