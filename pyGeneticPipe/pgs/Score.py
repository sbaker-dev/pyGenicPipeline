from pyGeneticPipe.geneticParsers.plinkObject import PlinkObject
from pyGeneticPipe.utils import error_codes as ec
from pyGeneticPipe.utils import misc as mc
from pyGeneticPipe.core.Input import Input
from csvObject import CsvObject
import numpy as np


class Score(Input):
    def __init__(self, args):
        super().__init__(args)

    def construct_pgs(self, sm_dict, core, load_path):

        sc_dict = self._construct_phenotype_dict(core, load_path)
        print(sm_dict.keys())
        print(sc_dict.keys())

        # Create a list of length number of individuals
        individual_prs = np.zeros(len(sc_dict[self.iid]))

        # Raw snps
        raw_snps = self.isolate_raw_snps(core, sm_dict)

        print(raw_snps[0])

        # todo multiple the 'beta' by the raw snp [(0, 1, 2) * -1]


        return

    def _construct_phenotype_dict(self, gen_file, load_path):
        # Load the genetic embedded information
        sc_dict = self._genetic_phenotypes(gen_file, load_path)

        # Extract the FIDs and IIDs from the sample
        if self.covariates:
            cov = CsvObject(self.covariates, set_columns=True)
            headers = {h.lower(): i for i, h in enumerate(cov.headers)}

            # Load sex if stored in covariates
            if self.sex.lower() in headers.keys():
                sc_dict[self.sex] = cov.column_data[headers[self.sex.lower()]]

            # Load PCs if they exist in the file
            pcs = [headers[h] for h in headers.keys() if h[:2] == self.pc.lower()]
            if len(pcs) > 0:
                sc_dict[self.pc] = np.array([[row[i] for i in pcs] for row in cov.row_data])
            else:
                sc_dict[self.pc] = np.array([-1 for _ in range(cov.column_length)])

        return sc_dict

    def _genetic_phenotypes(self, gen_file, load_path):
        """
        Load the full genetic data for this chromosome and isolate any information that can be isolated from it. In this
        case, .bed load types can access more than bgen due to the ability to extract sex from the .fam file.
        """

        sc_dict = {}
        # For plink files, load the fam file then extract the fid, iid and sex information
        if self.gen_type == ".bed":
            sc_dict[self.fam] = np.array(PlinkObject(load_path).get_family_identifiers())
            sc_dict[self.fid] = mc.variant_array(self.fid.lower(), sc_dict[self.fam])
            sc_dict[self.iid] = mc.variant_array(self.iid.lower(), sc_dict[self.fam])
            sc_dict[self.sex] = mc.variant_array(self.sex.lower(), sc_dict[self.fam])

        # Bgen doesn't have a fam equivalent, so just load the fid and iid
        elif self.gen_type == ".bgen":
            ids = gen_file.iid
            sc_dict[self.fid] = [fid for fid, iid in ids]
            sc_dict[self.iid] = [iid for fid, iid in ids]
            sc_dict[self.sex] = [-1 for _ in range(len(ids))]

        else:
            raise Exception("Unknown load type set")

        return sc_dict

    def _assert_construct_pgs(self):
        assert self.ld_radius, ec.missing_arg(self.operation, "LD_Radius")
