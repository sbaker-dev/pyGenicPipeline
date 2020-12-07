from pyGeneticPipe.geneticParsers.plinkObject import PlinkObject
from pyGeneticPipe.utils import error_codes as ec
from pyGeneticPipe.utils import misc as mc
from pyGeneticPipe.core.Input import Input
from csvObject import CsvObject
import numpy as np
import time


class Score(Input):
    def __init__(self, args):
        super().__init__(args)

    def construct_pgs(self, sm_dict, core, load_path):
        """
        This will construct the pgs from the weights construct with the Gibbs, the infinitesimal or gibbs estimated
        outcomes.
        """

        # Validation we have the necessary information for the scores
        self._assert_construct_pgs()

        # Construct a dict of arrays of our phenotype information
        ph_dict = self._construct_phenotype_dict(core, load_path)

        # Load the raw snps that have been isolated Raw snps
        raw_snps = self.isolate_raw_snps(core, sm_dict)

        # Calculate scores for the infinitesimal model
        self._calculate_score(sm_dict, ph_dict, self.inf_dec, raw_snps)

        # Now do the same for each model calculated by Gibbs
        for variant_fraction in sm_dict.keys():
            if self.gibbs in variant_fraction:
                self._calculate_score(sm_dict, ph_dict, variant_fraction, raw_snps)

        # Filter out individuals without defined sex, phenotypes or other invalidator information
        self._filter_ids(ph_dict)

        return

    def _construct_phenotype_dict(self, gen_file, load_path):
        """
        This will construct the phenotype dict that we will right out for the end user, storing individual level data to
        help us validate and clean individuals whom we do not have sufficient information to transfer the score too.
        """
        # Load the genetic embedded information
        ph_dict = self._genetic_phenotypes(gen_file, load_path)

        # Extract the FIDs and IIDs from the sample
        if self.covariates_file:
            cov = CsvObject(self.covariates_file, set_columns=True)
            headers = {h.lower(): i for i, h in enumerate(cov.headers)}

            # Load sex if stored in covariates
            if self.sex.lower() in headers.keys():
                ph_dict[self.sex] = cov.column_data[headers[self.sex.lower()]]

            # Load PCs if they exist in the file
            pcs = [headers[h] for h in headers.keys() if h[:2] == self.pc.lower()]
            if len(pcs) > 0:
                ph_dict[self.pc] = np.array([[row[i] for i in pcs] for row in cov.row_data])
            else:
                ph_dict[self.pc] = np.array([-1 for _ in range(cov.column_length)])

        return ph_dict

    def _genetic_phenotypes(self, gen_file, load_path):
        """
        Load the full genetic data for this chromosome and isolate any information that can be isolated from it. In this
        case, .bed load types can access more than bgen due to the ability to extract sex from the .fam file.
        """

        ph_dict = {}
        # For plink files, load the fam file then extract the fid, iid and sex information
        if self.gen_type == ".bed":
            ph_dict[self.fam] = np.array(PlinkObject(load_path).get_family_identifiers())
            ph_dict[self.fid] = mc.variant_array(self.fid.lower(), ph_dict[self.fam])
            ph_dict[self.iid] = mc.variant_array(self.iid.lower(), ph_dict[self.fam])
            ph_dict[self.sex] = mc.variant_array(self.sex.lower(), ph_dict[self.fam])

        # Bgen doesn't have a fam equivalent, so just load the fid and iid
        elif self.gen_type == ".bgen":
            ids = gen_file.iid
            ph_dict[self.fid] = np.array([fid for fid, iid in ids])
            ph_dict[self.iid] = np.array([iid for fid, iid in ids])
            ph_dict[self.sex] = np.array([-1 for _ in range(len(ids))])

        else:
            raise Exception("Unknown load type set")

        return ph_dict

    @staticmethod
    def _calculate_score(sm_dict, ph_dict, key, raw_snps):
        """
        Here we load the weights calculated and re-structure the 1 dimensional list to be a vector array. We then use
        this vector array of beta values from weights alongside the raw snps to calculate the effect of each snp based
        on the nucleotide of the individuals (0, 1 or 2) to compute a score for this key.
        """
        # Restructure weights to be a vector array
        weights = sm_dict[key]
        weights.shape = (len(sm_dict[key]), 1)

        # Calculate the PRS for the individuals
        ph_dict[key] = np.array([np.sum(row) for row in ((-1 * raw_snps) * weights).T])

    def _filter_ids(self, ph_dict):
        # Load the phenotype information
        phenotype = self._load_phenotype(ph_dict)

        if self.sex in ph_dict.keys():
            ph_dict[self.sex] = ph_dict[self.sex].astype(int)
            sex_filter = np.array([True if s != 0 else False for s in ph_dict[self.sex].astype(int)])
            mc.filter_array(ph_dict, sex_filter)

    def _load_phenotype(self, ph_dict):

        ids = []
        phenotypes = []
        with mc.open_setter(self.phenotype_file)(self.phenotype_file) as file:
            # Determine if we have custom headers or not via _loaded_sum_headers
            headers = [h.lower() for h in file.readline().split()]
            assert self.iid.lower() in headers, "FAILED TO FIND IID"

            id_index = headers.index(self.iid.lower())

            for line in file:
                line = line.split()
                ids.append(line[id_index])
                phenotypes.append(line[id_index + 1])

            file.close()

        phenotype_filter = [True if i in ids else False for i in ph_dict[self.iid]]
        mc.filter_array(ph_dict, phenotype_filter)

        return np.array(phenotypes)

    def _assert_construct_pgs(self):
        """Assert that the information required to run is present"""
        assert self.ld_radius, ec.missing_arg(self.operation, "LD_Radius")
        assert self.phenotype, ec.missing_arg(self.operation, "Phenotype")

        return time.time()
