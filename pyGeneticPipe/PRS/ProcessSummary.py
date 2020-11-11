from pyGeneticPipe.utils.Input import Input
from pyGeneticPipe.plink.plinkObject import PlinkObject
import numpy as np
import h5py


class ProcessSummary(Input):
    def __init__(self, args):
        super().__init__(args)

    def pre_process_plink(self, hdf5_path):

        # Extract the loci
        bim_loci = PlinkObject(self.args["LD_Reference_Genotype"]).bim_object()

        # Extract the unique sorted chromosomes
        valid_chromosomes = np.unique([int(loci.chromosome) for loci in bim_loci])
        valid_chromosomes.sort()

        # ssf = h5py.File(hdf5_path, "r")["Sum_Stats"]
        #
        # for chromosome in valid_chromosomes:
        #     print(chromosome)
        #     chromosome_group = ssf[str(chromosome)]
        #
        #     chromosome_snp_ids = (chromosome_group['snp_id'][...]).astype('<U30')
        #
        #     # A map from sid to index for summary stats
        #     ss_sid_dict = {}
        #     for i, sid in enumerate(chromosome_snp_ids):
        #         ss_sid_dict[sid] = i
        #
        #     print(chromosome_snp_ids)
        #     print(ss_sid_dict)
        #     break

    def pre_process_bgen(self):
        raise NotImplementedError("Bgen files not yet implemented")


