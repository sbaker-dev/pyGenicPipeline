from pyGeneticPipe.utils.Input import Input
from pyGeneticPipe.plink.plinkObject import PlinkObject
import numpy as np
import h5py


class ProcessSummary(Input):
    def __init__(self, args):
        super().__init__(args)

        self._bim_by_chromosomes = self._set_bim_by_chromosomes()
        # self._validation_bim_by_chromosomes = self._set_validation_file()

    def pre_process_plink(self, hdf5_path):

        # Open the cleaned summary statistics
        summary_file = h5py.File(hdf5_path, "r")["Sum_Stats"]

        for chrom_bim in self._bim_by_chromosomes:

            # Extract the summary statistics data for this chromosome, if we fail to extract it continue
            chrom_summary = self._summary_data_by_chromosome(summary_file, chrom_bim)
            if not chrom_summary:
                continue

            print(chrom_summary)
            break

    def pre_process_bgen(self):
        raise NotImplementedError("Bgen files not yet implemented")

    @staticmethod
    def _summary_data_by_chromosome(summary_file, chrom_bim):
        """
        Try to extract the data by chromosome, if we failed to extract it summary statistics file return None so it will
        not crash and warn the user this has happened.
        """
        try:
            return summary_file[chrom_bim.chromosome]
        except KeyError:
            print(f"WARNING: Failed to find {chrom_bim.chromosome} in summary statistics file")
            return None

    def _set_bim_by_chromosomes(self):
        """
        Extract the information in the bim file but by chromosomes
        """
        # todo: will need to be generalised for validation
        # Create a plink Objected
        plink_obj = PlinkObject(self.args["LD_Reference_Genotype"])

        # Extract the information by chromosomes
        return plink_obj.bim_by_chromosome()


