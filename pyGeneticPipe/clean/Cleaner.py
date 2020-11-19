from pyGeneticPipe.utils import error_codes as ec
from pyGeneticPipe.utils import misc as mc
from pyGeneticPipe.core.Input import Input
from pysnptools.distreader import Bgen
from pysnptools.snpreader import Bed
from pathlib import Path
import numpy as np


class Cleaner(Input):
    def __init__(self, args):
        super().__init__(args)

    def create_validation_group(self):
        """
        This will create a group which will contain two sets of data. The first one contains all the valid chromosomes
        for this project via accessing the numbers in the files, the second is all the valid snps found across the
        chromosomes for validation against summary statistics
        """
        # Check parameters and add meaningful error out if process attempts to repeat itself in append mode
        self._common_assert()
        assert self.load_directory, ec.missing_arg(self.operation, "Load_Directory")
        assert self.h5_validation not in self.project_file.keys(), ec.appending_error(self.project_name,
                                                                                      self.h5_validation)

        # Create the validation group
        validation_group = self.project_file.create_group(self.h5_validation)

        # Create a dataset of all the chromosomes we have to work with
        validation_group.create_dataset(self.h5_valid_chromosome, data=self._validation_chromosomes())

        # Create a dataset of all the validation snps
        validation_group.create_dataset(self.h5_valid_snps, data=self._validation_snps())
        self.project_file.close()
        print(f"Create Validation snps and chromosomes: {mc.terminal_time()}")

    def clean_summary_statistics(self):
        # Check parameters, validate that validation has been run, and that clean summary has not.
        assert self.h5_validation in self.project_file.keys(), ec.process_not_run(self.operation, self.project_name,
                                                                                  "create_validation_group")
        assert self.h5_summary not in self.project_file.keys(), ec.appending_error(self.project_name, self.h5_summary)
        assert self.load_file, ec.missing_arg(self.operation, "Load_File")
        self._common_assert()

        return

    def _validation_snps(self):
        """
        This will create a dataset for all the valid snps we need for validating summary statistics as an example.

        :return: A list from a set of all the snps that where found across the chromosomes for a given load_type
        :rtype: List
        """
        valid_snps = []
        for file in mc.directory_iterator(self.load_directory):
            if Path(self.load_directory, file).suffix == self.load_type:
                string_path = str(Path(self.load_directory, file).absolute())
                if self.load_type == ".bed":
                    valid_snps.append(Bed(string_path, count_A1=True).sid)
                elif self.load_type == ".bgen":
                    valid_snps.append([snp.split(",")[0] for snp in Bgen(string_path).sid])

        # Remove duplicates via set
        return np.array(list(set(mc.flatten(valid_snps))), dtype=self.h5_string_type)

    def _validation_chromosomes(self):
        """
        This will create a dataset of all the chromosomes that we have to work with with our validation group in the
        h5py file

        Note
        ------
        This pretty hard coded in how we access the number, should make it clear if people have split via other means

        :return: A list of valid chromosomes
        :rtype: list
        """
        valid_chromosomes = []
        for file in mc.directory_iterator(self.load_directory):
            if Path(self.load_directory, file).suffix == self.load_type:
                valid_chromosomes.append(int(Path(self.load_directory, file).stem.split("_")[-1]))
        valid_chromosomes.sort()
        return valid_chromosomes

    def _common_assert(self):
        """
        Multiple Cleaner methods need to have these variables, so we assert them via a common method call.
        """
        assert self.project_file, ec.missing_arg(self.operation, "Project_Name")
        assert self.load_type, ec.missing_arg(self.operation, "Load_Type")
