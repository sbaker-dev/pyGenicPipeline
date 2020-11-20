from pyGeneticPipe.utils import error_codes as ec
from pyGeneticPipe.utils import misc as mc
from pyGeneticPipe.core.Input import Input
from pysnptools.distreader import Bgen
from pysnptools.snpreader import Bed
from pathlib import Path
import numpy as np
import pickle
import gzip
import re


class Cleaner(Input):
    def __init__(self, args):
        super().__init__(args)

        self._error_dict = {"Invalid_Snps": [], "Chromosome": {}, "Position": {}, "Effect_Size": {}, "P_Value": {},
                            "Standard_Errors": {}, "Duplicate_Position": {}}

    def clean_summary_statistics(self):

        # Check for input arguments
        self._assert_clean_summary_statistics()

        with mc.open_setter(self.summary_file)(self.summary_file) as file:
            # Skip header row
            a = file.readline()
            print(a)
            file.close()

        return

    def _validation_snps(self):
        """
        This will create a dataset for all the valid snps we need for validating summary statistics as an example.

        :return: A list from a set of all the snps that where found across the chromosomes for a given load_type
        :rtype: List
        """

        hap_map_3 = self._load_hap_map_3()

        valid_snps = []
        for file in mc.directory_iterator(self.load_directory):
            if Path(self.load_directory, file).suffix == self.load_type:

                # pysnptools doesn't like path so construct a string of it
                string_path = str(Path(self.load_directory, file).absolute())

                # Load the snps based on load type
                if self.load_type == ".bed":
                    snps = Bed(string_path, count_A1=True).sid
                elif self.load_type == ".bgen":
                    snps = [snp.split(",")[0] for snp in Bgen(string_path).sid]
                else:
                    raise Exception("Unknown load type set")

                # If we only want the hap_map_3 snps then check each snp against the set of hap_map_3
                if hap_map_3:
                    valid_snps.append([s for s in snps if s in hap_map_3])
                else:
                    valid_snps.append(snps)

        # Remove duplicates via set
        return set(mc.flatten(valid_snps))

    def _validation_chromosomes(self):
        """
        This will create a dataset of all the chromosomes that we have to work with our validation group in the
        h5py file

        :return: A list of valid chromosomes
        :rtype: list
        """

        valid_chromosomes = []
        for file in mc.directory_iterator(self.load_directory):
            if Path(self.load_directory, file).suffix == self.load_type:
                valid_chromosomes.append(int(re.sub(r'[\D]', "", Path(self.load_directory, file).stem)))
        valid_chromosomes.sort()
        return valid_chromosomes

    def _load_hap_map_3(self):
        """
        Users may wish to limit valid snps to those found within HapMap3. If they do, they need to provide a path to the
        hapmap3 snp file which will be check that it exists, have the snps extracted and return. Otherwise set to none

        :return: The set of the valid HapMap3 snps or None
        :rtype: set | None
        """

        if self.hap_map_3_file:
            # Construct path as an object and check it exists

            # If the HapMap3 file exists, then extract the snp ids and return them
            f = gzip.open(self.hap_map_3_file, 'r')
            hm3_sids = pickle.load(f)
            f.close()
            return hm3_sids
        else:
            return None

    def _assert_clean_summary_statistics(self):
        """
        clean_summary_statistics requires

        The project file, for writing too
        That the cleaning has not already been undertaken
        That the summary file path exists
        That the load type for the genetic data exists
        That the load directory containing the chromosome split data exists
        """
        # Check parameters, validate that validation has been run, and that clean summary has not.
        assert self.project_file, ec.missing_arg(self.operation, "Project_Name")
        assert self.h5_summary not in self.project_file.keys(), ec.appending_error(self.project_name, self.h5_summary)
        assert self.summary_file, ec.missing_arg(self.operation, "Summary_Path")
        assert self.load_type, ec.missing_arg(self.operation, "Load_Type")
        assert self.load_directory, ec.missing_arg(self.operation, "Load_Directory")
