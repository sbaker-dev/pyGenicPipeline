from pyGeneticPipe.geneticParsers.plink.plinkObject import PlinkObject
from pyGeneticPipe.geneticParsers.bgen.bgenObject import BgenObject
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

        valid_chromosomes = self._validation_chromosomes()

        for chromosome in valid_chromosomes:
            print(chromosome)

            load_path = str(self._select_file(chromosome))
            validation, core, indexer = self._load_variants(load_path)

            with mc.open_setter(self.summary_file)(self.summary_file) as file:
                # Skip header row
                file.readline()

                # For each line in the GWAS Summary file
                for index, line in enumerate(file):
                    if index % 10000 == 0 and self.debug:
                        print(f"{index}")

            break


    def _select_file(self, chromosome):
        """
        For a given chromosome, get the respective file
        :param chromosome: Current chromosome to be loaded
        :return: Path to the current file as a Path from pathlib
        """
        for file in mc.directory_iterator(self.load_directory):
            if Path(self.load_directory, file).suffix == self.load_type:
                if int(re.sub(r'[\D]', "", Path(self.load_directory, file).stem)) == chromosome:
                    return Path(self.load_directory, file)

        raise Exception(f"Failed to find any relevant file for {chromosome} in {self.load_directory}")

    def _set_chromosome_dict(self):
        valid_chromosomes = np.load(f"{self.working_dir}/{self.h5_valid_chromosome}.npy", allow_pickle=True).tolist()

        chromosome_dict = {str(chromosome): {self.snp_id: [], self.bp_position: [], self.p_value: [], self.log_odds: [],
                                             self.beta: [], self.nucleotide: [], self.info: [], self.frequency: []}
                           for chromosome in valid_chromosomes}

        # Load the valid snps to test the summary stats against
        valid_snps = np.load(f"{self.working_dir}/{self.h5_valid_snps}.npy", allow_pickle=True).tolist()

        return valid_snps, chromosome_dict

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

    def _load_variants(self, load_path):
        """
        Load variants, for .bgen or plink files, as a set of snps that exist within the current chromosome. Uses the
        validation percentage to construct a validation group, and returns the set of snps for each group. If hap_map_3
        is enabled, it will strip out snps not in hap_map_3.

        We will also need a way to index out the variant information, so we set the indexer according to the load type

        :param load_path: Current Chromosome file
        :return: Set of the validation and core set
        """
        hap_map_3 = self._load_hap_map_3()

        #  Set validation and core sets of sids based on the load type
        if self.load_type == ".bed":
            validation_size = self._set_validation_sample_size(Bed(load_path, count_A1=True).iid_count)
            validation = Bed(load_path, count_A1=True)[:validation_size, :].sid
            core = Bed(load_path, count_A1=True)[validation_size:, :].sid
            indexer = None
            raise NotImplementedError("No indexer set for .bed yet")

        elif self.load_type == ".bgen":
            # Bgen files store [variant id, rsid], we just want the rsid hence the [1]; see https://bit.ly/2J0C1kC
            validation_size = self._set_validation_sample_size(Bgen(load_path).iid_count)
            validation = [snp.split(",")[1] for snp in Bgen(load_path)[:validation_size, :].sid]
            core = [snp.split(",")[1] for snp in Bgen(load_path)[validation_size:, :].sid]
            indexer = BgenObject(load_path)

        else:
            raise Exception("Unknown load type set")

        # If we only want the hap_map_3 snps then check each snp against the set of hap_map_3
        if hap_map_3:
            validation = [snp for snp in validation if snp in hap_map_3]
            core = [snp for snp in core if snp in hap_map_3]

        return set(validation), set(core), indexer

    def _load_hap_map_3(self):
        """
        Users may wish to limit valid snps to those found within HapMap3. If they do, they need to provide a path to the
        hapmap3 snp file which will be check that it exists, have the snps extracted and return. Otherwise set to none

        :return: The set of the valid HapMap3 snps or None
        :rtype: set | None
        """

        if self.hap_map_3_file:
            print("WARNING: THIS IS UNTESTED CODE FROM LDPRED")
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
        That the Validation_Size has been set
        """
        # Check parameters, validate that validation has been run, and that clean summary has not.
        assert self.project_file, ec.missing_arg(self.operation, "Project_Name")
        assert self.h5_summary not in self.project_file.keys(), ec.appending_error(self.project_name, self.h5_summary)
        assert self.summary_file, ec.missing_arg(self.operation, "Summary_Path")
        assert self.load_type, ec.missing_arg(self.operation, "Load_Type")
        assert self.load_directory, ec.missing_arg(self.operation, "Load_Directory")
        assert self.validation_size, ec.missing_arg(self.operation, "Validation_Size")

    def _set_validation_sample_size(self, full_sample_size):
        """
        This will return the value of the validation in terms of individuals rather than a percentage that the user
        specified

        :param full_sample_size: An integer of the number of samples in the full sample
        :type full_sample_size: int

        :return: Integer of the number of samples required in the validation sample
        :rtype: int
        """
        return int(full_sample_size * self.validation_size)

    def _validate_summary_line(self):
        pass
