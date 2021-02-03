from pyGenicPipeline.utils import errors as ec
from pyGenicPipeline.utils.misc import variant_array, validate_path
from .argsParser import ArgsParser

from miscSupports import directory_iterator
from bgen_reader import custom_meta_path
from pysnptools.distreader import Bgen
from pysnptools.snpreader import Bed
from csvObject import CsvObject
from collections import Counter
from pyGenicParser import *
from pathlib import Path
import numpy as np
import pickle
import gzip
import re


class CommonGenetic(ArgsParser):
    def __init__(self, args):
        """
        These attributes are or can be used by multiple processors so are not restricted to sub-loaders
        :param args:
        """
        super().__init__(args)

        self.target_chromosome = self.args["Target_Chromosome"]

        # Internal data brought across with package, loaded if requested
        self.hap_map_3 = self._load_local_data_path("HapMap3")
        self.lr_ld_path = self._load_local_data_path("Filter_Long_Range_LD")

        # Genetic files attributes
        self.gen_directory = validate_path(self.args["Load_Directory"])
        self.gen_type = self.args["Load_Type"]
        self._snp_tools = self.args["PySnpTools_Bgen"]

        # PySnpTools will by default write memory files to the directory of the Gen files, but this is often undesirable
        # on a sever environment where read and write permissions may not be universal.
        if self._snp_tools:
            self.make_sub_directory(None, "PySnpTools_Meta")
            custom_meta_path(Path(self.working_dir, "PySnpTools_Meta"))

        # A list of iid to be used for a reference panel if required
        self.ref_panel = self._set_reference_panel()

    def _load_local_data_path(self, access_key):
        """
        This will set a dataset path that has been embedded into the package as a non yaml file sourced from LDPred if
        the current arg is set to be true

        :param access_key: The key to access the data file
        :type access_key: str

        :return: Path to the relevant file if request is not equal to None, else None
        :rtype: Path | None
        """

        if self.args[access_key]:
            package_root = Path(__file__).parent.parent.parent

            if access_key == "HapMap3":
                access_path = Path(package_root, "Data", "hm3_sids.txt.gz")
                assert access_path, ec.path_invalid(access_path, "_load_local_data")
                return access_path
            elif access_key == "Filter_Long_Range_LD":
                access_path = Path(package_root, "Data", "long-range-ld-price-2008hg38.txt")
                assert access_path, ec.path_invalid(access_path, "_load_local_data")
                return access_path
            else:
                raise Exception(f"Unknown Key provided to _load_local_data: {access_key}")
        else:
            return None

    def load_hap_map_3(self):
        """
        Users may wish to limit valid snps to those found within HapMap3. If they do, we access them via the local file
        and return them as a set
        """

        # If the HapMap3 file exists, then extract the snp ids as a set and return them
        f = gzip.open(self.hap_map_3, 'r')
        hm3_sids = pickle.load(f)
        f.close()
        return set(hm3_sids)

    def load_lr_ld_dict(self):
        """
        This will read in the long rang ld dict from Price et al. AJHG 2008 long range LD tables taken from LDPred and
        then filter out the keys relevant to this chromosome.
        """
        long_dict = {chromosome_key: {} for chromosome_key in range(1, 24)}
        with open(str(self.lr_ld_path.absolute()), 'r') as f:
            for line in f:
                chromosome_line, start_pos, end_pos, hild = line.split()
                try:
                    long_dict[int(chromosome_line)][hild] = {'start_pos': int(start_pos), 'end_pos': int(end_pos)}
                except ValueError:
                    continue
        return long_dict

    def _set_reference_panel(self):
        """
        Many operations will need a reference panel of individuals that are genetically dis-similar/ Not related to each
        other. This will load a csv or similar text file with two columns of type FID - IID if set. Else will return
        None.

        Note
        -----
        This operation does NOT allow for headers, so do not set them!
        """
        if self.args["Reference_Panel"]:
            path_to_file = Path(self.args["Reference_Panel"])
            validate_path(path_to_file, False)
            return CsvObject(path_to_file, set_columns=True, file_headers=False).row_data

        else:
            return None

    def select_file_on_chromosome(self):
        """
        For a given chromosome, get the respective file from the genetic directory

        :return: Path to the current file as a Path from pathlib
        """
        for file in directory_iterator(self.gen_directory):
            if Path(self.gen_directory, file).suffix == self.gen_type:
                try:
                    if int(re.sub(r'[\D]', "", Path(self.gen_directory, file).stem)) == self.target_chromosome:
                        return str(Path(self.gen_directory, file).absolute())
                except (ValueError, TypeError):
                    continue

        raise Exception(f"Failed to find any relevant file for {self.target_chromosome} in {self.gen_directory}")

    def gen_reference(self, load_path):
        """Get the pysnptools reference via the load type"""
        if self.gen_type == ".bed":
            return Bed(load_path, count_A1=True)
        elif self.gen_type == ".bgen":
            if self._snp_tools:
                return Bgen(load_path)
            else:
                return BgenObject(load_path)
        else:
            raise Exception("Unknown load type set")

    def construct_reference_panel(self):
        """
        Take in a list of names if provide and index our gen file accordingly, else returns the first 10% of iid samples
        """
        gen_file = self.gen_reference(self.select_file_on_chromosome())
        if self.ref_panel:
            return gen_file[gen_file.iid_to_index(self.ref_panel), :]
        else:
            return gen_file[:int(gen_file.iid_count * 0.1), :]

    def load_variants(self):
        """
        Load variants, for .bgen or plink files, as a set of snps that exist within the current chromosome. Uses the
        validation percentage to construct a validation group, and returns the set of snps for each group. If hap_map_3
        is enabled, it will strip out snps not in hap_map_3.

        We will also need a way to index out the variant information, so we set the indexer according to the load type

        :return: Set of the validation and core set of snps, as well as an indexer to extract information from them
        """

        #  Set validation and core sets of sids based on the load type
        if self.gen_type == ".bed":
            ref = self.construct_reference_panel()
            indexer = PlinkObject(self.select_file_on_chromosome(), True)

        elif self.gen_type == ".bgen":
            indexer = BgenObject(self.select_file_on_chromosome())
            if self._snp_tools:
                print("Loading bgen with PySnpTools\n")
                # Bgen files store [variant id, rs_id], we just want the rs_id hence the [1]; see https://bit.ly/2J0C1kC
                ref = [snp.split(",")[1] for snp in self.construct_reference_panel().sid]

            else:
                print("Loading bgen with custom pybgen via pyGenicParser\n")
                ref = self.construct_reference_panel().sid_array()

        else:
            raise Exception("Unknown load type set")

        # If we only want the hap_map_3 snps then check each snp against the set of hap_map_3
        if self.hap_map_3:
            hap_map_3_snps = self.load_hap_map_3()
            ref = [snp for snp in ref if snp in hap_map_3_snps]

        # Check for duplicates that may be loaded later in the pipeline if we don't filter them out and will otherwise
        # not be detected due to returning a set
        r_count = len(ref)
        ref = [snp for (snp, count) in Counter(ref).items() if count == 1]

        # Count total duplicates
        duplicates = np.sum(r_count - len(ref))
        return set(ref), indexer, duplicates

    def snp_names(self, sm_dict):
        """Variant names differ in pysnptools bgen, so account for this and just return rs_id's"""
        return variant_array(self.snp_id.lower(), sm_dict[self.sm_variants])

    @staticmethod
    def revert_snp_names(snp_names, gen_file):
        """PySnpTools stores snp_names as (snp,snp) and we may need to restore the names for some processes"""
        v_dict = {snp[1]: snp[0] for snp in [snp.split(",") for snp in gen_file.sid]}
        return [f"{v_dict[rs_id]},{rs_id}" for rs_id in snp_names]

    def isolate_raw_snps(self, gen_file, variant_names):
        """
        This will isolate the raw snps for a given bed or bgen file

        :param gen_file: Genetic file you wish to load from
        :param variant_names: The snp names to isolate the dosage for
        :return: raw snps
        """

        # bed returns 2, 1, 0 rather than 0, 1, 2 although it says its 0, 1, 2; so this inverts it
        if self.gen_type == ".bed":
            ordered_common = gen_file[:, gen_file.sid_to_index(variant_names)]
            return np.array([abs(snp - 2) for snp in ordered_common.read(dtype=np.int8).val.T], dtype=np.int8)

        # We have a [1, 0, 0], [0, 1, 0], [0, 0, 1] array return for 0, 1, 2 respectively. So if we multiple the arrays
        # by their index position and then sum them we get [0, 1, 2]
        elif self.gen_type == ".bgen":
            if self._snp_tools:
                # Re-construct the variant_id-rs_id
                variant_names = self.revert_snp_names(variant_names, gen_file)
                ordered_common = gen_file[:, gen_file.sid_to_index(variant_names)]
                return sum(np.array([snp * i for i, snp in enumerate(ordered_common.read(dtype=np.int8).val.T)],
                                    dtype=np.int8))
            else:
                return gen_file.dosage_from_sid(variant_names)

        else:
            raise Exception(f"Critical Error: Unknown load type {self.gen_type} found in _isolate_dosage")

    def normalise_snps(self, gen_file, variant_names, std_return=False):
        """For gibbs we use normalised snps, this process will use the information we have filtered to construct it"""

        raw_snps = self.isolate_raw_snps(gen_file, variant_names)

        # Get the number of snps and individuals in the filtered dict
        n_snps, n_individuals = raw_snps.shape

        # Need to reformat the shape to construct the normalised snps
        raw_means = np.mean(raw_snps, 1, dtype='float32')
        raw_means.shape = (n_snps, 1)
        raw_stds = np.std(raw_snps, 1, dtype='float32')
        raw_stds.shape = (n_snps, 1)

        # Use this information to construct a normalised snps
        normalised_snps = np.array((raw_snps - raw_means) / raw_stds, dtype="float32")
        assert normalised_snps.shape == raw_snps.shape

        if std_return:
            return normalised_snps, raw_stds
        else:
            return normalised_snps, None
