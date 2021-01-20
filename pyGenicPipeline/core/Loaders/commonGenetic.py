from pyGenicPipeline.utils import errors as ec
from pyGenicPipeline.utils import misc as mc
from .argsParser import ArgsParser

from bgen_reader import custom_meta_path
from csvObject import CsvObject
from pathlib import Path


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
        self.gen_directory = mc.validate_path(self.args["Load_Directory"])
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
            package_root = Path(__file__).parent.parent

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

    def _set_reference_panel(self):
        if self.args["Reference_Panel"]:
            path_to_file = Path(self.args["Reference_Panel"])
            mc.validate_path(path_to_file, False)
            return CsvObject(path_to_file, set_columns=True).row_data

        else:
            return None
