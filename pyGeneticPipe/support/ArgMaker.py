from pyGeneticPipe.utils import error_codes as ec
from pyGeneticPipe.utils.misc import load_yaml
from pathlib import Path


class ArgMaker:
    def __init__(self):
        self._yaml_parameters = load_yaml(Path(Path(__file__).parent, "args.yaml"))

    def write_args(self, args_dict, write_directory, validation=True):
        """
        This will write a .yaml file contain all the parameters, split by mandatory and not required, for the current
        operation.
        """
        # If validation, because the args have been set within python via dict replacement, then check all mandatory
        # args have been set
        if validation:
            for key, value in zip(args_dict.keys(), args_dict.values()):
                assert value, ec.missing_arg(args_dict["Operation"], key)

        # Create the .yaml stub
        file = self._create_yaml_file(args_dict["Operation"], write_directory)

        # Mandatory args:
        file.write("# MANDATORY ARGS - IF ANY ARE NONE JOB WILL NOT RUN ###############################################"
                   "#####################\n\n")
        for key, value in zip(args_dict.keys(), args_dict.values()):
            file.write(f"# {self._arg_descriptions[key]}\n")
            if value:
                file.write(f"{key}: {value}\n\n")
            else:
                file.write(f"{key}: null\n\n")

        # todo need to allow for option args

        # Not Required for this job but will need to set via Input
        file.write("# NOT REQUIRED ARGUMENTS ##########################################################################"
                   "#####################\n")
        file.write("# These arguments are not required for this job, but are required to exist as None within this "
                   "file\n\n")

        for key in self.all_args:
            if key not in args_dict.keys():
                file.write(f"# {self._arg_descriptions[key]}\n")
                file.write(f"{key}: null\n\n")

        file.close()
        print(f"Constructed file for {args_dict['Operation']}")

    @property
    def operations(self):
        """Returns the current defined operations for processing"""
        return self._yaml_parameters["Operations"]

    @property
    def all_args(self):
        """Return all the arguments that are within the yaml file"""
        return self._yaml_parameters["all_args"]

    @property
    def split_bed_by_chromosome(self):
        """Returns the arguments the user needs to set for split_bed_by_chromosome"""
        return self._yaml_parameters["split_bed_by_chromosome"]

    @property
    def convert_to_bgen(self):
        """Returns the arguments the user needs to set for convert_to_bgen"""
        return self._yaml_parameters["convert_to_bgen"]

    @property
    def clean_summary_statistics(self):
        """Returns the arguments the user needs to set for clean_summary_statistics """
        return self._yaml_parameters["clean_summary_statistics"]

    @property
    def _arg_descriptions(self):
        """So we can replicate comments in write file"""
        return self._yaml_parameters["Arg_Descriptions"]

    @staticmethod
    def _create_yaml_file(operation, write_directory):
        """
        Create a number yaml file with a magic number for the current operation
        """
        path_to_file = Path(write_directory)
        assert path_to_file.exists(), ec
        file = open(Path(write_directory, f"{operation}.yaml"), "w")
        file.write("# 48656c6c73696e67\n")
        file.write("# Generated by pyGeneticPipe/support/ArgMaker.py\n\n")
        return file
