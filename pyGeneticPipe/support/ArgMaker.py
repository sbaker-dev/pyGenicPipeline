from pyGeneticPipe.utils import error_codes as ec
from pyGeneticPipe.utils.misc import load_yaml
from pathlib import Path
import textwrap


class ArgMaker:
    def __init__(self):
        self._yaml_parameters = load_yaml(Path(Path(__file__).parent, "args.yaml"))
        self.line_width = 120

    def write_yaml_args(self, args_dict, write_directory, validation=True):
        """
        This will write a .yaml file contain all the parameters, split by mandatory and not required, for the current
        operation.
        """
        # If validation, because the args have been set within python via dict replacement, then check all mandatory
        # args have been set
        if validation:
            for key, value in zip(args_dict["Mandatory"].keys(), args_dict["Mandatory"].values()):
                assert value, ec.missing_arg(args_dict["Mandatory"]["Operation"], key)

        # Create the .yaml stub
        file = self._create_yaml_file(args_dict["Mandatory"]["Operation"], write_directory)

        # Write Mandatory args
        self._write_header(file, "MANDATORY ARGS - IF ANY ARE NONE JOB WILL NOT RUN")
        self._write_args(file, args_dict, "Mandatory")

        # Write Optional args
        self._write_header(file, "OPTIONAL ARGS - NOT STRICTLY REQUIRED TO RUN BUT ALTERS OPERATION")
        self._write_args(file, args_dict, "Optional")

        # Write require to run the system, but needed for this operation
        self._write_header(file, "SYSTEM ARGS - WILL NOT BE USED FOR THIS OPERATION BUT REQUIRED TO EXIST IN FILE")
        loaded_keys = list(args_dict["Mandatory"].keys()) + list(args_dict["Optional"].keys())
        args_dict["System"] = {key: None for key in self.all_args if key not in loaded_keys}
        self._write_args(file, args_dict, "System")
        file.close()
        print(f"Constructed file for {args_dict['Mandatory']['Operation']}")

    def write_yaml_config_dict(self, config_dict, write_directory, operation):
        """Write a dict of one level of separation"""
        file = self._create_yaml_file(operation, write_directory)
        for key, values in zip(config_dict.keys(), config_dict.values()):
            if "Description" in values.keys():
                self._write_header(file, values["Description"])

            file.write(f"{key}: \n")
            self._write_args(file, config_dict, key, spacing=1, write_descriptions=False)

    def _make_working_dict(self, key):
        """This will construct the working dict of args that the user needs to submit for a given operation"""
        working_dict = {"Mandatory": self._get_operation_dict(self._yaml_parameters[f"{key}_M"]),
                        "Optional": self._get_operation_dict(self._yaml_parameters[f"{key}_O"])}
        working_dict["Mandatory"]["Operation"] = key
        return working_dict

    def _get_operation_dict(self, operation_keys):
        """Isolate the keys from all operational args that we need for this operation"""
        return {key: value for key, value in zip(self.all_args.keys(), self.all_args.values()) if key in operation_keys}

    def _write_header(self, file, message):
        """Write a header with a message and trailing # up to the line length"""
        file.write(f"# {message} " + f"".join(["#" for _ in range(self.line_width - (len(message) + 3))]) + "\n\n")

    def _write_args(self, file, args_dict, args_type, spacing=0, write_descriptions=True):
        """For a given sub type of args, write the value if set of null otherwise"""
        for key, value in zip(args_dict[args_type].keys(), args_dict[args_type].values()):
            if write_descriptions:
                self._write_description(file, self._arg_descriptions[key])

            if value:
                file.write(f"{''.join([' ' for _ in range(spacing)])}" + f"{key}: {value}\n\n")
            else:
                file.write(f"{''.join([' ' for _ in range(spacing)])}" + f"{key}: null\n\n")

        file.write("\n\n")

    def _write_description(self, file, description_text):
        """Some descriptions need to be wrapped"""
        text_lines = textwrap.wrap(description_text, self.line_width - 2)
        for line in text_lines:
            file.write(f"# {line}\n")

    @staticmethod
    def _create_yaml_file(operation, write_directory):
        """
        Create a number yaml file with a magic number for the current operation
        """
        path_to_file = Path(write_directory)
        assert path_to_file.exists(), ec
        file = open(Path(write_directory, f"{operation}.yaml"), "w")
        file.write("# Generated by pyGeneticPipe/support/ArgMaker.py\n\n")
        return file

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
        return self._make_working_dict("split_bed_by_chromosome")

    @property
    def convert_to_bgen(self):
        """Returns the arguments the user needs to set for convert_to_bgen"""
        return self._make_working_dict("convert_to_bgen")

    @property
    def _arg_descriptions(self):
        """So we can replicate comments in write file"""
        return self._yaml_parameters["Arg_Descriptions"]
