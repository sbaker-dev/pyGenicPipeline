from pyGenicPipeline.utils import errors as ec

from miscSupports import load_yaml
from pathlib import Path
import textwrap


class ArgMaker:
    def __init__(self, defaults=None):
        self._defaults = self._set_defaults(defaults)
        self.yaml_parameters = load_yaml(Path(Path(__file__).parent, "args.yaml"))
        self.line_width = 120

    def __repr__(self):
        """Human readable print out of class object"""
        return f"ArgMaker with {len(self._defaults.keys())} Default args"

    def write_yaml_args(self, args_dict, write_directory, validation=True, set_write_name=False):
        """
        This will write a .yaml file contain all the parameters, split by mandatory and not required, for the current
        operation.
        """
        # Args dict may come from user or from process. If from process it will not have mandatory or optional keys so
        # we need to set them
        args_dict = self._configure_read_dict(args_dict)

        # If validation, because the args have been set within python via dict replacement, then check all mandatory
        # args have been set
        if validation:
            for key, value in zip(args_dict["Mandatory"].keys(), args_dict["Mandatory"].values()):
                assert value is not None, ec.missing_arg(args_dict["Mandatory"]["Operation"], key)

        # Create the .yaml stub called the operation if no name is set, else the set name
        if set_write_name:
            write_name = set_write_name
        else:
            write_name = args_dict["Mandatory"]["Operation"]
        file = self._create_yaml_file(write_name, write_directory)

        # Write Mandatory args
        self._write_header(file, "MANDATORY ARGS - IF ANY ARE NONE JOB WILL NOT RUN")
        self._write_args(file, args_dict["Mandatory"])

        # Write Optional args
        self._write_header(file, "OPTIONAL ARGS - NOT STRICTLY REQUIRED TO RUN BUT ALTERS OPERATION")
        self._write_args(file, args_dict["Optional"])

        # Write require to run the system, but needed for this operation
        self._write_header(file, "SYSTEM ARGS - WILL NOT BE USED FOR THIS OPERATION AND CAN STAY AS NULL BUT REQUIRED "
                                 "TO EXIST IN FILE")
        loaded_keys = list(args_dict["Mandatory"].keys()) + list(args_dict["Optional"].keys())
        args_dict["System"] = {key: None for key in self.all_args if key not in loaded_keys}
        self._write_args(file, args_dict["System"])

        # Finish the file
        file.close()
        print(f"Constructed yaml file for {args_dict['Mandatory']['Operation']}")

    def _configure_read_dict(self, read_dict):
        """
        If the file was parsed in as yaml then it will not have the Mandatory Optional parameters and these will need to
        be re-set for write_yaml_args. If it was the input from the properties of ArgMaker, we can just return the dict
        that will have been passed

        :param read_dict: Dict of Dicts from ArgMaker properties or Dict from Yaml create from ArgMaker
        :type read_dict: dict

        :return: Configured Dict with Mandatory and Optional Tags
        """

        if ("Mandatory" in read_dict.keys()) and ("Optional" in read_dict.keys()):
            return read_dict
        else:
            # Yaml Messes with False values by setting to null which will trip the Validation
            # Extract the formatted file for this operation
            formatted = self._make_working_dict(read_dict["Operation"])

            # Update any keys based on what was submitted
            for variable_type in formatted.keys():
                for parameter in formatted[variable_type].keys():
                    for attribute in read_dict.keys():
                        if attribute == parameter:
                            formatted[variable_type][parameter] = read_dict[attribute]

            return formatted

    def _make_working_dict(self, operation):
        """This will construct the working dict of args that the user needs to submit for a given operation"""

        # isolate the job data from our args dict for this operation
        data = self._get_job_location(operation)

        # Set the mandatory and optional args from the yaml_parameters
        working_dict = {"Mandatory": self._get_operation_dict(data["Mandatory"]),
                        "Optional": self._get_operation_dict(data["Optional"])}

        # Set the operation to be the operation key
        working_dict["Mandatory"]["Operation"] = operation

        # Update any keys based on the defaults
        for key in working_dict.keys():
            for attribute in working_dict[key].keys():
                if attribute in self._defaults.keys():
                    working_dict[key][attribute] = self._defaults[attribute]

        return working_dict

    def _get_operation_dict(self, operation_keys):
        """Isolate the keys from all operational args that we need for this operation"""
        return {key: value for key, value in zip(self.all_args.keys(), self.all_args.values()) if key in operation_keys}

    def _write_args(self, file, args_dict, spacing=0, write_descriptions=True):
        """
        For a given sub type of args, write the value if set of null otherwise

        Args are group into mandatory, optional, and system. This takes the file object, and the dict of keys of a given
        group to iterate though. If write descriptions are set, then for each key in the dict a set of lines explaining
        what the key means. Then, if a value has been set that is not None, the value is writen to the file else null
        which is yaml for None.

        A Space in yaml represents a tab in most other languages, and spacing can be set to values above zero if you
        want to create items within dict for example.
        """
        for key, value in zip(args_dict.keys(), args_dict.values()):
            if write_descriptions:
                self._write_description(file, self._arg_descriptions[key])

            # To ensure Boolean False is not set to null
            if value is not None:
                file.write(f"{''.join([' ' for _ in range(spacing)])}" + f"{key}: {value}\n\n")
            else:
                file.write(f"{''.join([' ' for _ in range(spacing)])}" + f"{key}: null\n\n")

        file.write("\n\n")

    def _write_description(self, file, description_text):
        """Some descriptions need to be wrapped"""
        text_lines = textwrap.wrap(description_text, self.line_width - 2)
        for line in text_lines:
            file.write(f"# {line}\n")

    def _write_header(self, file, message):
        """Write a header with a message and trailing # up to the line length"""
        file.write(f"# {message} " + f"".join(["#" for _ in range(self.line_width - (len(message) + 3))]) + "\n\n")

    def write_yaml_group_dict(self, config_dict, write_directory, file_name):
        """
         Writes a yaml file called file_name from config_dict of one level of separation to write_directory

        Some operations need to have common attributes stored of multiple parameters. This takes a dict where each key
        in the dict is a parameter you wish to set, and each value is a dict of attributes you wish to set to this key.
        """
        file = self._create_yaml_file(file_name, write_directory)

        # For each key, set the keys parameters by the dict assigned as the value
        for key, values in zip(config_dict.keys(), config_dict.values()):
            if "Description" in values.keys():
                self._write_header(file, values["Description"])

            file.write(f"{key}: \n")
            self._write_args(file, config_dict[key], spacing=1, write_descriptions=False)

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

    def _get_job_location(self, job):
        """Isolate the parameters of the job"""
        for job_type in self.yaml_parameters["Operations"]:
            for operation in self.yaml_parameters["Operations"][job_type]:
                if operation == job:
                    return self.yaml_parameters["Operations"][job_type][operation]

    @staticmethod
    def _set_defaults(defaults):
        """Set the default args dict"""
        if isinstance(defaults, dict):
            return defaults
        elif defaults is None:
            return {}
        else:
            return load_yaml(defaults)

    @property
    def _arg_descriptions(self):
        """So we can replicate comments in write file"""
        return self.yaml_parameters["Arg_Descriptions"]

    @property
    def operations(self):
        """Returns the current defined operations for processing"""
        return self.yaml_parameters["Operations"]

    @property
    def all_args(self):
        """Return all the arguments that are within the yaml file"""
        return self.yaml_parameters["all_args"]

    # Shell Commands
    @property
    def split_bed_by_chromosome(self):
        """Returns the arguments the user needs to set for split_bed_by_chromosome"""
        return self._make_working_dict("split_bed_by_chromosome")

    @property
    def convert_to_bgen(self):
        """Returns the arguments the user needs to set for convert_to_bgen"""
        return self._make_working_dict("convert_to_bgen")

    # PGS by chromosome Commands
    @property
    def pgs_clean_summary_stats(self):
        """Returns the arguments the user needs to set for pgs_clean_summary_stats"""
        return self._make_working_dict("pgs_clean_summary_stats")

    @property
    def pgs_filter_snps(self):
        """Returns the arguments the user needs to set for pgs_filter_snps"""
        return self._make_working_dict("pgs_filter_snps")

    @property
    def pgs_calculate_ld(self):
        """Returns the arguments the user needs to set for pgs_calculate_ld_score"""
        return self._make_working_dict("pgs_calculate_ld")

    @property
    def pgs_infinitesimal(self):
        """Returns the arguments the user needs to set for pgs_infinitesimal"""
        return self._make_working_dict("pgs_infinitesimal")

    @property
    def pgs_scores(self):
        """Returns the arguments the user needs to set for pgs_scores"""
        return self._make_working_dict("pgs_scores")

    # Not by Chromosome PGS Commands
    @property
    def suggest_ld_radius(self):
        """Returns the arguments the user needs to set for suggest_ld_radius"""
        return self._make_working_dict("suggest_ld_radius")

    @property
    def calculate_genome_wide_heritability(self):
        """Returns the arguments the user needs to set for calculate_genome_wide_heritability"""
        return self._make_working_dict("calculate_genome_wide_heritability")

    @property
    def aggregate_scores(self):
        """Returns the arguments the user needs to set for aggregate_scores"""
        return self._make_working_dict("aggregate_scores")
