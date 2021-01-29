from pyGenicPipeline.utils.misc import set_args, set_current_job
from pyGenicPipeline.utils import errors as ec
from .ArgMaker import ArgMaker

from miscSupports import terminal_time
from pathlib import Path
import os


class Shell:
    def __init__(self, args, write_path):

        # Load the args, defaults for by chromosome operations, and the operation that is being performed.
        self.args = set_args(args)
        self.default_args = self.args["Default_Args"]
        self.operation = set_current_job(self.args["Operation"])

        # Make a sub directory so we don't include the yaml master within the folder needed for upload
        try:
            os.mkdir(Path(write_path, self.operation))
        except FileExistsError:
            pass

        self.local_directory = Path(write_path, self.operation)

    def split_bed_by_chromosome(self):
        """
        Create a script to split the current master .bed into separate chromosomes
        """
        # Create the file
        assert self.operation == "split_bed_by_chromosome", f"Job {self.operation} does not match the shell constructor"
        file = self._create_shell_file("split_bed_by_chromosome")
        self._create_batch_header(file)
        file.write("# Source: zx8754 - https://www.biostars.org/p/387132/ \n\n")

        # Validate the args
        self._construct_args(file, [self._plink_key, self._load_key])
        file.write("# Load shell modules required to run\n")
        file.write(f"module load ${self._plink_key}\n\n")
        file.write('# For each chromosome within the file specified create a new file\n')

        # Iterate through chromosomes and generate a .bed file for each one
        file.write("for chr in {1..23}; do \\\n")
        file.write(f"plink2 --bfile ${self._load_key} --chr $chr --make-bed "
                   f"--out ${{{self._load_key}}}_${{chr}}; \\\n")
        file.write(f"done\n")
        self._close(file)

    def convert_to_bgen(self):
        """
        Creates a script to use QCtoolv2 to convert .bed files to .bgen v1.2 files
        """
        # Create the file
        assert self.operation == "convert_to_bgen", f"Job {self.operation} does not match the shell constructor"
        file = self._create_shell_file("convert_to_bgen")
        self._create_batch_header(file)
        file.write("# See other conversions from "
                   "https://www.well.ox.ac.uk/~gav/qctool/documentation/examples/converting.html \n\n")

        # Validate the args
        self._construct_args(file, [self._qc_key, self._load_key, self._bgenix_key])
        file.write("# Load shell modules required to run\n")
        file.write(f"module load ${self._qc_key}\n\n")
        file.write('# For each .bed chromosome file within this directory convert it to .bgen v1.2\n')

        # Iterate through chromosomes and create a .bgen file for each one
        file.write("for chr in {1..23}; do \\\n")
        file.write(f"qctool -g ${{{self._load_key}}}_${{chr}}.bed -og ${{{self._load_key}}}_${{chr}}.bgen; \\\n")
        file.write(f"${{BGENIX_Path}} -g ${{{self._load_key}}}_${{chr}}.bgen -index; \\\n")
        file.write(f"done\n")
        self._close(file)

    def chromosome_controllers(self):

        # Load the python controller
        py_file = self._create_py_file("py_controller")
        py_file.write("from pyGenicPipeline import Main\n")
        py_file.write("import sys\n")
        py_file.write("Main(sys.argv[1])\n")
        py_file.close()

        # Set the Submission bash file is going to be batch submitted via sbatch for example
        file = self._create_shell_file(f"bash_batch_submission")
        self._create_batch_header(file)

        file.write("echo 'Starting Job'\n\n")
        file.write("# first source bashrc (with conda.sh), then conda can be used\n")
        file.write("source ~/.bashrc\n\n")

        file.write("# make sure conda base is activated\n")
        # Todo we need to exposure the environment name
        file.write("conda activate pyGenicPipe\n\n")

        file.write("# Set the file_path name as the first argument passed from central bash controller\n")
        file.write("FILE_PATH=$1\n\n")

        file.write("# Run the script\n")
        file.write("python py_controller.py $FILE_PATH\n")
        file.close()

        # Set the central base controller, set the operation to the current operation to be run
        master = self._create_shell_file("bash_controller")
        master.write(f"operation={self.operation}\n")

        # Iterate through files in the given job folder, submitting the file_name as the argument to the bash controller
        # If the file_name is the same as the operator
        master.write(f"for file_name in *\n")
        master.write(f"do\n\n")
        master.write("if [[ $file_name == *$operation* ]]\n")
        master.write("then\n")
        master.write('echo "Submitting $file_name"\n')
        master.write('sbatch bash_batch_submission.sh "`pwd`/$file_name"\n')
        master.write("fi\n\n")
        master.write(f"done\n")
        master.close()

        # Set the yaml files to be by chromosome in this local directory
        args_maker = ArgMaker(self.default_args)

        for i in range(1, 23):
            self.args["Target_Chromosome"] = i
            args_maker.write_yaml_args(self.args, self.local_directory, False, f"{self.operation}_{i}")

    def _close(self, file):
        """
        Close the file and log the process is finished
        """
        file.close()
        print(f"Created {self.operation}.sh script {terminal_time()}")

    def _create_shell_file(self, file_name):
        """
        Create an sh file called file_name
        """
        file = open(Path(self.local_directory, f"{file_name}.sh"), "w")
        file.write("#!/bin/bash\n\n")
        file.write("# Generate by pyGeneticPipe/Utils/Shell.py\n\n")
        return file

    def _create_py_file(self, file_name):
        """
        Create an py file called file_name
        """
        file = open(Path(self.local_directory, f"{file_name}.py"), "w")
        file.write("#!/usr/bin/env python\n\n")
        return file

    def _construct_args(self, file, args_list):
        """
        Add the args that have been used for this file to make debugging easier, allow quick switch and editing for
        those who know what went wrong.
        """
        file.write("# This pyGeneticPipe job takes the following direct args:\n")
        for arg in args_list:
            assert self.args[arg], ec.missing_arg(self.operation, arg)
            file.write(f"{arg}={self.args[arg]}\n")
        file.write("\n")

    def _create_batch_header(self, file):
        """
        Different submission types require a different header for batch control, set the header according to the type
        specified in Job_Scheduler
        """
        if self.args["Job_Scheduler"] == "slurm":
            self._set_slurm_header(file)

        elif self.args["Job_Scheduler"] == "qsub":
            raise NotImplementedError("Qsub not yet supported")

        elif self.args["Job_Scheduler"] == "local":
            print("No batch job requested, skipping header")

        else:
            raise Exception(f"Job_Scheduler takes the arguments slurm, qsub, and local yet was passed "
                            f"{self.args['Job_Scheduler']}")

    def _set_slurm_header(self, file):
        """
        Set the header for slurm for the current shell script
        """
        file.write(f"# Setup batch control for slurm\n")
        file.write(f"#SBATCH --job-name={self.operation}\n")
        file.write(f"#SBATCH --partition={self.args['Partition']}\n")
        file.write(f"#SBATCH --nodes={self.args['Nodes']}\n")

        # todo Expose these
        file.write(f"#SBATCH --ntasks-per-node=1\n")
        file.write(f"#SBATCH --cpus-per-task=1\n")

        file.write(f"#SBATCH --time={self.args['Job_Time']}\n")
        file.write(f"#SBATCH --mem={self.args['Job_Memory']}\n\n")

    @property
    def _qc_key(self):
        """
        Key for qc path from args
        """
        return "QCTool_Path"

    @property
    def _load_key(self):
        """
        Key for Load_File from args
        """
        return "Load_File"

    @property
    def _plink_key(self):
        """
        Key for plink_key from args
        """
        return "Plink_Path"

    @property
    def _bgenix_key(self):
        """
        Key for working with bgenix
        """
        return "BGENIX_Path"
