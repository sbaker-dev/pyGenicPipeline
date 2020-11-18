from pyGeneticPipe.core.Input import Input
from pathlib import Path


class ShellMaker(Input):
    def __init__(self, args):
        super().__init__(args)
        self.file = self._create_shell_file()
        self._create_batch_header()

    def split_bed_by_chromosome(self):
        """
        Create a script to split the current master .bed into separate chromosomes
        """
        self.file.write("# Source: zx8754 - https://www.biostars.org/p/387132/ \n\n")
        self._construct_args([self.plink_key, self.load_key])
        self.file.write("# Load shell modules required to run\n")
        self.file.write(f"module load ${self.plink_key}\n\n")
        self.file.write('# For each chromosome within the file specified create a new file\n')

        # Iterate through chromosomes and generate a .bed file for each one
        self.file.write("for chr in {1..23}; do \\\n")
        self.file.write(f"plink2 --bfile ${self.load_key} --chr $chr --make-bed "
                        f"--out ${{{self.load_key}}}_${{chr}}; \\\n")
        self.file.write(f"done\n")
        self._close()

    def convert_to_bgen(self):
        """
        Creates a script to use QCtoolv2 to convert .bed files to .bgen v1.2 files
        """
        self.file.write("# See other conversions from "
                        "https://www.well.ox.ac.uk/~gav/qctool/documentation/examples/converting.html \n\n")
        self._construct_args([self.qc_key, self.load_key])
        self.file.write("# Load shell modules required to run\n")
        self.file.write(f"module load ${self.qc_key}\n\n")
        self.file.write('# For each .bed chromosome file within this directory convert it to .bgen v1.2\n')

        # Iterate through chromosomes and create a .bgen file for each one
        self.file.write("for chr in {1..23}; do \\\n")
        self.file.write(f"qctool -g ${{{self.load_key}}}_${{chr}}.bed -og ${{{self.load_key}}}_${{chr}}.bgen; \\\n")
        self.file.write(f"done\n")
        self._close()

    def _close(self):
        """
        Close the file and log the process is finished
        """
        self.file.close()
        print(f"Created {self.operation}.sh script")

    def _create_shell_file(self):
        """
        Create an sh file
        """
        file = open(Path(self.working_dir, f"{self.operation}.sh"), "w")
        file.write("#!/bin/bash\n\n")
        file.write("# Generate by pyGeneticPipe/support/ShellMaker.py\n\n")
        return file

    def _construct_args(self, args_list):
        """
        Add the args that have been used for this file to make debugging easier, allow quick switch and editing for
        those who know what went wrong.
        """
        self.file.write("# This pyGeneticPipe job takes the following direct args:\n")
        for arg in args_list:
            self.file.write(f"{arg}={self.args[arg]}\n")
        self.file.write("\n")

    def _create_batch_header(self):
        """
        Different submission types require a different header for batch control, set the header according to the type
        specified in Job_Scheduler
        """
        if self.args["Job_Scheduler"] == "slurm":
            self._set_slurm_header()

        elif self.args["Job_Scheduler"] == "qsub":
            raise NotImplementedError("Qsub not yet supported")

        elif self.args["Job_Scheduler"] == "local":
            print("No batch job requested, skipping header")

        else:
            raise Exception(f"Job_Scheduler takes the arguments slurm, qsub, and local yet was passed "
                            f"{self.args['Job_Scheduler']}")

    def _set_slurm_header(self):
        """
        Set the header for slurm for the current shell script
        """
        self.file.write(f"# Setup batch control for slurm\n")
        self.file.write(f"#SBATCH --job-name={self.operation}\n")
        self.file.write(f"#SBATCH --partition={self.args['Partition']}\n")
        self.file.write(f"#SBATCH --nodes={self.args['Nodes']}\n")

        # todo Expose these
        self.file.write(f"#SBATCH --ntasks-per-node=1\n")
        self.file.write(f"#SBATCH --cpus-per-task=1\n")

        self.file.write(f"#SBATCH --time={self.args['Job_Time']}\n")
        self.file.write(f"#SBATCH --mem={self.args['Job_Memory']}\n\n")

    @property
    def qc_key(self):
        """
        Key for qc path from args
        """
        return "QCTool_Path"

    @property
    def load_key(self):
        """
        Key for Load_File from args
        """
        return "Load_File"

    @property
    def plink_key(self):
        """
        Key for plink_key from args
        """
        return "Plink_Path"
