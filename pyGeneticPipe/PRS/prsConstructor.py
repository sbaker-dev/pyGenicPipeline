from pyGeneticPipe.PRS.Summary import Summary
from pyGeneticPipe.utils.Input import Input
from pathlib import Path
import h5py


class PRSConstructor(Summary, Input):
    def __init__(self, json_args):
        super().__init__(json_args)

        self.hdf5_file = h5py.File(Path(self.working_dir, self.project_name), "w")

    def prs_from_summary(self):
        self.clean_summary_data()
