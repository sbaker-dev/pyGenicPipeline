from pyGeneticPipe.PRS.Summary import Summary
from pyGeneticPipe.utils.Input import Input


class PRSConstructor(Summary, Input):
    def __init__(self, json_args):
        super().__init__(json_args)

    def prs_from_summary(self):
        self.clean_summary_data()
