from pyGeneticPipe.utils.misc import terminal_time
from pyGeneticPipe.PRS.CleanSummary import CleanSummary
from pyGeneticPipe.PRS.ProcessSummary import ProcessSummary
from pyGeneticPipe.utils.Input import Input
from pathlib import Path


class PRSConstructor(ProcessSummary, Input):
    def __init__(self, json_args):
        super().__init__(json_args)

        self._clean_summary_name = "Summary_stats"

    def prs_from_summary(self):
        print(f"Starting PRS construction {terminal_time()}")

        # First clean the summary statistics
        if not Path(self.working_dir, f"{self.project_name}_{self._clean_summary_name}").exists():
            print("Cleaning GWAS summary statistics")
            # self.clean_summary_data()
            CleanSummary(self.args).clean_summary_data(self._clean_summary_name)
        else:
            print(f"Found, so will use, pre-existing Summary statistics for project: {self.project_name}")

        # Now we pre-process the summary statistics file
        if self.ld_ref_mode == "plink":
            self.pre_process_plink(Path(self.working_dir, f"{self.project_name}_{self._clean_summary_name}"))
        else:
            self.pre_process_bgen()

