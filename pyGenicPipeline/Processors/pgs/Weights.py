from pyGenicPipeline.utils import misc as mc
from pyGenicPipeline.core.Input import Input

from miscSupports import load_pickle, terminal_time, flip_list
from csvObject import write_csv
import numpy as np
import time


class Weights(Input):
    def __init__(self, args):
        super().__init__(args)
        self.start_time = 0

    def pgs_infinitesimal(self):
        """
        This will weight the beta values conditional on LD using infinitesimal shrinkage
        """
        # Check to see if the ld radius is sufficient (sufficiency inferred from LD-pred)
        t0 = time.time()
        if self.genomic_data[self.genome][self.avg_ld] > self.ld_radius / 10.0:
            raise ValueError("WARNING LD Radius appears small compare to the genome wide average LD estimate\n"
                             "Increase the LD radius or use less snps")

        # Re load the filtered snps, LD information for this chromosome, and genome-wide data into memory
        sm_dict = self.sm_dict_from_csv(self.filter_directory, f"Filtered_{self.target_chromosome}.csv")
        ld_data = load_pickle(self.ld_directory, f"LD{self.target_chromosome}")
        gen_data = self.genomic_data[self.target_chromosome]

        # Isolate the information we need to construct the inf-decimal weights for this chromosome
        sid_count, iid_count, herit = gen_data[self.sid_count], gen_data[self.iid_count], gen_data[self.herit]
        normalised_snps, snp_stds = ld_data[self.norm_snps], ld_data[self.raw_stds]

        # Update the betas via infinitesimal shrinkage, weight by standard errors
        updated_betas = self._infinitesimal_betas(herit, iid_count, normalised_snps, sid_count, sm_dict)
        infinitesimal = updated_betas / snp_stds.flatten()

        # Write the snp name - constructed betas to a csv
        write_out = flip_list([self.snp_names(sm_dict), infinitesimal])
        write_csv(self.inf_directory, f"Inf_{self.target_chromosome}", [self.snp_id, self.inf_beta], write_out)
        print(f"Constructed infinitesimal weights in {time.time() - t0} at {terminal_time()}")

    def _infinitesimal_betas(self, herit, iid_count, normalised_snps, sid_count, sm_dict):
        """
        This uses a rolling window, twice the length of the ld_radius, to weight the snp betas via infinitesimal
        weighting
        """
        updated_betas = np.empty(sid_count)
        for index, window in enumerate(range(0, sid_count, self.ld_radius * 2)):
            # Set the window parameters
            start_i = window
            stop_i = min(sid_count, window + self.ld_radius * 2)
            current_window = stop_i - start_i

            # Isolate the snps in the current window
            snp_window = self.window_values(normalised_snps, sid_count, window)

            # Calculate the disequilibrium, then append the shrunk r2 matrix to the ld_matrix
            disequilibrium_dp = np.dot(snp_window, snp_window.T) / iid_count
            dp_shrink = mc.shrink_r2_matrix(disequilibrium_dp, iid_count)

            # numpy.eye is just an identity matrix, don't know what 'a' standards for in LDPred
            a = np.array(((sid_count / herit) * np.eye(current_window)) + (self.sample_size / 1.0) * dp_shrink)

            # Update the betas
            beta_updated = np.dot(np.linalg.pinv(a) * self.sample_size, sm_dict[self.beta][start_i: stop_i])
            updated_betas[start_i: stop_i] = beta_updated
            print(f"Window {window} / {sid_count} Constructed at {terminal_time()}")

        return updated_betas
