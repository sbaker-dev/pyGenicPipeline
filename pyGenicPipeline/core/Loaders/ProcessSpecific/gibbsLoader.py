from pyGenicPipeline.core.Loaders.argsParser import ArgsParser

from pathlib import Path


class GibbLoader(ArgsParser):
    def __init__(self, args):
        super().__init__(args)

        # Set the write information for the weight data
        self.make_sub_directory("PGS", "Weights")
        self.weights_directory = Path(self.working_dir, "PGS", "Weights")

        # Setup the Gibbs processor Config args
        self.gibbs_causal_fractions = self._set_causal_fractions()
        self.gibbs_run = self.args["Gibbs_Run"]
        self.gibbs_iter = self.args["Gibbs_Iter"]
        self.gibbs_burn_in = self.args["Gibbs_Burn_In"]
        self.gibbs_shrink = self.args["Gibbs_Shrink"]
        self.gibbs_zero_jump = self.args["Gibbs_Zero_Jump"]
        self.gibbs_random_seed = self.args["Gibbs_Random_Seed"]
        self.gibbs_tight = self.args["Gibbs_Tight"]
        self.gibbs_breaker = self.args["Gibbs_Breaker"]

        # Set the headers of the write file
        self.gibbs_headers, self._gibbs_header_dict = self._set_gibbs_headers()

    def _set_causal_fractions(self):
        """
        If the user has provided a set of causal fractions of variants to use for the gibbs sampler then use those, else
        use the default that LDPred used.
        """
        if self.args["Gibbs_Causal_Fractions"]:
            return self.args["Gibbs_Causal_Fractions"]
        else:
            return [1, 0.3, 0.1, 0.03, 0.01, 0.003, 0.001]

    def _set_gibbs_headers(self):
        """Construct the headers that will be used in the writing of weights"""

        gibbs_headers = [self.chromosome, self.bp_position, self.snp_id, self.effect_allele, self.alt_allele,
                         self.beta, self.log_odds, self.ld_scores, self.gibbs_beta, self.effect_size]

        gibbs_dict = {header: i for i, header in enumerate(gibbs_headers)}

        return gibbs_headers, gibbs_dict

