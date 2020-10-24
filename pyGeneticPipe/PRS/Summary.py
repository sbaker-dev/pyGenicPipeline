"""
This code is based on LDpred available at https://github.com/bvilhjal/ldpred
"""
from pyGeneticPipe.utils.Input import Input
from pyGeneticPipe.utils import misc as mc
import numpy as np
from scipy import stats


class Summary(Input):
    def __init__(self, args):
        super().__init__(args)
        self._error_dict = {"Invalid_Snps": [], "Chromosome": {}, "Position": {}, "Effect_Size": {}, "P_Value": {},
                            "Standard_Errors": {}, "Duplicate_Position": {}}

        self._chromosome_dict = {chromosome: {"snp_id": [], "position": [], "p_value": [], "log_odds": [], "beta": [],
                                              "nucleotide": [], "info": [], "frequency": []}
                                 for chromosome in self.valid_chromosomes}

    def clean_summary_data(self):

        with mc.open_setter(self.summary_path)(self.summary_path) as file:
            # Skip header row
            file.readline()

            # For each line in the GWAS Summary file
            for index, line in enumerate(file):
                if index % 10000 == 0:
                    print(f"{index}")
                    if index != 0:
                        break

                # Decode the line and extract the snp_id
                line = mc.decode_line(line, self.zipped)
                snp_id = line[self.snp_id]

                if snp_id in self.valid_snps:
                    self._validate_line(snp_id, line, self.snp_map)
                else:
                    self._error_dict["Invalid_Snps"].append(snp_id)


        return 0

    def _validate_line(self, snp_id, line, snp_pos_map):

        # If we have chromosomes in our summary statistics check the chromosome of the snp against the validation
        chromosome = snp_pos_map[snp_id]['Chromosome']
        if (self.chromosome is not None) and (line[self.chromosome] != chromosome):
            self._error_dict["Chromosome"][snp_id] = {"summary_chromosome": line[self.chromosome],
                                                      "valid_chromosome": snp_pos_map[snp_id]["Chromosome"]}
            return None

        # If we have base pair position in our summary then validate the base par
        position = snp_pos_map[snp_id]['Position']
        if (self.bp_position is not None) and (line[self.bp_position] != position):
            self._error_dict["Position"][snp_id] = {"summary_position": line[self.bp_position],
                                                    "valid_position": snp_pos_map[snp_id]["Position"]}
            return None

        # Set beta unless it is not a finite number
        beta = float(line[self.effect_size])
        if not np.isfinite(beta):
            self._error_dict["Effect_Size"][snp_id] = {"effect_size": line[self.effect_size]}
            return None

        # Set p value as long as it is not zero or a non finite number
        p_value = float(line[self.p_value])
        if not np.isfinite(p_value) or p_value == 0:
            self._error_dict["P_Value"][snp_id] = {"p_value": line[self.p_value]}
            return None

        # calculate the beta and beta odds
        beta, beta_odds = self._calculate_beta(beta, line, snp_id, p_value)
        if not beta or not beta_odds:
            return None

        # Set the nucleotides
        nucleotides = [line[self.effect_allele].upper(), line[self.alt_allele].upper()]

        # All necessary information has been found, but check for info and frequency setting to -1 otherwise
        # Get the INFO score if it exists
        if self.info is not None:
            info = float(line[self.info])
        else:
            info = -1

        if self.frequencies:
            frequency = self._sum_stats_frequencies()
        else:
            frequency = -1

        # remove duplicated positions
        if position in self._chromosome_dict[chromosome]["position"]:
            self._error_dict["Duplicate_Position"][snp_id] = {"position": position}
            return None
        else:
            self._chromosome_dict[chromosome]["snp_id"].append(snp_id)
            self._chromosome_dict[chromosome]["position"].append(position)
            self._chromosome_dict[chromosome]["p_value"].append(p_value)
            self._chromosome_dict[chromosome]["log_odds"].append(beta_odds)
            self._chromosome_dict[chromosome]["beta"].append(beta)
            self._chromosome_dict[chromosome]["nucleotide"].append(nucleotides)
            self._chromosome_dict[chromosome]["info"].append(info)
            self._chromosome_dict[chromosome]["frequency"].append(frequency)

    def _sum_stats_frequencies(self):
        raise NotImplementedError("Frequencies are not yet implemented")

    def _beta_by_type(self, beta):
        """
        If we are working with ods rations we need to take the log of the read beta
        """
        if self.effect_type == "OR":
            return np.log(beta)
        else:
            return beta

    def _beta_from_se(self, beta, beta_typed, se):
        """Calculate z score with standard error"""
        if self.effect_type == "OR":
            abs_beta = np.absolute(1 - beta) / se
        else:
            abs_beta = np.absolute(beta) / se
        return np.sign(beta_typed) * (abs_beta / np.sqrt(self.sample_size))

    def _calculate_beta(self, beta, line, snp_id, p_value):
        """
        Calculate both the beta, and the beta odds depending on the effect_type and if the user wants to constructed a
        standardised z score or not

        :param beta: beta from summary stats
        :type beta: float

        :param line: The line of the current file
        :type line: list

        :param snp_id: the current snp_id for writing errors if z_scores selected and se is invalid
        :type snp_id: str

        :param p_value: p value from summary stats
        :type p_value: float

        :return: If succesfful, return the beta and beta odds otherwise None and None
        """

        beta_odds = self._beta_by_type(beta)

        # If effect type is Best linear unbiased prediction (BLUP) return beta
        if self.effect_type == "BLUP":
            return beta, beta_odds

        # If we want to compute z scores, compute them as long as standard errors are valid
        elif self.z_scores:
            se = float(line[self.standard_errors])
            if not np.isfinite(se) or se == 0:
                self._error_dict["Standard_Errors"][snp_id] = {line[self.standard_errors]}
                return None, None
            else:
                return self._beta_from_se(beta, beta_odds, se), beta_odds

        # Otherwise compute beta from p values
        else:
            return np.sign(beta) * (stats.norm.ppf(p_value / 2.0) / np.sqrt(self.sample_size)), beta_odds
