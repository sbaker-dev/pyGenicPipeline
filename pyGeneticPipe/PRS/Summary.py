"""
This code is based on LDpred available at https://github.com/bvilhjal/ldpred
"""
from pyGeneticPipe.utils.Input import Input
from pyGeneticPipe.utils import misc as mc
import math


class Summary(Input):
    def __init__(self, args):
        super().__init__(args)
        self._error_dict = {"Chromosome": {}, "Position": {}, "Effect_Size": {}, "P_Value": {}}

    def clean_summary_data(self):

        with mc.open_setter(self.summary_path)(self.summary_path) as file:
            # Skip header row
            file.readline()

            # For each line in the GWAS Summary file
            for index, line in enumerate(file):
                if index % 10000 == 0:
                    print(f"{index}")

                line = mc.decode_line(line, self.zipped)
                snp_id = line[self.snp_id]
                print(snp_id)

                if snp_id in self.valid_snps:
                    data = self._sum_stats(snp_id, line, self.snp_map)
                    if data:
                        print("Then we add the information")

                break

        return 0

    def _sum_stats(self, snp_id, line, snp_pos_map):

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
        if not math.isfinite(beta):
            self._error_dict["Effect_Size"][snp_id] = {"effect_size": line[self.effect_size]}
            return None

        # Set p value as long as it is not zero or a non finite number
        p_value = float(line[self.p_value])
        if not math.isfinite(p_value) or p_value == 0:
            self._error_dict["P_Value"][snp_id] = {"p_value": line[self.p_value]}
            return None

        if self.frequencies:
            self._sum_stats_frequencies()

        return 0

    def _sum_stats_frequencies(self):
        raise NotImplementedError("Frequencies are not yet implemented")


