"""
This code is based on LDpred available at https://github.com/bvilhjal/ldpred
"""
from pyGeneticPipe.utils import error_codes as ec
from pyGeneticPipe.utils.Input import Input
from pyGeneticPipe.utils import misc as mc
import math


class Summary(Input):
    def __init__(self, args):
        super().__init__(args)
        self._error_dict = {"Chromosome": {}, "Position": {}, "Effect_Size": {}, "P_Value": {}}

    def clean_summary_data(self):

        print(self.summary_headers)
        # with mc.open_setter(self.summary_path)(self.summary_path) as file:
        #
        #     # Determine if we have custom headers or not via _loaded_sum_headers
        #     raw_headers = file.readline()
        #     headers = {header: self._check_header(header, mc.decode_line(raw_headers, self.zipped))
        #                for header in self._loaded_sum_headers}
        #
        #     for index, line in enumerate(file):
        #         if index % 10000 == 0:
        #             print(f"{index}")
        #
        #         line = mc.decode_line(line, self.zipped)
        #         snp_id = line[headers["SNP_ID"]]
        #
        #         if snp_id in self.valid_snps:
        #             data = self._sum_stats(snp_id, line, headers, self.snp_map)
        #             if data:
        #                 print("Then we add the information")
        #
        #             break

        return 0

    def _sum_stats(self, snp_id, line, headers, snp_pos_map):

        # If we have chromosomes in our summary statistics check the chromosome of the snp against the validation
        chromosome = snp_pos_map[snp_id]['Chromosome']
        if (headers["Chromosome"] is not None) and (line[headers["Chromosome"]] != chromosome):
            self._error_dict["Chromosome"][snp_id] = {"summary_chromosome": line[headers["Chromosome"]],
                                                      "valid_chromosome": snp_pos_map[snp_id]["Chromosome"]}
            return None

        # If we have base pair position in our summary then validate the base par
        position = snp_pos_map[snp_id]['Position']
        if (headers["Position"] is not None) and (line[headers["Position"]] != position):
            self._error_dict["Position"][snp_id] = {"summary_position": line[headers["Position"]],
                                                    "valid_position": snp_pos_map[snp_id]["Position"]}
            return None

        # Set beta unless it is not a finite number
        beta = float(line[headers["Effect_size"]])
        if not math.isfinite(beta):
            self._error_dict["Effect_Size"][snp_id] = {"effect_size": line[headers["Effect_size"]]}
            return None

        # Set p value as long as it is not zero or a non finite number
        p_value = float(line[headers["P_Value"]])
        if not math.isfinite(p_value) or p_value == 0:
            self._error_dict["P_Value"][snp_id] = {"p_value": line[headers["P_Value"]]}
            return None

        if self.frequencies:
            self._sum_stats_frequencies()

        return 0

    def _sum_stats_frequencies(self):
        raise NotImplementedError("Frequencies are not yet implemented")


