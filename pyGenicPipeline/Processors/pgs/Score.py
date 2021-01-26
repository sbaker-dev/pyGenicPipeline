from pyGenicPipeline.utils import misc as mc
from pyGenicPipeline.core.Input import Input

from miscSupports import directory_iterator, terminal_time
from csvObject import CsvObject, write_csv
from pyGenicParser import *
from pathlib import Path
import numpy as np


class Score(Input):
    def __init__(self, args):
        super().__init__(args)
        self._score_error_dict = {"Missing IID": 0, "Missing Phenotype": 0}
        self._validation = {}

    def pgs_scores(self):
        """
        This will construct the pgs for a given type, such as infinitesimal, within this chromosome
        """

        # Load the reference to the full sample of ID's, and use it to extract genetic phenotype information
        core = self.gen_reference(self.select_file_on_chromosome())
        ph_dict = self.genetic_phenotypes(core)

        # Load the betas based on the weighted beta type specified by the user
        score_file = f"{self.score_type.split('_')[0]}_{self.target_chromosome}.csv"
        weights = CsvObject(Path(self.working_dir, "PGS", self.score_type, score_file), [str, float],
                            set_columns=True)

        # Chunk the data into memory chunks to be processed
        chunked_snps, chunks = self.chunked_snp_names(weights[self.snp_id], True)
        chunked_weights = np.array_split(weights[self.inf_beta], chunks)

        # Weight the dosage data to construct the scores
        scores = self._weight_dosage(chunked_snps, chunked_weights, core, ph_dict)

        # Combine the FID/IID, genetic phenotype information, and the score for this chromosome
        scores.shape = (len(scores), 1)
        iid_fid = np.array([[v[i] for v in ph_dict.values()] for i in range(core.iid_count)])
        write_out = np.hstack([iid_fid, scores]).tolist()

        # Write this information to a csv
        headers = list(ph_dict.keys()) + ["Scores"]
        write_csv(self.scores_directory, f"Scores_{self.target_chromosome}", headers, write_out)
        print(f"Finished Constructing scores for Chromosome {self.target_chromosome} {terminal_time()}")

    def genetic_phenotypes(self, gen_file):
        """
        Load the full genetic data for this chromosome and isolate any information that can be isolated from it. In this
        case, .bed load types can access more than bgen due to the ability to extract sex from the .fam file.
        """

        ph_dict = {}
        # For plink files, load the fam file then extract the fid, iid and sex information
        if self.gen_type == ".bed":
            ph_dict[self.fam] = np.array(PlinkObject(self.select_file_on_chromosome()).get_family_identifiers())
            ph_dict[self.fid] = mc.variant_array(self.fid.lower(), ph_dict[self.fam])
            ph_dict[self.iid] = mc.variant_array(self.iid.lower(), ph_dict[self.fam])
            ph_dict.pop(self.fam, None)

        # Bgen doesn't have a fam equivalent, so just load the fid and iid
        # todo update to allow for sex and missing if we have loaded .sample
        elif self.gen_type == ".bgen":
            if self._snp_tools:
                ids = gen_file.iid
            else:
                ids = gen_file.iid_array()

            ph_dict[self.fid] = np.array([fid for fid, iid in ids])
            ph_dict[self.iid] = np.array([iid for fid, iid in ids])

        else:
            raise Exception("Unknown load type set")

        return ph_dict

    def _weight_dosage(self, chunked_snps, chunked_weights, core, ph_dict):
        """
        Use the weighted beta values and the dosages from our genetic file to construct a score

        The weights calculated and re-structure the 1 dimensional list to be a vector array. We then use this vector
        array of beta values from weights alongside the raw snps to calculate the effect of each snp based on the
        nucleotide of the individuals (0, 1 or 2) to compute a score for this key.
        """
        scores = np.zeros(len(ph_dict[self.iid]))
        for index, (snp_names, score_weights) in enumerate(zip(chunked_snps, chunked_weights)):
            print(f"Processing Scores Chunk {index} out of {len(chunked_snps)} - {terminal_time()}")

            # Load the raw snps that have been isolated
            raw_snps = self.isolate_raw_snps(core, snp_names)

            # Restructure the scores beta weights to be a vector array
            weights = np.array(score_weights)
            weights.shape = (len(weights), 1)

            # Cumulative the PRS for the individuals
            scores += np.array([np.sum(row) for row in ((-1 * raw_snps) * weights).T])
        return scores

    def aggregate_scores(self):
        """
        This will check the headers where we have full information across all our chromosomes and return the names of
        those headers as a list of strings. Failures will be printed to the terminal.

        The successful headers that where isolated will then be aggregated
        """

        combined_array = []
        for index, file in enumerate(directory_iterator(self.scores_directory)):
            score_file = CsvObject(Path(self.scores_directory, file), set_columns=True)

            # If its the first file we want to extract the iid and fid values as well
            if index == 0:
                fid, iid, score = score_file[self.fid], score_file[self.iid], score_file["Scores"]
                fid, iid, score = np.array(fid), np.array(iid), np.array(score).astype(float)
                fid.shape, iid.shape, score.shape = (len(fid), 1), (len(iid), 1), (len(score), 1)
                combined_array = [fid, iid, score]

            # Else extract the scores and append it to the array
            else:
                score = np.array(score_file["Scores"]).astype(float)
                score.shape = (len(score), 1)
                combined_array.append(score)

        # Combine the (IID_Count, 1) * (chromosome count + 2) arrays into a single (IDD_Count, chromosome_count)
        iid_array = np.hstack(combined_array[:2])
        score_array = np.sum(np.hstack(combined_array[2:]), axis=1)
        score_array.shape = (len(score_array), 1)

        # Write the scores to the working directory
        write_rows = np.hstack([iid_array, score_array]).tolist()
        write_csv(Path(self.working_dir, "PGS"), "PolyGenicScores", ["FID", "IID", "Scores"], write_rows)


