from pyGenicPipeline.utils import error_codes as ec
from pyGenicPipeline.utils import misc as mc
from pyGenicPipeline.core.Input import Input

from csvObject import CsvObject, write_csv
from collections import Counter
from pathlib import Path
from scipy import linalg
import numpy as np


class Score(Input):
    def __init__(self, args):
        super().__init__(args)
        self._score_error_dict = {"Missing Phenotype": 0}
        self._validation = {}

    def construct_chromosome_pgs(self, sm_dict, load_path, chromosome):
        """
        This will construct the pgs from the weights construct with the Gibbs, the infinitesimal or gibbs estimated
        outcomes.
        """

        # Validation we have the necessary information for the scores
        self._assert_construct_pgs()

        # Load the reference to the full sample of ID's
        core = self.gen_reference(load_path)

        # Construct a dict of arrays of our ID's
        ph_dict = self.genetic_phenotypes(core, load_path)

        scores_dict = {h: [] for h in sm_dict.keys() if (self.gibbs in h) or (h == self.inf_dec)}
        chunked_snps = self.chunked_snp_names(sm_dict)
        for index, snp_names in enumerate(chunked_snps):
            print(f"Processing Scores Chunk {index} out of {len(chunked_snps)}")
            # Load the raw snps that have been isolated Raw snps
            raw_snps = self.isolate_raw_snps(core, snp_names)

            # Calculate scores for the infinitesimal model
            self._calculate_score(sm_dict, scores_dict, self.inf_dec, raw_snps, index)

            # Now do the same for each model calculated by Gibbs
            for variant_fraction in sm_dict.keys():
                if self.gibbs in variant_fraction:
                    self._calculate_score(sm_dict, scores_dict, variant_fraction, raw_snps, index)

        scores_dict = {h: mc.flatten(values) for h, values in scores_dict.items()}
        ph_dict = {**ph_dict, **scores_dict}
        # Validate we have all elements of equal length to the number of id's in core
        for key in ph_dict.keys():
            assert len(ph_dict[key]) == core.iid_count

        # Write to file
        rows = [[v[i] for v in ph_dict.values()] for i in range(core.iid_count)]
        write_csv(self.scores_directory, f"Scores_{chromosome}", list(ph_dict.keys()), rows)
        print(f"Finished Constructing scores for Chromosome {chromosome}")

    @staticmethod
    def _calculate_score(sm_dict, score_array, key, raw_snps, index):
        """
        Here we load the weights calculated and re-structure the 1 dimensional list to be a vector array. We then use
        this vector array of beta values from weights alongside the raw snps to calculate the effect of each snp based
        on the nucleotide of the individuals (0, 1 or 2) to compute a score for this key.
        """
        # Restructure weights to be a vector array
        weights = np.array(list([np.array(sm_dict[key]).tolist()])[index])
        weights.shape = (len(sm_dict[key]), 1)

        # Calculate the PRS for the individuals
        score_array[key].append(np.array([np.sum(row) for row in ((-1 * raw_snps) * weights).T]))

    def compile_pgs(self):

        # todo WARNING - STRANGE BEHAVIOUR WHEN COMPILING DIFFERENT FILE TYPES

        # Check we have necessary args for this method
        self._assert_compile_pgs()

        # Get the file names for output from pgs_chromosome_scores
        score_files = mc.directory_iterator(self.scores_directory)

        # Isolate the headers to be aggregated, then aggregate successful scores
        ph_dict, score_keys = self.aggregated_scores(score_files)

        # Load in the raw phenotype, and filter out anyone without one.
        self._load_phenotype(ph_dict)

        # If set, load and filter on any covariant's in the covariates_file
        if self.covariates_file:
            self._load_and_clean_covariants(ph_dict)

        # Calculate the direct effect
        for key in score_keys:
            self._adjust_prs(ph_dict, key, "Direct")

            # Load possible base adjustments, and if we have more than 1 calculation and run all possible combinations
            combinations = [c for c in [self.sex, self.pc, self.covariants] if c in ph_dict.keys()]
            if len(combinations) > 0:
                for combination in mc.possible_combinations(combinations):
                    self._adjust_prs(ph_dict, key, combination, combination)

        # Write out the scores to the working directory
        rows = [[v[i][0] if len(v[i]) == 1 else v[i] for v in ph_dict.values()] for i in range(len(ph_dict[self.iid]))]
        write_csv(self.working_dir, f"Cumulative_Scores", list(ph_dict.keys()), rows)

    def aggregated_scores(self, score_files):
        """
        This will check the headers where we have full information across all our chromosomes and return the names of
        those headers as a list of strings. Failures will be printed to the terminal.

        The successful headers that where isolated will then be aggregated
        """

        # Load the ids, and use this to setup the dict of values
        load_ids = CsvObject(Path(self.scores_directory, score_files[0]), set_columns=True)
        ph_dict = {self.fid: np.array(load_ids[self.fid]), self.iid: np.array(load_ids[self.iid])}

        # Count each header which isn't an id identifier
        headers = Counter(mc.flatten([CsvObject(Path(self.scores_directory, file)).headers for file in score_files]))
        headers.pop(self.fid, None)
        headers.pop(self.iid, None)

        # Isolate headers that are complete
        isolates = {h: np.zeros(len(ph_dict[self.iid])) for h, v in zip(headers.keys(), headers.values())
                    if (v == len(score_files))}

        # Assert we have found any successful headers, and tell users what was dropped, then join our two dicts together
        ec.scores_valid_headers(isolates.keys(), headers, score_files)
        ph_dict = {**ph_dict, **isolates}

        # For each chromosome file
        for file in score_files:
            load_file = CsvObject(Path(self.scores_directory, file), set_columns=True)

            # For each header slice the column from the load file and save it as a float32, then add this to our dict
            for header in isolates.keys():
                ph_dict[header] = np.add(ph_dict[header], np.array(load_file[header], np.float32))

        return ph_dict, list(isolates.keys())

    def _load_phenotype(self, ph_dict):
        """Load the raw phenotype values, and filter anyone out who does have a value in phenotype array"""

        ids = []
        phenotypes = []
        with mc.open_setter(self.phenotype_file)(self.phenotype_file) as file:
            # Determine if we have custom headers or not via _loaded_sum_headers
            headers = [h.lower() for h in file.readline().split()]
            assert self.iid.lower() in headers, "FAILED TO FIND IID"

            id_index = headers.index(self.iid.lower())

            # For each line of phenotype information
            for line in file:
                line = line.split()
                try:
                    # If the phenotype is numeric, not zero, -9 (plink error code, and is a finite value then we will
                    # keep these individuals otherwise we will filter them out.
                    value = float(line[id_index + 1])
                    if (value != 0) and (value != -9) and np.isfinite(value):
                        # We do tell users to put the phenotype next to the iid but maybe change this?
                        phenotypes.append(float(line[id_index + 1]))
                        ids.append(line[id_index])

                # Value may be non-numeric (Ie NA) so skip if a type error occurs
                except TypeError:
                    pass

            file.close()

        # Keep individuals within our phenotype dict if they are within the phenotype file, filter out otherwise
        phenotype_filter = [True if i in ids else False for i in ph_dict[self.iid]]
        mc.filter_array(ph_dict, phenotype_filter)
        self._score_error_dict["Missing Phenotype"] += len(phenotype_filter) - np.sum(phenotype_filter)

        # Construct an array of the phenotypes we found
        phenotypes = np.array(phenotypes)

        # Filter it, as id may exist in phenotype file that is not in the loaded phenotype dict
        phenotype_filter = [True if i in ph_dict[self.iid] else False for i in ids]
        phenotypes = phenotypes[phenotype_filter]
        self._score_error_dict["Missing Phenotype"] += len(phenotype_filter) - np.sum(phenotype_filter)

        # Check everything is of equal length then log the phenotype information to dict
        assert len(Counter([len(v) for v in ph_dict.values()])) == 1, "Filtering on phenotype failed"
        ph_dict[self.phenotype] = phenotypes

        # Reshape our phenotype array to be array to be (N, 1)
        mc.reshape_dict_array(ph_dict, self.phenotype)

    def _load_and_clean_covariants(self, ph_dict):
        """
        This will construct the phenotype dict that we will right out for the end user, storing individual level data to
        help us validate and clean individuals whom we do not have sufficient information to transfer the score too.
        """

        # Load the covariants file
        cov = CsvObject(self.covariates_file, set_columns=True)
        headers = {h.lower(): i for i, h in enumerate(cov.headers)}

        # Load sex if stored in covariates
        if self.sex.lower() in headers.keys():
            ph_dict[self.sex] = np.array(cov.column_data[headers[self.sex.lower()]], dtype=int)

            # Filter out any sex that is not 1 or 2
            sex_filter = np.array([True if (s == 1) or (s == 2) else False for s in ph_dict[self.sex].astype(int)])
            self._score_error_dict["Miss Matched Sex"] = len(sex_filter) - np.sum(sex_filter)
            mc.filter_array(ph_dict, sex_filter)
            mc.reshape_dict_array(ph_dict, self.sex)

        # Load PCs if they exist in the file
        pcs = [headers[h] for h in headers.keys() if h[:2] == self.pc.lower()]
        if len(pcs) > 0:
            self._filter_covariant(ph_dict, self.pc, pcs, cov, "Invalid PCs")

        # If there is anything else, assume it is a covariant
        covariant = [headers[h] for h in headers.keys()
                     if (h[:2] != self.pc.lower()) and (h != self.sex.lower()) and (h not in (self.iid, self.fid))]
        if len(covariant) > 0:
            self._filter_covariant(ph_dict, self.covariants, covariant, cov, "Invalid Covariants")

    def _filter_covariant(self, ph_dict, key, key_indexes, cov, error_dict_key):
        """
        For numeric continuous values we check if they each row contains finite values, otherwise we remove them
        """

        # Save to dict the values from the key
        ph_dict[key] = np.array([[row[i] for i in key_indexes] for row in cov.row_data])

        # Filter out anything that is not finite
        pc_filter = np.array([True if np.isfinite(pc).all() else False for pc in ph_dict[key]])
        self._score_error_dict[error_dict_key] = len(pc_filter) - np.sum(pc_filter)

    def _adjust_prs(self, ph_dict, prs_key, adjust_type, adjust_key=None):
        """
        This will adjust the prs of a given prs key by a given adjustment list provided to adjust_key if it is set,
        otherwise it will just do the direct effect
        """

        # Reshape the current score type
        mc.reshape_dict_array(ph_dict, prs_key)

        # Direct effect
        if not adjust_key:
            (betas, rss00, r, s) = linalg.lstsq(np.ones((len(ph_dict[self.phenotype]), 1)), ph_dict[self.phenotype])

        # Make Adjustment on this key, such as Sex or PC.
        else:
            (betas, rss00, r, s) = linalg.lstsq(np.hstack(self._set_stack(adjust_key, ph_dict, prs_key)))

        x = np.hstack(self._set_stack(adjust_key, ph_dict, prs_key))
        (betas, rss, r, s) = linalg.lstsq(x, ph_dict[self.phenotype])

        # Calculate basic information to print to terminal
        pred_r2 = 1 - rss / rss00
        intercept, effect = betas
        print(F"For {adjust_type}: Predict R2: {pred_r2[0]}, Intercept {intercept[0]}, Effect {effect[0]}\n")

        # Adjust our prs by the betas then store than in ph_dict under this adjust_type
        ph_dict[f"{prs_key}_{adjust_type}"] = np.dot(x, betas)

    def _assert_construct_pgs(self):
        """Assert that the information required to run is present"""
        assert self.ld_radius, ec.missing_arg(self.operation, "LD_Radius")

    def _assert_compile_pgs(self):
        assert self.phenotype, ec.missing_arg(self.operation, "Phenotype")

    def _set_stack(self, adjust_key, ph_dict, prs_key):
        """If covariants are supplied, return a list of base + each value of each k in adjust_keys"""
        base = [ph_dict[prs_key], np.ones((len(ph_dict[self.phenotype]), 1))]
        if adjust_key:
            return base + [ph_dict[k] for k in adjust_key]
        else:
            return base
