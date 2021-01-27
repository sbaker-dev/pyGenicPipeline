from pyGenicPipeline.utils import errors as ec
from pyGenicPipeline.utils import misc as mc
from ..argsParser import ArgsParser

from miscSupports import string_to_bool
from csvObject import write_csv
from pathlib import Path


class SummaryLoader(ArgsParser):
    def __init__(self, args):
        """
        These arguments are specifically to assist with PGS SummaryCLeaner
        """
        super().__init__(args)
        self._config = self.config["Summary"]

        # Set the write directory for Summary Files
        self.make_sub_directory("PGS", "Cleaned")
        self.summary_directory = Path(self.working_dir, "PGS", "Cleaned")

        # Validate the summary file, and set the load type, sample size and header indexing
        self.summary_file = mc.validate_path(self.args["Summary_Path"])
        self.zipped, self.sample_size = self._set_summary_stats()
        self._summary_headers = self._set_summary_headers()

        # Determine effect type, and load the configuration of alleles, set z scores if requried
        self.effect_type = self._set_effect_type(self.args["Summary_Effect_Type"])
        self.ambiguous_snps, self.allowed_alleles, self.allele_flip = self._configure_alleles()
        self.z_scores = self._set_z_scores(self.args["Z_Scores"])

        self.summary_headers, self.summary_dict = self._set_summary_write_headers()

    def _set_summary_stats(self):
        """
        If we are reading in the summary statistics file then validate its path, construct a valid set of snps in a set
        and map the chromosome and bp_position to each valid snp in a dict. Validate the type of zipped structure and
        the sample size of the study

        :return: summary_path, snp_pos_map, valid_snp, gz_status, sample_size
        """

        # Stop if not required
        if not self.summary_file:
            return None, None

        # Check the sample size from this study has been provided
        sample_size = self.args["Summary_Sample_Size"]
        assert (sample_size is not None) and (sample_size > 0), ec.sample_size()

        # Determine if the summary file is g-zipped
        gz_status = (self.summary_file.suffix == ".gz")
        return gz_status, sample_size

    def _set_summary_headers(self):
        """
        We may have users using custom headers, or they may be using a format we already have covered

        Note
        -----
        In GUI we basically want to call the header check to make sure we align columns correctly. If not they can set
        it themselves

        :return: The headers to validate
        """
        if not self.summary_file:
            return None

        custom_headers = self.args["Custom_Summary_Header"]
        if custom_headers:
            # Recast so that the values are in a list so they can be checked by the same method as defaults
            header_sets = {key: [v] for key, v in zip(custom_headers.keys(), custom_headers.values())}
        else:
            # Based on known summary statistics from LDPred sum_stats_parsers.py
            header_sets = self._config["header_keys"]

        with mc.open_setter(self.summary_file)(self.summary_file) as file:

            # Determine if we have custom headers or not via _loaded_sum_headers
            raw_headers = file.readline()

            headers = {header: self._check_header(header, mc.decode_line(raw_headers, self.zipped), header_sets)
                       for header in header_sets}

            file.close()
            return headers

    def _check_header(self, sum_header, headers, summary_headers):
        """
        We need to standardise our headers, and locate where the current header is in our summary file in terms of a
        base zero index.

        :param sum_header: Standardised header
        :param headers: summary statistics headers

        :return: None if not found else the index of the header in our file for this standardised header
        :rtype: None | int
        """
        header_indexes = [i for i, h in enumerate(headers) if h in summary_headers[sum_header]]

        assert len(header_indexes) < 2, ec.ambiguous_header(sum_header, headers, summary_headers[sum_header])
        if len(header_indexes) == 0:
            assert sum_header not in self._config["Mandatory_Headers"], ec.mandatory_header(
                sum_header, headers, summary_headers[sum_header])
            return None
        else:
            return header_indexes[0]

    def _set_effect_type(self, effect_type):
        """
        Set the effect type of the betas for GWAS summary stats if set
        """
        if effect_type:
            assert effect_type in self._config["effect_types"], ec.invalid_effect_type(
                self._config["effect_types"], effect_type)
            return effect_type
        else:
            return None

    def _configure_alleles(self):
        """
        Yaml storage of tuples/ sets didn't work so this configures a list of lists into a set of tuples for ambiguous,
        sets a set of allow alleles, and also returns the dict of allele_flip
        """
        ambiguous = set(tuple(ambiguous) for ambiguous in self._config["ambiguous_snps"])
        return ambiguous, set(self._config["allowed_alleles"]), self._config["allele_flip"]

    def _set_z_scores(self, set_z_scores):
        """
        If the user wants to compute z scores, then standard_errors most be set but otherwise it isn't a mandatory
        header.

        :param set_z_scores: A bool of if z scores should be calculated or not
        :type set_z_scores: bool

        :return: True if assertion of standard errors is also True, if set_z_scores is None then return None
        :rtype: None | bool
        """
        if set_z_scores is None:
            set_z_scores = False

        if string_to_bool(set_z_scores):
            assert self._summary_headers[self.standard_errors] is not None, ec.z_scores_with_standard_errors
            return True
        else:
            return None

    def _set_summary_write_headers(self):
        """Construct headers to be used for writing and reading cleaned files"""
        write_headers = [self.chromosome, self.bp_position, self.snp_id, self.effect_allele, self.alt_allele,
                         self.log_odds, self.beta, self.freq]

        write_dict = {header: i for i, header in enumerate(write_headers)}
        return write_headers, write_dict

    def write_summary_files(self, sm_dict, write, chromosome, summary_type, directory):
        """Write out the information from sm_dict into a csv"""
        if write:
            rows_out = []
            for v, log_odds, beta, freq, in zip(sm_dict[self.sm_variants], sm_dict[self.log_odds], sm_dict[self.beta],
                                                sm_dict[self.freq]):
                rows_out.append(v.items() + [log_odds, beta, freq])

            write_csv(directory, f"{summary_type}_{chromosome}", self.summary_headers, rows_out)
        return sm_dict

    @property
    def cleaned_types(self):
        """The types of each column in the cleaned data"""
        return [int, int, str, str, str, float, float, float]

    @property
    def sm_chromosome(self):
        """Summary stat chromosome header index in GWAS summary file"""
        return self._summary_headers[self.chromosome]

    @property
    def sm_snp_id(self):
        """Snp/variant id header index in GWAS summary file"""
        return self._summary_headers[self.snp_id]

    @property
    def sm_effect_allele(self):
        """Effect allele header index in GWAS summary file"""
        return self._summary_headers[self.effect_allele]

    @property
    def sm_alt_allele(self):
        """Alt allele header index in GWAS summary file"""
        return self._summary_headers[self.alt_allele]

    @property
    def sm_bp_position(self):
        """Base pair position header index in GWAS summary file"""
        return self._summary_headers[self.bp_position]

    @property
    def sm_p_value(self):
        """P value position header index in GWAS summary file"""
        return self._summary_headers[self.p_value]

    @property
    def sm_effect_size(self):
        """Effect size header index in GWAS summary file"""
        return self._summary_headers[self.effect_size]

    @property
    def sm_minor_allele_freq(self):
        """Minor allele Frequency header index in GWAS summary file"""
        return self._summary_headers[self.minor_allele_freq]

    @property
    def sm_info(self):
        """Info score header index in GWAS summary file"""
        return self._summary_headers[self.info]

    @property
    def sm_standard_errors(self):
        """Standard error header index in GWAS summary file"""
        return self._summary_headers[self.standard_errors]

    @property
    def sm_case_freq(self):
        """Case_Freq header index in GWAS summary file"""
        return self._summary_headers[self.case_freq]

    @property
    def sm_case_n(self):
        """Case_N header index in GWAS summary file"""
        return self._summary_headers[self.case_freq]

    @property
    def sm_control_freq(self):
        """Control_Freq header index in GWAS summary file"""
        return self._summary_headers[self.case_freq]

    @property
    def sm_control_n(self):
        """Control_N header index in GWAS summary file"""
        return self._summary_headers[self.case_freq]
