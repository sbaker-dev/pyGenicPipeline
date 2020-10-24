from pyGeneticPipe.utils import error_codes as ec
from pyGeneticPipe.utils import misc as mc
from pathlib import Path
import pickle
import gzip


class Input:
    def __init__(self, args):

        self.args = args
        self.debug = args["Debug"]

        # Summary setters
        self._hap_map_3 = self._set_hap_map_3()
        self._mandatory_headers = ["SNP_ID", "Effect_Allele", "Alt_Allele", "Effect_size", "P_Value"]
        self._effect_types = ["OR", "LINREG", "LOGOR", "BLUP"]

        self.ld_ref_mode, self.bgen, self.bed, self.bim, self.fam = self._set_ld_ref(args["LD_Reference_Genotype"])
        self.summary_path, self.snp_map, self.valid_snps, self.zipped, self.sample_size, self.valid_chromosomes = \
            self._set_summary_stats(args["Summary_Stats"])
        self._summary_headers = self._set_summary_headers(args["Summary_Headers"], args["Summary_Stats"])
        self.frequencies = args["Summary_Frequency"]
        self.effect_type = self._set_effect_type(args["Effect_type"])
        self.z_scores = self._set_z_scores(args["Z_Scores"])

    @staticmethod
    def _set_ld_ref(ref_path):
        """
        When cleaning summary statistics we need an ld-reference-genotype file to do so. This method will attempt to
        load ethier a .bgen file and return it as an object or load a .bed, .bim, and .fam plink file.

        :param ref_path: path to ld_ref_file
        :type ref_path: None | str

        :return: The mode we are working in (bgen or plink), the bgen object if loaded else None, and the three plink
            files that where load else None.
        """
        # If there was no path for ld_ref, then just return None for all 5 parameters
        if not ref_path:
            return None, None, None, None, None

        # Construct path as an object
        ld_path = Path(ref_path)

        # Check file home directory can be reached
        assert ld_path.parent.exists(), ec.path_invalid(ld_path.parent, "_set_ld_ref")

        # Check the mode we are running in, and return args of mode, bgenObject and the 3 plink files accordingly
        if ld_path.suffix == ".bgen":
            # return "bgen", bgenObject(ld_path), None, None, None
            raise NotImplementedError("Reading bgen files not yet implemented")

        # If a file has a plink suffix take the stem of the name otherwise just take the name
        if ld_path.suffix == (".bed" or ".bim" or ".fam"):
            bed = Path(f"{str(ld_path.parent)}/{ld_path.stem}.bed")
            bim = Path(f"{str(ld_path.parent)}/{ld_path.stem}.bim")
            fam = Path(f"{str(ld_path.parent)}/{ld_path.stem}.fam")
        else:
            bed = Path(f"{str(ld_path.parent)}/{ld_path.name}.bed")
            bim = Path(f"{str(ld_path.parent)}/{ld_path.name}.bim")
            fam = Path(f"{str(ld_path.parent)}/{ld_path.name}.fam")

        # Check the files exists then return with mode of plink, no bgen object and a bed, bim and fam file.
        assert bed.exists(), ec.path_invalid(bed, "_set_ld_ref")
        assert bim.exists(), ec.path_invalid(bim, "_set_ld_ref")
        assert fam.exists(), ec.path_invalid(fam, "_set_ld_ref")
        return "plink", None, bed, bim, fam

    def _set_summary_stats(self, summary_stats_path):
        """
        If we are reading in the summary statistics file then validate its path, construct a valid set of snps in a set
        and map the chromosome and bp_position to each valid snp in a dict. Validate the type of zipped structure and
        the sample size of the study
        :param summary_stats_path: Path to the summary stats file or None if we don't need to read one in
        :type summary_stats_path: str | None

        :return: summary_path, snp_pos_map, valid_snp, gz_status, sample_size
        """

        # Stop if not required
        if not summary_stats_path:
            return None

        # Construct path as an object and that check it exists and the sample size from this study has been provided
        summary_path = Path(summary_stats_path)
        sample_size = self.args["Sample_Size"]
        assert summary_path.exists(), ec.path_invalid(summary_path, "_set_summary_stats")
        assert (sample_size is not None) and (sample_size > 0), ec.sample_size()

        # Construct the valid snp list
        snp_pos_map, valid_snp, valid_chromosomes = self._validation_snp_list()

        # Determine if the summary file is g-zipped
        gz_status = (summary_path.suffix == ".gz")
        return summary_path, snp_pos_map, valid_snp, gz_status, sample_size, valid_chromosomes

    def _validation_snp_list(self):
        """
        Create a set of valid snps and a position map to them via plink or bgen with option censuring of chromosome and
        snps via HapMap3

        :return: A dict of the snp_id: morgan positioning and chromosome
        """
        if self.ld_ref_mode == "plink":
            accepted_chromosomes = self._set_accepted_chromosomes()
            snp_pos_map = {}
            valid_snp = set()
            valid_chromosomes = set()

            with open(self.bim) as f:
                for line in f:
                    chromosome, variant_id, morgan_pos, bp_cord, a1, a2 = line.split()
                    # If the user has specified certain chromosomes check that this snps chromosome is in accepted_list
                    if accepted_chromosomes and (chromosome not in accepted_chromosomes):
                        continue

                    # If the user has specified only to use snp id's from HapMap3 then check this condition
                    if self._hap_map_3 and (variant_id in self._hap_map_3):
                        valid_snp.add(variant_id)
                        snp_pos_map[variant_id] = {"Position": bp_cord, "Chromosome": chromosome}
                        valid_chromosomes.add(chromosome)

                    # Otherwise add to valid snps / snp_pos_map
                    else:
                        valid_snp.add(variant_id)
                        snp_pos_map[variant_id] = {"Position": bp_cord, "Chromosome": chromosome}
                        valid_chromosomes.add(chromosome)

            assert len(valid_snp) > 0, ec.no_valid_snps(self.bim, accepted_chromosomes, self._hap_map_3)
            return snp_pos_map, valid_snp, valid_chromosomes

        elif self.ld_ref_mode == "bgen":
            raise NotImplementedError("Bgen mode not yet implemented")

        else:
            raise Exception(f"CRITICAL ERROR: ld_ref_mode takes the value 'plink' or 'bgen' yet found"
                            f" {self.ld_ref_mode}")

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
            assert sum_header not in self._mandatory_headers, ec.mandatory_header(
                sum_header, headers, summary_headers[sum_header])
            return None
        else:
            return header_indexes[0]

    def _set_summary_headers(self, header_arg, construct_summary_stats):
        """
        We may have users using custom headers, or they may be using a format we already have covered

        Note
        -----
        In GUI we basically want to call the header check to make sure we align columns correctly. If not they can set
        it themselves

        :return: The headers to validate
        """
        if not construct_summary_stats:
            return None

        if header_arg:
            # Recast so that the values are in a list so they can be checked by the same method as defaults
            header_sets = {key: [v] for key, v in zip(header_arg.keys(), header_arg.values())}
        else:
            # Based on known summary statistics from LDPred sum_stats_parsers.py
            header_sets = {
                "Chromosome": ["CHR", "chr", "hg19chrc"],
                "SNP_ID": ["SNP_ID", "rs", "snpid", "MarkerName", "SNP"],
                "Effect_Allele": ["REF", "ref", "a1", "A1", "Allele1"],
                "Alt_Allele": ["ALT", "alt", "a2", "A2", "Allele2"],
                "Position": ["POS", "pos", "bp", "BP"],
                "P_Value": ["PVAL", "pval", "p", "P", "P.2gc"],
                "Effect_size": ["LINREG", "OR"],
                "Minor_Allele_Freq": ["REF_FRQ", "reffrq", "Freq.Hapmap.Ceu", "Freq.Allele1.HapMapCEU", "MAF"
                                      "Freq.Hapmap.Ceu"],
                "Info": ["Info", "info"],
                "Standard_Errors": ["SE", "se", "SE.2gc"]
            }

        with mc.open_setter(self.summary_path)(self.summary_path) as file:

            # Determine if we have custom headers or not via _loaded_sum_headers
            raw_headers = file.readline()
            headers = {header: self._check_header(header, mc.decode_line(raw_headers, self.zipped), header_sets)
                       for header in header_sets}

            file.close()
            return headers

    def _set_accepted_chromosomes(self):
        """
        Individuals may wish to censure Y or Mitochondrial genes, and can do so via Custom_Chromosome. Else all
        chromosomes are assumed to be valid

        :return: A list of accept chromosome names
        """
        if self.args["Custom_Chromosome"]:
            return self.args["Custom_Chromosome"]
        else:
            return None

    def _set_hap_map_3(self):
        """
        Users may wish to limit valid snps to those found within HapMap3. If they do, they need to provide a path to the
        hapmap3 snp file which will be check that it exists, have the snps extracted and return. Otherwise set to none
        :return: The valid HapMap3 snps or None
        """

        if self.args["Only_HapMap3"]:
            # Construct path as an object and check it exists
            hap_map_path = Path(self.args["Only_HapMap3"])
            assert hap_map_path.exists(), ec.path_invalid(hap_map_path, "_set_hap_map_3")

            # If the HapMap3 file exists, then extract the snp ids and return them
            f = gzip.open(hap_map_path, 'r')
            hm3_sids = pickle.load(f)
            f.close()
            return hm3_sids
        else:
            return None

    def _set_effect_type(self, effect_type):
        """
        Set the effect type of the betas for GWAS summary stats if set
        """
        if effect_type:
            assert effect_type in self._effect_types, ec.invalid_effect_type(self._effect_types, effect_type)
            return effect_type
        else:
            return None

    def _set_z_scores(self, set_z_scores):
        """
        If the user wants to compute z scores, then standard_errors most be set but otherwise it isn't a mandatory
        header.
        :param set_z_scores: A bool of if z scores should be calculated or not
        :type set_z_scores: bool

        :return: True if assertion of standard errors is also True, if set_z_scores is None then return None
        :rtype: None | bool
        """
        if set_z_scores:
            assert self.standard_errors is not None, ec.z_scores_with_standard_errors
            return True
        else:
            return None

    @property
    def chromosome(self):
        """Summary stat chromosome header index in GWAS summary file"""
        return self._summary_headers["Chromosome"]

    @property
    def snp_id(self):
        """Snp/variant id header index in GWAS summary file"""
        return self._summary_headers["SNP_ID"]

    @property
    def effect_allele(self):
        """Effect allele header index in GWAS summary file"""
        return self._summary_headers["Effect_Allele"]

    @property
    def alt_allele(self):
        """Alt allele header index in GWAS summary file"""
        return self._summary_headers["Alt_Allele"]

    @property
    def bp_position(self):
        """Base pair position header index in GWAS summary file"""
        return self._summary_headers["Position"]

    @property
    def p_value(self):
        """P value position header index in GWAS summary file"""
        return self._summary_headers["P_Value"]

    @property
    def effect_size(self):
        """Effect size header index in GWAS summary file"""
        return self._summary_headers["Effect_size"]

    @property
    def minor_allele_freq(self):
        """Minor allele Frequency header index in GWAS summary file"""
        return self._summary_headers["Minor_Allele_Freq"]

    @property
    def info(self):
        """Info score header index in GWAS summary file"""
        return self._summary_headers["Info"]

    @property
    def standard_errors(self):
        """Standard error header index in GWAS summary file"""
        return self._summary_headers["Standard_Errors"]
