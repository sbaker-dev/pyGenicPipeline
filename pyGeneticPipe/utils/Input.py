from pyGeneticPipe.utils import error_codes as ec
from pyGeneticPipe.utils import misc as mc
from pathlib import Path


class Input:
    def __init__(self, args):

        self._args = args
        self._frequencies = args["Summary_Frequency"]
        self._mandatory_headers = ["SNP_ID", "Effect_Allele", "Alt_Allele", "Effect_size", "P_Value"]
        self._summary_headers = self._set_summary_headers()

        self.debug = args["Debug"]
        self.ld_ref_mode, self.bgen, self.bed, self.bim, self.fam = self._set_ld_ref(args["LD_Reference_Genotype"])
        self.summary_stats = self._set_summary_stats(args["Summary_Stats"])

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
        if not summary_stats_path:
            return None

        # Construct path as an object
        summary_path = Path(summary_stats_path)

        # Check file exists
        assert summary_path.exists(), ec.path_invalid(summary_path, "_set_summary_stats")

        # Determine if file is g-zipped
        gz_status = (summary_path.suffix == ".gz")

        with mc.open_setter(summary_path)(summary_path) as file:
            raw_headers = file.readline()
            # Determine if we have custom headers or not via _summary_headers
            headers = {header: self._check_header(header, mc.decode_line(raw_headers, gz_status))
                       for header in self._summary_headers}

            if self._frequencies:
                self._sum_stats_frequencies()
            else:
                self._sum_stats(file, gz_status, headers)

        return 0

    def _sum_stats(self, file, gz_status, headers):
        for index, line in enumerate(file):
            print(index, line)
            line = mc.decode_line(line, gz_status)
            print(line)
            print(line[headers["SNP_ID"]])
            print(self.fam)

            break

    def _sum_stats_frequencies(self):
        raise NotImplementedError("Frequencies are not yet implemented")

    def _check_header(self, sum_header, headers):
        """
        We need to standardise our headers, and locate where the current header is in our summary file in terms of a
        base zero index.

        :param sum_header: Standardised header
        :param headers: summary statistics headers
        :return: None if not found else the index of the header in our file for this standardised header
        :rtype: None | int
        """
        header_indexes = [i for i, h in enumerate(headers) if h in self._summary_headers[sum_header]]

        assert len(header_indexes) < 2, ec.ambiguous_header(sum_header, headers, self._summary_headers[sum_header])
        if len(header_indexes) == 0:
            assert sum_header not in self._mandatory_headers, ec.mandatory_header(
                sum_header, headers, self._summary_headers[sum_header])
            return None
        else:
            return header_indexes[0]

    def _set_summary_headers(self):
        """
        We may have users using custom headers, or they may be using a format we already have covered

        Note
        -----
        In GUI we basically want to call the header check to make sure we align columns correctly. If not they can set
        it themselves

        :return: The headers to validate
        """
        if self._args["Summary_Headers"]:
            # Recast so that the values are in a list so they can be checked by the same method as defaults
            header_sets = self._args["Summary_Headers"]
            return {key: [v] for key, v in zip(header_sets.keys(), header_sets.values())}
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
                "Info": ["Info", "info"]
            }
            return header_sets
