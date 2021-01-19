from pyGenicPipeline.utils import misc as mc

from miscSupports import load_yaml
from pathlib import Path
import os


class ArgsParser:
    def __init__(self, args):
        """
        Loads args into memory to be used by other Loaders to initialise Input

        Contains common names as attributes so that they can be used as arguments to avoid spelling/alternative naming
        where-ever possible
        """
        self.args = mc.set_args(args)
        self.working_dir = mc.validate_path(self.args["Working_Directory"], False)
        self.config = load_yaml(Path(Path(__file__).parent, "Config.yaml"))

    def make_sub_directory(self, job, name):
        """
        Make a sub-directory within the working directory

        :param job: If this sub directory is part of a job make a job sub directory first, otherwise don't.
        :type job: None | str

        :param name: The name of the end directory you want
        :type name: str
        """
        if job:
            try:
                os.mkdir(Path(self.working_dir, job))
                try:
                    os.mkdir(Path(self.working_dir, job, name))
                except FileExistsError:
                    pass
            except FileExistsError:
                pass
        else:
            try:
                os.mkdir(Path(self.working_dir, name))
            except FileExistsError:
                pass

    @property
    def chromosome(self):
        """Key used for accessing Chromosome headers, groups or other attributes"""
        return "Chromosome"

    @property
    def snp_id(self):
        """Key used for accessing SNP_ID headers, groups or other attributes"""
        return "SNP_ID"

    @property
    def bp_position(self):
        """Key used for accessing bp_position headers, groups or other attributes"""
        return "BP_Position"

    @property
    def effect_allele(self):
        """Key used for accessing Effect_Allele headers, groups or other attributes"""
        return "Effect_Allele"

    @property
    def alt_allele(self):
        """Key used for accessing Alt_Allele headers, groups or other attributes"""
        return "Alt_Allele"

    @property
    def effect_size(self):
        """Key used for accessing Effect_size headers, groups or other attributes"""
        return "Effect_Size"

    @property
    def p_value(self):
        """Key used for accessing P_Value headers, groups or other attributes"""
        return "P_Value"

    @property
    def minor_allele_freq(self):
        """Key used for accessing Minor_Allele_Freq headers, groups or other attributes"""
        return "Minor_Allele_Freq"

    @property
    def info(self):
        """Key used for accessing Minor_Allele_Freq headers, groups or other attributes"""
        return "Info"

    @property
    def standard_errors(self):
        """Key used for accessing Minor_Allele_Freq headers, groups or other attributes"""
        return "Standard_Errors"

    @property
    def case_freq(self):
        """Key used for accessing Case_Freq headers, groups or other attributes"""
        return "Case_Freq"

    @property
    def case_n(self):
        """Key used for accessing Case_N headers, groups or other attributes"""
        return "Case_N"

    @property
    def control_freq(self):
        """Key used for accessing Control_Freq headers, groups or other attributes"""
        return "Control_Freq"

    @property
    def control_n(self):
        """Key used for accessing Control_N headers, groups or other attributes"""
        return "Control_N"

    @property
    def sm_lines(self):
        """Key for accessing lines in the summary stats file"""
        return "SM_Lines"

    @property
    def sm_variants(self):
        """Key for accessing variants associated with the snps found in the summary stats file"""
        return "SM_Variants"

    @property
    def beta(self):
        """Key used for accessing Beta headers, groups or other attributes"""
        return "Beta"

    @property
    def log_odds(self):
        """Key used for accessing Log_Odds headers, groups or other attributes"""
        return "Log_Odds"

    @property
    def nucleotide(self):
        """Key used for accessing Nucleotide headers, groups or other attributes"""
        return "Nucleotide"

    @property
    def freq(self):
        """Key used for accessing Frequency headers, groups or other attributes"""
        return "Frequency"

    @property
    def stds(self):
        """Key used for accessing Standard_Deviations headers, groups or other attributes"""
        return "Standard_Deviations"

    @property
    def raw_snps(self):
        """Key used for accessing Raw_Snps headers, groups or other attributes"""
        return "Raw_Snps"

    @property
    def filter_key(self):
        """Key used for accessing a Filter in headers, groups or other attributes"""
        return "Filter"

    @property
    def norm_snps(self):
        """Key used for accessing Raw_Snps headers, groups or other attributes"""
        return "Normalised_Snps"

    @property
    def count_snp(self):
        """Key used for accessing the number of snps in headers, groups or other attributes"""
        return "Snp_Count"

    @property
    def count_iid(self):
        """Key used for accessing the number of individuals in headers, groups or other attributes"""
        return "IID_Count"

    @property
    def avg_ld(self):
        """Key used for accessing Average LD in headers, groups or other attributes"""
        return "Avg_LD"

    @property
    def herit(self):
        """Key used for accessing Heritability in headers, groups or other attributes"""
        return "Heritability"

    @property
    def ld_scores(self):
        """Key used for accessing LD_Scores in headers, groups or other attributes"""
        return "LD_Scores"

    @property
    def ld_dict(self):
        """Key used for accessing LD_Dict in headers, groups or other attributes"""
        return "LD_Dict"

    @property
    def gibbs_beta(self):
        """Key used for accessing Gibbs_Beta in headers, groups or other attributes"""
        return "Gibbs_Beta"

    @property
    def genome_key(self):
        """Key used for accessing Genome wide data in headers, groups or other attributes"""
        return "Genome"

    @property
    def inf_dec(self):
        """Key used for accessing Infinitesimal data in headers, groups or other attributes"""
        return "Infinitesimal"

    @property
    def gibbs(self):
        """Key used for accessing Gibbs data in headers, groups or other attributes"""
        return "Gibbs"

    @property
    def fam(self):
        """Key used for accessing the FamID Variants in headers, groups or other attributes"""
        return "Family"

    @property
    def iid(self):
        """Key used for accessing the Individual Identifiers in headers, groups or other attributes"""
        return "IID"

    @property
    def fid(self):
        """Key used for accessing the Family Identifiers in headers, groups or other attributes"""
        return "FID"

    @property
    def f_id(self):
        """Key used for accessing the Fathers Identifiers in headers, groups or other attributes"""
        return "F_ID"

    @property
    def m_id(self):
        """Key used for accessing the Mothers Identifiers in headers, groups or other attributes"""
        return "F_ID"

    @property
    def sex(self):
        """Key used for accessing the Sex in headers, groups or other attributes"""
        return "Sex"

    @property
    def pc(self):
        """Key used for accessing the Principle components in headers, groups or other attributes"""
        return "PC"

    @property
    def phenotype(self):
        """Key used for accessing the raw phenotype in headers, groups or other attributes"""
        return "Phenotype"

    @property
    def covariants(self):
        """Key used for accessing the covariants in headers, groups or other attributes"""
        return "Covariants"
