class BimLoci:
    def __init__(self, bim_line):
        """
        Holds information of a bim file in Object notation. Bim files contain a set 6 categories. Full specification:

        https://www.cog-genomics.org/plink2/formats#bim

        """
        self.chromosome, self.variant_id, self.morgan_pos, self.bp_position, self.a1, self.a2 = bim_line.split()

    def __repr__(self):
        """
        Allow for human readable printing

        :return: Each column in order printed as a string warped in square brackets
        """
        return f"[{self.chromosome}, {self.variant_id}, {self.morgan_pos}, {self.bp_position}, {self.a1}, {self.a2}]"


class BimByChromosome:
    def __init__(self, chromosome, sids, indexes, pos, nts):
        """
        Sometimes we may need to have bim information by chromosome, this will contain all the information required from
        the bim files, but for a given chromosome.
        """
        self.chromosome = chromosome
        self.chromosome_snps = sids
        self.indexes = indexes
        self.chromosome_positions = pos
        self.nucleotides = nts

    def __repr__(self):
        """
        Human readable printing
        :return: Just the chromosome number
        """
        return f"BimByChromosome {self.chromosome}"
