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
        self.snps = sids
        self.indexes = indexes
        self.positions = pos
        self.nucleotides = nts
        self.indexed_snps = {snp_id: index for snp_id, index in zip(self.snps, self.indexes)}

    def __repr__(self):
        """
        Human readable printing
        :return: Just the chromosome number
        """
        return f"BimByChromosome {self.chromosome}"


class Nucleotide:
    def __init__(self, nucleotide):
        """
        Sometimes we need to work with saved nucleotides, so we can construct an object that allows us to extract the
        first and second alleles based on what they written as; for example effect allele alt allele.
        :param nucleotide:
        """
        self.a1 = nucleotide[0]
        self.a2 = nucleotide[1]

    def __repr__(self):
        """
        Human readable printing
        :return: String of the tuple of the attributes
        """
        return f"{self.to_tuple()}"

    def to_tuple(self):
        """
        Allow return as tuple
        :return: Tuple of attributes
        """
        return self.a1, self.a2

    def to_list(self):
        """
        Allow return as list
        :return: List of attributes
        """
        return [self.a1, self.a2]
