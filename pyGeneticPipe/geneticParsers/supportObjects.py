class Variant:
    __slots__ = ["chromosome", "bp_position", "variant_id", "a1", "a2"]

    def __init__(self, chromosome, bp_position, variant_id, a1, a2):
        """
        Contains variant information for bgen files, unlike Bim does not contain morgan position
        """
        self.chromosome = chromosome
        self.bp_position = bp_position
        self.variant_id = variant_id
        self.a1 = a1
        self.a2 = a2

    def __repr__(self):
        return f"Variant {self.variant_id} - {self.chromosome}:{self.bp_position} A1{self.a1}:A2{self.a2}"


class BimVariant:
    __slots__ = ["chromosome", "variant_id", "morgan_pos", "bp_position", "a1", "a2"]

    def __init__(self, chromosome, variant_id, morgan_pos, bp_position, a1, a2):
        """
        Holds information of a bim file in Object notation. Bim files contain a set 6 categories. Full specification:

        https://www.cog-genomics.org/plink2/formats#bim
        """
        self.chromosome = chromosome
        self.variant_id = variant_id
        self.morgan_pos = morgan_pos
        self.bp_position = bp_position
        self.a1 = a1
        self.a2 = a2

    def __repr__(self):
        return f"Variant {self.variant_id} - {self.chromosome}:{self.bp_position} A1{self.a1}:A2{self.a2}"

    def to_variant(self):
        """
        Bgen doesn't use morgan position, and to allow for standard variant operations rather than having optional
        attributes, this will return the Variant Class containing all the attributes of the bim minus the morgan_pos.

        :return: Variant
        :rtype: Variant
        """
        return Variant(self.chromosome, self.bp_position, self.variant_id, self.a1, self.a2)


class BimLoci:
    def __init__(self, bim_line):
        """
        Holds information of a bim file in Object notation. Bim files contain a set 6 categories. Full specification:

        https://www.cog-genomics.org/plink2/formats#bim

        """
        chromosome, variant_id, morgan_pos, bp_position, a1, a2 = bim_line.split()

        self.chromosome = chromosome
        self.variant_id = variant_id
        self.morgan_pos = float(morgan_pos)
        self.bp_position = int(bp_position)
        self.a1 = str(a1)
        self.a2 = str(a2)

    def __repr__(self):
        """
        Allow for human readable printing

        :return: Each column in order printed as a string warped in square brackets
        """
        return f"[{self.chromosome}, {self.variant_id}, {self.morgan_pos}, {self.bp_position}, {self.a1}, {self.a2}]"


class BimByChromosome:
    def __init__(self, chromosome, variant_ids, indexes, bp_positions, nucleotides):
        """
        Sometimes we may need to have bim information by chromosome, this will contain all the information required from
        the bim files, but for a given chromosome.
        """
        self.chromosome = chromosome
        self.snps = variant_ids  # Will still need this as we want to check for common snps across types
        self.indexes = indexes
        self.positions = bp_positions
        self.nucleotides = nucleotides
        self.indexed_snps = {snp_id: index for snp_id, index in zip(self.snps, self.indexes)}

        self.snp_information = {sids: SNPIndex(sids, i, position, nts)
                                for sids, i, position, nts in zip(variant_ids, indexes, bp_positions, nucleotides)}

    def __repr__(self):
        """
        Human readable printing
        :return: Just the chromosome number
        """
        return f"BimByChromosome {self.chromosome}"

    def extract_snp(self, variant):
        """
        Return a snp from snp information

        :param variant: A snp, such as rs3131962
        :return: A SNPIndex containing the information for that snp
        """
        return self.snp_information[variant]


class SNPIndex:
    def __init__(self, sids, index, bp_position, nts):

        self.snp = sids
        self.index = index
        self.bp_position = bp_position
        self.nucleotide = Nucleotide(nts)

    def __repr__(self):
        return f"{self.snp}: Index: {self.index} BP_Position: {self.bp_position} Nucleotide: {self.nucleotide}"


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
