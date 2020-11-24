class Variant:
    __slots__ = ["chromosome", "bp_position", "variant_id", "a1", "a2"]

    def __init__(self, chromosome, bp_position, variant_id, a1, a2):
        """
        Contains variant information for bgen files, unlike Bim does not contain morgan position
        """
        self.chromosome = chromosome
        self.bp_position = int(bp_position)
        self.variant_id = variant_id
        self.a1 = a1
        self.a2 = a2

    def __repr__(self):
        return f"{self.variant_id} - CHR{self.chromosome} - POS{self.bp_position} - A1-{self.a1}:A2-{self.a2}"

    def nucleotide(self, as_list=False):
        """Return a tuple of a1 a2 as the nucleotide"""
        if as_list:
            return [self.a1, self.a2]
        else:
            return self.a1, self.a2

    def bgen_variant_id(self):
        """Bgen for pysnptools requires the variant and rs-id but in this case we just submit the same for both"""
        return f"{self.variant_id},{self.variant_id}"


class BimVariant:
    __slots__ = ["chromosome", "variant_id", "morgan_pos", "bp_position", "a1", "a2"]

    def __init__(self, chromosome, variant_id, morgan_pos, bp_position, a1, a2):
        """
        Holds information of a bim file in Object notation. Bim files contain a set 6 categories. Full specification:

        https://www.cog-genomics.org/plink2/formats#bim
        """
        self.chromosome = chromosome
        self.variant_id = variant_id
        self.morgan_pos = float(morgan_pos)
        self.bp_position = int(bp_position)
        self.a1 = a1
        self.a2 = a2

    def __repr__(self):
        return f"{self.variant_id} - CHR{self.chromosome} - BP:{self.bp_position}/MP:{self.morgan_pos} - " \
               f"A1-{self.a1}:A2-{self.a2}"

    def to_variant(self):
        """
        Bgen doesn't use morgan position, and to allow for standard variant operations rather than having optional
        attributes, this will return the Variant Class containing all the attributes of the bim minus the morgan_pos.

        :return: Variant
        :rtype: Variant
        """
        return Variant(self.chromosome, self.bp_position, self.variant_id, self.a1, self.a2)


class SMVariant:
    __slots__ = ["chromosome", "variant_id", "bp_position", "a1", "a2", "beta", "beta_odds", "p_value", "info",
                 "frequency"]

    def __init__(self, chromosome, variant_id, bp_position, a1, a2, beta, beta_odds, p_value, info, frequency):
        """
        Contains Summary variant information
        """
        self.chromosome = chromosome
        self.variant_id = variant_id
        self.bp_position = bp_position
        self.a1 = a1
        self.a2 = a2
        self.beta = beta
        self.beta_odds = beta_odds
        self.p_value = p_value
        self.info = info
        self.frequency = frequency


class Nucleotide:
    def __init__(self, a1, a2):
        """
        Nucleotides, so we can construct an object that allows us to extract the first and second alleles based on what
        they written as; for example effect allele alt allele.
        """
        self.a1 = a1
        self.a2 = a2

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
