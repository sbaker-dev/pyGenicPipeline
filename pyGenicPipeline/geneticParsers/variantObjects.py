class Variant:
    __slots__ = ["chromosome", "bp_position", "snp_id", "a1", "a2"]

    def __init__(self, chromosome, bp_position, snp_id, a1, a2):
        """
        Contains variant information for bgen files, unlike Bim does not contain morgan position
        """
        self.chromosome = str(int(chromosome))
        self.bp_position = int(bp_position)
        self.snp_id = snp_id
        self.a1 = a1
        self.a2 = a2

    def __repr__(self):
        """Human Readable print"""
        return f"{self.snp_id} - CHR{self.chromosome} - POS{self.bp_position} - A1-{self.a1}:A2-{self.a2}"

    def __getitem__(self, item):
        """Get an item from Variant"""
        return getattr(self, item)

    def items(self):
        """Return all the items in a list for this variant"""
        return [self.chromosome, self.bp_position, self.snp_id, self.a1, self.a2]

    def nucleotide(self, as_list=False):
        """Return a tuple of a1 a2 as the nucleotide"""
        if as_list:
            return [self.a1, self.a2]
        else:
            return self.a1, self.a2

    def bgen_snp_id(self):
        """Bgen for pysnptools requires the variant and rs-id but in this case we just submit the same for both"""
        return f"{self.snp_id},{self.snp_id}"


class BimVariant:
    __slots__ = ["chromosome", "snp_id", "morgan_pos", "bp_position", "a1", "a2"]

    def __init__(self, chromosome, snp_id, morgan_pos, bp_position, a1, a2):
        """
        Holds information of a bim file in Object notation. Bim files contain a set 6 categories. Full specification:

        https://www.cog-genomics.org/plink2/formats#bim
        """
        self.chromosome = chromosome
        self.snp_id = snp_id
        self.morgan_pos = float(morgan_pos)
        self.bp_position = int(bp_position)
        self.a1 = a1
        self.a2 = a2

    def __repr__(self):
        return f"{self.snp_id} - CHR{self.chromosome} - BP:{self.bp_position}/MP:{self.morgan_pos} - " \
               f"A1-{self.a1}:A2-{self.a2}"

    def to_variant(self):
        """
        Bgen doesn't use morgan position, and to allow for standard variant operations rather than having optional
        attributes, this will return the Variant Class containing all the attributes of the bim minus the morgan_pos.

        :return: Variant
        :rtype: Variant
        """
        return Variant(self.chromosome, self.bp_position, self.snp_id, self.a1, self.a2)


class FamId:
    __slots__ = ["fid", "iid", "f_id", "m_id", "sex", "phenotype"]

    def __init__(self, fid, iid, father_id, mother_id, sex, phenotype):
        """
        Holds information About family Identifiers
        """
        self.fid = fid
        self.iid = iid
        self.f_id = father_id
        self.m_id = mother_id
        self.sex = sex
        self.phenotype = phenotype

    def __repr__(self):
        """Add some basic print information"""
        return f"{self.iid}: Sex - {self.sex}"

    def __getitem__(self, item):
        """Get an item from Variant"""
        return getattr(self, item)

    @property
    def valid(self):
        """If we have valid information for this ID"""
        if (self.phenotype == -9) or (self.sex == 0):
            return False
        else:
            return True


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
