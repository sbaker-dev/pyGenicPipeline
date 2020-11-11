class BimObject:
    def __init__(self, bim_line):
        """
        Holds information of a bim file in Object notation. Bim files contain a set 6 categories. Full specification:

        https://www.cog-genomics.org/plink2/formats#bim

        """
        self.chromosome, self.variant_id, self.morgan_pos, self.bp_cord, self.a1, self.a2 = bim_line.split()

