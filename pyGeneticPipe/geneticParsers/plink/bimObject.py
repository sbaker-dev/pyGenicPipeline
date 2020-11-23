from pyGeneticPipe.geneticParsers.supportObjects import BimVariant


class BimObject:
    def __init__(self, bim_path):
        self._bim_file = open(bim_path, "r")

    def construct_index(self):
        """
        Bim files need to be index via seek, so we can extract a given snp loci wihout having to store all of this of
        then in memory
        """
        cumulative_seek = 0
        indexer = {}
        for _ in self._bim_file:
            line = self._bim_file.readline()
            if len(line) > 0:
                _, rs_id, _, _, _, _ = line.split()
                indexer[rs_id] = cumulative_seek
                cumulative_seek += len(line)
        self._bim_file.close()
        return indexer

    def extract_loci(self, seek, as_variant=False):
        """
        Extract a given variants loci based on the seek index from construct_index
        :param seek: The seek index
        :param as_variant: If you want it as a standardised across parameter variant, or a Bim Variant with morgan pos
        :return: The line
        """
        self._bim_file.seek(seek)
        chromosome, variant_id, morgan_pos, bp_position, a1, a2 = self._bim_file.readline().split()

        if as_variant:
            return BimVariant(chromosome, variant_id, morgan_pos, bp_position, a1, a2).to_variant()
        else:
            return BimVariant(chromosome, variant_id, morgan_pos, bp_position, a1, a2)
