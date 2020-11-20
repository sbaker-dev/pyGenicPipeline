"""
This is a modified version of pandas plink found at https://github.com/limix/pandas-plink to act more like plinkio
available at https://github.com/mfranberg/libplinkio
"""
from pyGeneticPipe.geneticParsers.supportObjects import BimLoci, BimByChromosome
from pyGeneticPipe.geneticParsers.plink.bedObject import BedObject
from pyGeneticPipe.utils import error_codes as ec
from pathlib import Path
import numpy as np


class PlinkObject:
    def __init__(self, plink_path):
        self.base_path = plink_path
        self.bed_path, self.bim_path, self.fam_path = self.validate_paths()

    def bed_object(self):
        print("Starting bed")
        variant_number, sample_number = self.get_dimensions()
        BedObject(self.bed_path, variant_number, sample_number)

    def bim_object(self):
        """
        Turns the bim file into a set of objects to be called.

        :return: A BimObject for each line in the bim file
        """
        with open(self.bim_path) as f:
            return [BimLoci(line) for line in f]

    def bim_by_chromosome(self):
        """
        When validating, we may need to cross compare by chromosome. THis will take the bim information and reformat it
        into a by chromosome rather than by snp.
        """

        # Extract the loci from the bim file
        bim_loci = self.bim_object()

        # Extract the unique sorted chromosomes
        valid_chromosomes = np.unique([int(loci.chromosome) for loci in bim_loci])
        valid_chromosomes.sort()

        # Create a dict, where each key is chromosome with its snps ids, indexes of those ids, base pair positions of
        # those ids and the nucleotides of those ids
        chr_dict = {}
        for chromosome in valid_chromosomes:
            chr_dict[str(chromosome)] = {'sids': [], 'snp_indices': [], 'positions': [], 'nts': []}

        # For everything found within the loci via bim, append it to the sub dict within the master dict
        for i, loci in enumerate(bim_loci):
            chr_dict[loci.chromosome]['sids'].append(loci.variant_id)
            chr_dict[loci.chromosome]['snp_indices'].append(i)
            chr_dict[loci.chromosome]['positions'].append(loci.bp_position)
            chr_dict[loci.chromosome]['nts'].append([loci.a1, loci.a2])

        # Return a ChromosomeBimObj for each chromosome
        return [BimByChromosome(chromosome, values["sids"], values["snp_indices"], values["positions"], values["nts"])
                for chromosome, values in zip(chr_dict.keys(), chr_dict.values())]

    def validation_snps(self, accepted_chromosomes, hap_map_3):
        """
        Create a set of valid snps and a position map to them via plink or bgen with option censuring of chromosome and
        snps via HapMap3

        :param accepted_chromosomes: Allowed chromosomes, if a chromosome is found for a given snp that is not in this
         list then this snp will be skipped

        :param hap_map_3: If you only want to use hap_map_3 snps then you need to provide a path to that file

        :return: A dict of snps where each snp has a key of base pair position (Position) and Chromosome AND
                 A set of all the valid snps AND
                 A set of all the valid chromosomes
        """
        snp_pos_map = {}
        valid_snp = set()
        valid_chromosomes = set()

        with open(self.bim_path) as f:
            for line in f:
                bim = BimLoci(line)
                # If the user has specified certain chromosomes check that this snps chromosome is in accepted_list
                if accepted_chromosomes and (bim.chromosome not in accepted_chromosomes):
                    continue

                # If the user has specified only to use snp id's from HapMap3 then check this condition
                if hap_map_3 and (bim.variant_id in hap_map_3):
                    valid_snp.add(bim.variant_id)
                    snp_pos_map[bim.variant_id] = {"Position": int(bim.bp_position), "Chromosome": bim.chromosome}
                    valid_chromosomes.add(bim.chromosome)

                # Otherwise add to valid snps / snp_pos_map
                else:
                    valid_snp.add(bim.variant_id)
                    snp_pos_map[bim.variant_id] = {"Position": int(bim.bp_position), "Chromosome": bim.chromosome}
                    valid_chromosomes.add(bim.chromosome)

        return snp_pos_map, valid_snp, valid_chromosomes

    def validate_paths(self):
        """
        Users may submit a path to a specific file within plink, such as a .bed/.bim/.fam or they just provide the root
        name. This method validates and returns the paths.

        :return: The path to the bed, bim, and fam file in that order
        """
        # Construct path as an object
        ld_path = Path(self.base_path)

        # Check file home directory can be reached
        assert ld_path.parent.exists(), ec.path_invalid(ld_path.parent, "_set_ld_ref")

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

        return bed, bim, fam

    def get_dimensions(self):
        """
        When working with Bed files we need to know the number of variants, the length of the bim file, and the number
        of individuals, the length of the fam file.

        :return: The length of the variants and number of samples
        """
        # Get the number of variants in the bim file
        with open(self.bim_path, "r") as f:
            variant_length = len(list(f))
        f.close()

        # Get the number of samples in the fam file
        with open(self.fam_path, "r") as f:
            number_of_samples = len(list(f))
        f.close()
        return variant_length, number_of_samples
