from pyGeneticPipe.plink.supportObjects import Nucleotide
from pyGeneticPipe.plink.plinkObject import PlinkObject
from pyGeneticPipe.utils.Input import Input
import numpy as np
import h5py


class ProcessSummary(Input):
    def __init__(self, args):
        super().__init__(args)

        self._bim_by_chromosomes = self._set_bim_by_chromosomes()
        self.error_dict = {"Ambiguous_SNP": [], "Non_Allowed_Allele": [], "Non_Matching": []}


    def pre_process_plink(self, hdf5_path):

        # Open the cleaned summary statistics
        summary_file = h5py.File(hdf5_path, "r")["Sum_Stats"]

        for chrom_bim in self._bim_by_chromosomes:

            # Extract the summary statistics data for this chromosome, if we fail to extract it continue
            chrom_summary = self._summary_data_by_chromosome(summary_file, chrom_bim)
            if not chrom_summary:
                continue

            # Validate the snp ids
            summary_sids, common_snps, ss_sid_dict = self._validate_snp_ids(chrom_summary, chrom_bim)

            # Ordering by positioning
            g_snp_map, ss_snp_map = self.order_by_position(chrom_bim, common_snps, ss_sid_dict)

            # Coordinate the snps
            ok_indices, ok_nts = self.coordinate(chrom_bim, chrom_summary, g_snp_map, ss_snp_map, summary_sids)


            break

    def order_by_position(self, chrom_bim, common_snps, ss_sid_dict):
        g_snp_map = [chrom_bim.indexed_snps[sid] for sid in common_snps]
        g_positions = np.array(chrom_bim.positions)[g_snp_map]
        order = np.argsort(g_positions)
        g_snp_map = (np.array(g_snp_map)[order]).tolist()
        common_sids = np.array(common_snps)[order]
        ss_snp_map = [ss_sid_dict[sid] for sid in common_sids]

        # Don't know what this for, or if its more part of validation ##############################################
        # g_ss_nt_concord_count = np.sum(
        #     g_nts[g_snp_map] == ss_nts[ss_snp_map]) / 2.0
        ############################################################################################################
        return g_snp_map, ss_snp_map

    def coordinate(self, chrom_bim, chrom_summary, g_snp_map, ss_snp_map, summary_sids):

        # load the data from the chrom_bim or chrom_summary accordingly
        g_nts = np.array(chrom_bim.nucleotides)
        ss_nts = (chrom_summary['nucleotide'][...]).astype("|S1")
        betas = chrom_summary['beta'][...]
        log_odds = chrom_summary['log_odds'][...]
        # todo freq

        # Identifying which SNPs have nucleotides that are ok..
        ok_nts = []
        ok_indices = {'g': [], 'ss': []}
        for g_i, ss_i in zip(g_snp_map, ss_snp_map):

            assert chrom_bim.snps[g_i] == summary_sids[ss_i], 'Some issues with coordinating the genotypes.'
            ss_nt = Nucleotide(ss_nts[ss_i].astype("U13"))  # Summary stats have been saved in byte code

            # Is the nucleotide ambiguous?
            # We have Major and Minor alleles but if its AT or CG then its the base pair i guess
            g_nt = Nucleotide(g_nts[g_i])
            if g_nt.to_tuple() in self.ambiguous_snps:
                self.error_dict["Ambiguous_SNP"].append(g_nt.to_list())
                continue

            # First check if nucleotide is sane (Ie it is only a t c and g)?
            if (g_nt.a1 not in self.allowed_alleles) or (g_nt.a2 not in self.allowed_alleles):
                self.error_dict["Non_Allowed_Allele"].append(g_nt.to_list())
                continue

            flip_status = self.flip_nucleotide(betas, g_nt, log_odds, ss_i, ss_nt)
            if not flip_status:
                continue

            ok_indices['g'].append(g_i)
            ok_indices['ss'].append(ss_i)
            ok_nts.append(g_nt)
        return ok_indices, ok_nts

    def flip_nucleotide(self, betas, g_nt, log_odds, ss_i, ss_nt):
        # Opposite strands for flip checking if required
        os_g_nt = Nucleotide([self.allele_flip[g_nt.a1], self.allele_flip[g_nt.a2]])
        if not (np.all(g_nt.to_list() == ss_nt.to_list()) or np.all(os_g_nt.to_list() == ss_nt.to_list())):

            flip_nts = (g_nt.a2 == ss_nt.a1 and g_nt.a1 == ss_nt.a2) or (
                    os_g_nt.a2 == ss_nt.a1 and os_g_nt.a1 == ss_nt.a2)
            if flip_nts:
                betas[ss_i] = -betas[ss_i]
                log_odds[ss_i] = -log_odds[ss_i]
                # todo freq
                return True
            else:
                self.error_dict["Non_Matching"].append([g_nt.to_list(), ss_nt.to_list(), os_g_nt.to_list()])
                print("Failed to flip")
                return None

    def _validate_snp_ids(self, chrom_summary, chrom_bim):
        """
        This will extract a common list of snps when compared to a validation file if it is set, else will just return
        the bim files snps

        :param chrom_summary:
        :param chrom_bim:
        :return:
        """
        # Extract the sids from summary
        summary_sids = (chrom_summary['snp_id'][...]).astype('<U30')
        if self.validation_file:
            common_snps = None
            raise NotImplementedError("Not yet implement")
        else:
            common_snps = summary_sids

        # A map from sid to index for summary stats
        ss_sid_dict = {sid: i for i, sid in enumerate(summary_sids)}

        # The indices to retain for the LD reference genotypes
        common_snps = np.intersect1d(common_snps, chrom_bim.snps)
        return summary_sids, common_snps, ss_sid_dict

    def pre_process_bgen(self):
        raise NotImplementedError("Bgen files not yet implemented")

    @staticmethod
    def _summary_data_by_chromosome(summary_file, chrom_bim):
        """
        Try to extract the data by chromosome, if we failed to extract it summary statistics file return None so it will
        not crash and warn the user this has happened.
        """
        try:
            return summary_file[chrom_bim.chromosome]
        except KeyError:
            print(f"WARNING: Failed to find {chrom_bim.chromosome} in summary statistics file")
            return None

    def _set_bim_by_chromosomes(self):
        """
        Extract the information in the bim file but by chromosomes
        """
        # todo: will need to be generalised for validation
        # Create a plink Objected
        plink_obj = PlinkObject(self.args["LD_Reference_Genotype"])

        # Extract the information by chromosomes
        return plink_obj.bim_by_chromosome()


