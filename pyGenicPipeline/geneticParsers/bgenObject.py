"""
This is a modified version of the pybgen project available at https://github.com/lemieuxl/pybgen
"""
from pyGenicPipeline.geneticParsers.variantObjects import Variant
from pyGenicPipeline.utils import error_codes as ec
from pyGenicPipeline.utils import misc as mc
from pathlib import Path
import numpy as np
import sqlite3
import struct
import zlib
import zstd
import os


class BgenObject:
    def __init__(self, file_path, bgi_present=True, probability_return=None, iter_array_size=1000, prob=0.9,
                 iid_index=slice(None, None, None), sid_index=slice(None, None, None), sample_path=None):
        """
        :param file_path:

        :param bgi_present: Takes a value of True if the .bgi is in the same directory and named file_path.bgi
            otherwise can ec passed as a path if it is in a different directory.
        :type bgi_present: bool | str

        :param probability_return:
        """

        self.file_path = Path(file_path)
        self._bgen_binary = open(file_path, "rb")
        self.iid_index = iid_index
        self.sid_index = sid_index

        self.offset, self.headers, self._variant_number, self._sample_number, self.compression, self.compressed, \
            self.layout, self.sample_identifiers, self._variant_start = self.parse_header()

        # Index our sid and iid values if we have indexes
        self.sid_count = len(np.arange(self._variant_number)[self.sid_index])
        self.iid_count = len(np.arange(self._sample_number)[self.iid_index])
        self.iid = self._set_iid(sample_path)

        self.probability_return = probability_return
        self.probability = prob
        self.iter_array_size = iter_array_size

        self.bgi_present = bgi_present
        self.bgi_file = self._set_bgi()
        if self.bgi_file:
            self.bgen_file, self.bgen_index, self.last_variant_block = self._connect_index
        else:
            self.bgen_file, self.bgen_index, self.last_variant_block = None, None, None

    def __getitem__(self, item):
        """Return a new BgenObject with slicing set."""
        if isinstance(item, slice):
            return BgenObject(self.file_path, self.bgi_present, self.probability_return, self.iter_array_size,
                              self.probability, sid_index=item)

        elif isinstance(item, tuple):
            assert np.sum([isinstance(s, slice) for s in item]) == 2, ec
            return BgenObject(self.file_path, self.bgi_present, self.probability_return, self.iter_array_size,
                              self.probability, iid_index=item[0], sid_index=item[1])
        else:
            raise Exception("Sc")

    def parse_header(self):
        """
        Extract information from the header of the bgen file.

        Spec at https://www.well.ox.ac.uk/~gav/bgen_format/spec/latest.html

        :return: offset, headers, variant_number, sample_number, compression, layout, and sample_identifiers
        """

        # Check the header block is not larger than offset
        offset = self.unpack("<I", 4)
        headers = self.unpack("<I", 4)
        assert headers <= offset, ec.offset_violation(self._bgen_binary.name, offset, headers)
        variant_start = offset + 4

        # Extract the number of variants and samples
        variant_number = self.unpack("<I", 4)
        sample_number = self.unpack("<I", 4)

        # Check the file is valid
        magic = self.unpack("4s", 4)
        assert (magic == b'bgen') or (struct.unpack("<I", magic)[0] == 0), ec.magic_violation(self._bgen_binary.name)

        # Skip the free data area
        self._bgen_binary.read(headers - 20)

        # Extract the flag, then set compression layout and sample identifiers from it
        compression, compressed, layout, sample_identifiers = self._header_flag()
        return (offset, headers, variant_number, sample_number, compression, compressed, layout,
                sample_identifiers, variant_start)

    def _header_flag(self):
        """
        The flag represents a 4 byte unsigned int, where the bits relates to the compressedSNPBlock at bit 0-1, Layout
        at 2-5, and sampleIdentifiers at 31

        Spec at https://www.well.ox.ac.uk/~gav/bgen_format/spec/latest.html

        :return: Compression, layout, sampleIdentifiers
        """
        # Reading the flag
        flag = np.frombuffer(self._bgen_binary.read(4), dtype=np.uint8)
        flag = np.unpackbits(flag.reshape(1, flag.shape[0]), bitorder="little")

        # [N1] Bytes are stored right to left hence the reverse, see shorturl.at/cOU78
        # Check the compression of the data
        compression_flag = mc.bits_to_int(flag[0: 2][::-1])
        assert 0 <= compression_flag < 3, ec.compression_violation(self._bgen_binary.name, compression_flag)
        if compression_flag == 0:
            compressed = False
            compression = self._no_decompress
        elif compression_flag == 1:
            compressed = True
            compression = zlib.decompress
        else:
            compressed = True
            compression = zstd.decompress

        # Check the layout is either 1 or 2, see [N1]
        layout = mc.bits_to_int(flag[2:6][::-1])
        assert 1 <= layout < 3, ec.layout_violation(self._bgen_binary.name, layout)

        # Check if the sample identifiers are in the file or not, then return
        assert flag[31] == 0 or flag[31] == 1, ec.sample_identifier_violation(self._bgen_binary.name, flag[31])
        if flag[31] == 0:
            return compression, compressed, layout, False
        else:
            return compression, compressed, layout, True

    def unpack(self, struct_format, size, list_return=False):
        """
        Use a given struct formatting to unpack a byte code

        Struct formatting
        ------------------
        https://docs.python.org/3/library/struct.html

        :param struct_format: The string representation of the format to use in struct format. See struct formatting for
            a list of example codes.
        :type struct_format: str

        :param size: The byte size
        :type size: int

        :key list_return: If we expect multiple values then we return a tuple of what was unpacked, however if there is
            only one element then we often just index the first element to return it directly. Defaults to false.
        :type list_return: bool

        :return: Whatever was unpacked
        :rtype: Any
        """

        if list_return:
            return struct.unpack(struct_format, self._bgen_binary.read(size))
        else:
            return struct.unpack(struct_format, self._bgen_binary.read(size))[0]

    def _read_bgen(self, struct_format, size):
        """
        Sometimes we need to read the number of bytes read via unpack

        :param struct_format: The string representation of the format to use in struct format. See struct formatting for
            a list of example codes.
        :type struct_format: str

        :param size: The byte size
        :type size: int

        :return: Decoded bytes that where read
        """
        return self._bgen_binary.read(self.unpack(struct_format, size)).decode()

    @property
    def _connect_index(self):
        """Connect to the index (which is an SQLITE database)."""
        bgen_file = sqlite3.connect(str(self.file_path.absolute()) + ".bgi")
        bgen_index = bgen_file.cursor()

        # Fetching the number of variants and the first and last seek position
        bgen_index.execute(
            "SELECT COUNT (rsid), "
            "       MIN (file_start_position), "
            "       MAX (file_start_position) "
            "FROM Variant"
        )
        nb_markers, first_variant_block, last_variant_block = bgen_index.fetchone()

        # Check the number of markers are the same across bgen and bgi, and that they start in the same block
        assert nb_markers == self._variant_number, ec
        assert first_variant_block == self._variant_start

        # Checking the number of markers
        if nb_markers != self._variant_number:
            raise ValueError("Number of markers different between headers of bgen and bgi")

        # Checking the first variant seek position
        if first_variant_block != self._variant_start:
            raise ValueError(f"{self.file_path.name}: invalid index")

        return bgen_file, bgen_index, last_variant_block

    def iter_variant_info(self):
        """Iterate over marker information."""
        assert self.bgen_index, ec.bgen_index_violation("iter_variant_info")

        self.bgen_index.execute(
            "SELECT chromosome, position, rsid, allele1, allele2 FROM Variant",
        )

        # Fetching the results
        results = self.bgen_index.fetchmany(self.iter_array_size)
        while results:
            for chromosome, position, variant_id, a1, a2 in results:
                yield Variant(chromosome, position, variant_id, a1, a2)
            results = self.bgen_index.fetchmany(self.iter_array_size)

    def get_variant(self, seek, dosage=False):
        """
        Use the index of seek to move to the location of the variant in the file, then return the variant as Variant
        """
        self._bgen_binary.seek(seek)
        return self._read_current_variant(dosage)

    def _read_current_variant(self, dosage):
        """Reads the current variant."""
        # Getting the variant's information
        variant = self._get_curr_variant_info()

        if dosage:
            return self._get_curr_variant_data()
        else:
            return variant

    def _get_curr_variant_info(self):
        """Gets the current variant's information."""

        if self.layout == 1:
            assert self.unpack("<I", 4) == self._sample_number, ec

        # Reading the variant id (may be in form chr1:8045045:A:G or just a duplicate of rsid and not used currently)
        self._read_bgen("<H", 2)

        # Reading the variant rsid
        rs_id = self._read_bgen("<H", 2)

        # Reading the chromosome
        chromosome = self._read_bgen("<H", 2)

        # Reading the position
        pos = self.unpack("<I", 4)

        # Getting the alleles
        alleles = [self._read_bgen("<I", 4) for _ in range(self._set_number_of_alleles())]

        # Return the Variant - currently only supports first two alleles
        return Variant(chromosome, pos, rs_id, alleles[0], alleles[1])

    def sid_array(self):
        """
        Construct an array of all the snps that exist in this file
        """
        assert self.bgen_index, ec.bgen_index_violation("get_variant")

        # Fetching all the seek positions
        self.bgen_index.execute("SELECT rsid FROM Variant")

        return np.array([name for name in self.bgen_index.fetchall()[self.sid_index]]).flatten()

    def sid_indexer(self):
        """
        Construct the seek index for all the snps in the file
        """
        assert self.bgen_index, ec.bgen_index_violation("get_variant")

        # Fetching all the seek positions
        self.bgen_index.execute("SELECT file_start_position, rsid FROM Variant")

        # Return a dict of type {Name: seek}
        return {name: seek for seek, name in self.bgen_index.fetchall()[self.sid_index]}

    def dosage_from_sid(self, snp_names):
        """
        Construct the seek index for all snps provide as a list or tuple of snp_names
        """
        assert self.bgen_index, ec.bgen_index_violation("get_variant")

        # Select all the variants where the rsid is in the names provided
        self.bgen_index.execute("SELECT file_start_position FROM Variant WHERE rsid IN {}".format(tuple(snp_names)))

        return np.array([self.get_variant(seek[0], True)[self.iid_index] for seek in self.bgen_index.fetchall()])

    def _set_number_of_alleles(self):
        """
        Bgen version 2 can allow for more than 2 alleles, so if it is version 2 then unpack the number stored else
        return 2
        :return: number of alleles for this snp
        :rtype: int
        """
        if self.layout == 2:
            return self.unpack("<H", 2)
        else:
            return 2

    def _get_curr_variant_data(self):
        """Gets the current variant's dosage or probabilities."""

        if self.layout == 1:
            print("WARNING - UNTESTED CODE FROM PY-BGEN")
            # Getting the probabilities
            probs = self._get_curr_variant_probs_layout_1()

            if self.probability:
                # Returning the probabilities
                return probs

            else:
                # Returning the dosage
                return self._layout_1_probs_to_dosage(probs)

        else:
            # Getting the probabilities
            probs, missing_data = self._get_curr_variant_probs_layout_2()

            if self.probability_return:
                # Getting the alternative allele homozygous probabilities
                last_probs = self._get_layout_2_last_probs(probs)

                # Stacking the probabilities
                last_probs.shape = (last_probs.shape[0], 1)
                full_probs = np.hstack((probs, last_probs))

                # Setting the missing to NaN
                full_probs[missing_data] = np.nan

                # Returning the probabilities
                return full_probs

            else:
                # Computing the dosage
                dosage = self._layout_2_probs_to_dosage(probs)

                # Setting the missing to NaN
                dosage[missing_data] = np.nan

                # Returning the dosage
                return dosage

    def _get_curr_variant_probs_layout_1(self):
        """Gets the current variant's probabilities (layout 1)."""
        c = self._sample_number
        if self.compressed:
            c = self.unpack("<I", 4)

        # Getting the probabilities
        probs = np.frombuffer(
            self.compression(self._bgen_binary.read(c)),
            dtype="u2",
        ) / 32768
        probs.shape = (self._sample_number, 3)

        return probs

    def _layout_1_probs_to_dosage(self, probs):
        """Transforms probability values to dosage (from layout 1)"""
        # Constructing the dosage
        dosage = 2 * probs[:, 2] + probs[:, 1]
        if self.probability > 0:
            dosage[~np.any(probs >= self.probability, axis=1)] = np.nan

        return dosage

    def _get_curr_variant_probs_layout_2(self):
        """Gets the current variant's probabilities (layout 2)."""
        # The total length C of the rest of the data for this variant
        c = self.unpack("<I", 4)

        # The number of bytes to read
        to_read = c

        # D = C if no compression
        d = c
        if self.compressed:
            # The total length D of the probability data after
            # decompression
            d = self.unpack("<I", 4)
            to_read = c - 4

        # Reading the data and checking
        data = self.compression(self._bgen_binary.read(to_read))
        assert len(data) == d, "INVALID HERE"

        # Checking the number of samples
        n = mc.struct_unpack("<I", data[:4])
        assert n == self._sample_number, "INVALID HERE"

        data = data[4:]

        # Checking the number of alleles (we only accept 2 alleles)
        nb_alleles = mc.struct_unpack("<H", data[:2])
        assert nb_alleles == 2, "INVALID HERE"
        data = data[2:]

        # TODO: Check ploidy for sexual chromosomes
        # The minimum and maximum for ploidy (we only accept ploidy of 2)
        min_ploidy = mc.byte_to_int(data[0])
        max_ploidy = mc.byte_to_int(data[1])
        if min_ploidy != 2 and max_ploidy != 2:
            raise ValueError("INVALID HERE")

        data = data[2:]

        # Check the list of N bytes for missingness (since we assume only
        # diploid values for each sample)
        ploidy_info = np.frombuffer(data[:n], dtype=np.uint8)
        ploidy_info = np.unpackbits(
            ploidy_info.reshape(1, ploidy_info.shape[0]).T,
            axis=1,
        )
        missing_data = ploidy_info[:, 0] == 1
        data = data[n:]

        # TODO: Permit phased data
        # Is the data phased?
        is_phased = data[0] == 1
        if is_phased:
            raise ValueError(
                "{}: only accepting unphased data".format("INVALID")
            )
        data = data[1:]

        # The number of bits used to encode each probabilities
        b = mc.byte_to_int(data[0])
        data = data[1:]

        # Reading the probabilities (don't forget we allow only for diploid
        # values)
        if b == 8:
            probs = np.frombuffer(data, dtype=np.uint8)

        elif b == 16:
            probs = np.frombuffer(data, dtype=np.uint16)

        elif b == 32:
            probs = np.frombuffer(data, dtype=np.uint32)

        else:
            probs = mc.pack_bits(data, b)

        # Changing shape and computing dosage
        probs.shape = (self._sample_number, 2)

        return probs / (2 ** b - 1), missing_data

    @staticmethod
    def _get_layout_2_last_probs(probs):
        """
        Gets the layout 2 last probabilities (homo alternative).
        :rtype: np.ndarray
        """
        return 1 - np.sum(probs, axis=1)

    def _layout_2_probs_to_dosage(self, probs):
        """Transforms probability values to dosage (from layout 2)."""
        # Computing the last genotype's probabilities
        last_probs = self._get_layout_2_last_probs(probs)

        # Constructing the dosage
        dosage = 2 * last_probs + probs[:, 1]

        # Setting low quality to NaN
        if self.probability > 0:
            good_probs = (
                    np.any(probs >= self.probability, axis=1) |
                    (last_probs >= self.probability)
            )
            dosage[~good_probs] = np.nan

        return dosage

    def _set_bgi(self):
        """
        Connect to the index file either via a bgi file in the same directory or in another directory.

        """

        if not self.bgi_present:
            return False
        elif self.bgi_present:
            if not os.path.isfile(f"{self.file_path}.bgi"):
                raise IOError(f"{self.file_path}.bgi was not found")
            else:
                return True
        elif isinstance(self.bgi_present, str):
            if not os.path.isfile(f"{self.bgi_present}"):
                raise IOError(f"{self.bgi_present} was not found")
            else:
                return True
        else:
            raise TypeError(ec.bgi_path_violation(self.bgi_present))

    @staticmethod
    def _no_decompress(data):
        """Don't decompress"""
        return data

    def _set_iid(self, sample_path):
        """
        If sample identifiers are within the bgen file then these can be extracted and set, however this has not yet
        been tested and am unsure if sex and missing stored in bgen as is with .sample files?

        If a path to the samples has been provided, then we can load the information within it. Sample files contain
        both the FID and the IID as well as missing and sex allowing us more options akin to .fam files.

        If Nothing is provided, and nothing is embedded, we create a list of id's on then number of id's after indexing.

        :param sample_path: Path to sample file
        :type sample_path: str | None

        :return: An array of id information
        """
        if self.sample_identifiers:
            raise NotImplementedError("Sorry, i haven't found a bgen file with id's within it yet to test")
        elif sample_path:
            raise NotImplementedError("Sorry this needs to be tested")
        else:
            return np.array([[i, i] for i in np.arange(self._sample_number)[self.iid_index]])
