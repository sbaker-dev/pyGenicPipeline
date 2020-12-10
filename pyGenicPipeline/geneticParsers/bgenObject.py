"""
This is a modified version of the pybgen project available at https://github.com/lemieuxl/pybgen
"""
from pyGenicPipeline.geneticParsers.variantObjects import Variant
from pyGenicPipeline.utils import error_codes as ec
from pyGenicPipeline.utils.misc import bits_to_int
from pathlib import Path
import numpy as np
import sqlite3
import struct
import zlib
import zstd
import os


class BgenObject:
    def __init__(self, file_path, bgi_present=True, probability=None, iter_array_size=1000,
                 sid_index=slice(None, None, None), iid_index=slice(None, None, None)):
        """

        :param file_path:

        :param bgi_present: Takes a value of True if the .bgi is in the same directory and named file_path.bgi
            otherwise can ec passed as a path if it is in a different directory.
        :type bgi_present: bool | str

        :param probability:
        """

        self.file_path = Path(file_path)
        self._bgen_binary = open(file_path, "rb")

        self.offset, self.headers, self.sid_count, self.iid_count, self.compression, self.layout, \
            self.sample_identifiers, self._variant_start = self.parse_header()

        self.probability = probability
        self.iter_array_size = iter_array_size

        self.bgi_present = bgi_present
        self.bgi_file = self._set_bgi()
        if self.bgi_file:
            self.bgen_file, self.bgen_index, self.last_variant_block = self._connect_index
        else:
            self.bgen_file, self.bgen_index, self.last_variant_block = None, None, None

        self.iid_index = iid_index
        self.sid_index = sid_index

    def __getitem__(self, item):
        if isinstance(item, slice):
            return BgenObject(self.file_path, self.bgi_present, self.probability, self.iter_array_size, sid_index=item)

        elif isinstance(item, tuple):
            assert np.sum([isinstance(s, slice) for s in item]) == 2, ec
            return BgenObject(self.file_path, self.bgi_present, self.probability, self.iter_array_size,
                              sid_index=item[0], iid_index=item[1])
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
        compression, layout, sample_identifiers = self._header_flag()
        return offset, headers, variant_number, sample_number, compression, layout, sample_identifiers, variant_start

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
        compression_flag = bits_to_int(flag[0: 2][::-1])
        assert 0 <= compression_flag < 3, ec.compression_violation(self._bgen_binary.name, compression_flag)
        if compression_flag == 0:
            compression = self._no_decompress
        elif compression_flag == 1:
            compression = zlib.decompress
        else:
            compression = zstd.decompress

        # Check the layout is either 1 or 2, see [N1]
        layout = bits_to_int(flag[2:6][::-1])
        assert 1 <= layout < 3, ec.layout_violation(self._bgen_binary.name, layout)

        # Check if the sample identifiers are in the file or not, then return
        assert flag[31] == 0 or flag[31] == 1, ec.sample_identifier_violation(self._bgen_binary.name, flag[31])
        if flag[31] == 0:
            return compression, layout, False
        else:
            return compression, layout, True

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
        assert nb_markers == self.sid_count, ec
        assert first_variant_block == self._variant_start

        # Checking the number of markers
        if nb_markers != self.sid_count:
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
            # will return variant and dosage
            raise NotImplementedError("Dosage Not yet implemented")
        else:
            return variant

    def _get_curr_variant_info(self):
        """Gets the current variant's information."""

        if self.layout == 1:
            assert self.unpack("<I", 4) == self.iid_count, ec

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
        return {name: seek for seek, name in self.bgen_index.fetchall()}

    def index_from_sid(self, snp_names):
        """
        Construct the seek index for all snps provide as a list or tuple of snp_names
        """
        assert self.bgen_index, ec.bgen_index_violation("get_variant")

        # Select all the variants where the rsid is in the names provided
        self.bgen_index.execute("SELECT file_start_position FROM Variant WHERE rsid IN {}".format(tuple(snp_names)))

        # Fetching all the seek positions
        seek_positions = [index[0] for index in self.bgen_index.fetchall()]

        # Return a dict of type {Name: seek}
        return {name: seek for name, seek in zip(snp_names, seek_positions)}

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
