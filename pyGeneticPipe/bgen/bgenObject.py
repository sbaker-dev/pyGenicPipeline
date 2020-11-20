"""
This is a modified version of the pybgen project available at https://github.com/lemieuxl/pybgen
"""

from pyGeneticPipe.utils import error_codes as ec
from pyGeneticPipe.utils.misc import bits_to_int
import numpy as np
import struct
import zlib
import zstd
import os


class BgenObject:
    def __init__(self, file_path, bgi_file_path=False, probability=None):
        """

        :param file_path:

        :param bgi_file_path: Takes a value of True if the .bgi is in the same directory and named file_path.bgi
            otherwise can ec passed as a path if it is in a different directory.
        :type bgi_file_path: bool | str

        :param probability:
        """

        self.file_path = Path(file_path)
        self._bgen_binary = open(file_path, "rb")

        self.offset, self.headers, self.variant_number, self.sample_number, self.compression, self.layout, \
            self.sample_identifiers = self.parse_header()

        self.probability = probability
        self.bgi_file = self._set_bgi(file_path, bgi_file_path)

        print(self.offset, self.headers, self.variant_number, self.sample_number, self.compression, self.layout,
              self.sample_identifiers)

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
        return offset, headers, variant_number, sample_number, compression, layout, sample_identifiers

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

    @staticmethod
    def _set_bgi(bgen_path, bgi_file_path):
        """
        Connect to the index file either via a bgi file in the same directory or in another directory.

        :param bgen_path: The path to bgen_path
        :param bgi_file_path: The bgi_position
        """

        if not bgi_file_path:
            return False
        elif bgi_file_path:
            if not os.path.isfile(f"{bgen_path}.bgi"):
                raise IOError(f"{bgen_path}.bgi was not found")
            else:
                pass  # This is where we set the .bgi file
        elif isinstance(bgi_file_path, str):
            if not os.path.isfile(f"{bgi_file_path}"):
                raise IOError(f"{bgi_file_path} was not found")
            else:
                pass  # This is where we set the .bgi file
        else:
            raise TypeError(ec.bgi_path_violation(bgi_file_path))

    @staticmethod
    def _no_decompress(data):
        return data
