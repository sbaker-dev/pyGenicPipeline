"""
This is a modified version of pandas plink: https://github.com/limix/pandas-plink
"""

from pyGeneticPipe.utils import error_codes as ec
import struct


class BedObject:
    def __init__(self, bed_path, variant_number, sample_number):
        # Required to load
        self._bed_path = bed_path
        self.variant_number = variant_number
        self.sample_number = sample_number

        # Open the file
        self._bed_binary = open(self._bed_path, "rb")
        self.sample_order = self._validate()

    def close(self):
        """
        When we are done with the file, close it
        :return: Nothing, just close the opened file
        """
        self._bed_binary.close()

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
            return struct.unpack(struct_format, self._bed_binary.read(size))
        else:
            return struct.unpack(struct_format, self._bed_binary.read(size))[0]

    def _validate(self):
        """
        Check the magic number as well as the order of the bed file before continuing

        :return: The order of the bed file, if it is variant order is sample, then return True else False where false
            will referer to variant ordering
        :rtype: bool
        """

        # Check the magic numbers
        magic = self.unpack("2s", 2).hex()
        assert magic == "6c1b", ec.bed_magic_violation(self._bed_path, magic)

        # Validate the order of the matrix where "00" is sample layout and "01" variant
        order = self.unpack("1s", 1).hex()
        assert order == "00" or order == "01", ec.bed_matrix_order(self._bed_path, order)
        if order == "00":
            return True
        else:
            return False
