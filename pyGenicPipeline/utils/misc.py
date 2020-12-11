from itertools import combinations
from datetime import datetime
from itertools import chain
from colorama import Fore
from pathlib import Path
from math import ceil
import numpy as np
import struct
import time
import yaml
import gzip
import os


def bits_to_int(bits):
    """Converts bits to int."""
    result = 0
    for bit in bits:
        result = (result << 1) | bit

    return result


def byte_to_int(byte):
    """Converts a byte to a int for python 3."""
    return byte


def pack_bits(data, b):
    """Unpacks BGEN probabilities (as bits)."""
    # Getting the data from the bytes
    data = np.fromiter(
        ((byte_to_int(byte) >> i) & 1 for byte in data for i in range(8)),
        dtype=bool,
    )
    data.shape = (data.shape[0] // b, b)

    # Finding the closest full bytes (if required)
    # TODO: Improve this so that it is more efficient
    full_bytes = data[:, ::-1]
    if data.shape[1] % 8 != 0:
        nb_bits = int(ceil(b / 8)) * 8
        full_bytes = np.zeros((data.shape[0], nb_bits), dtype=bool)
        full_bytes[:, -b:] += data[:, ::-1]

    # Packing the bits
    packed = np.packbits(full_bytes, axis=1)

    # Left-shifting for final value
    final = packed[:, 0]
    for i in range(1, packed.shape[1]):
        final = np.left_shift(final, 8, dtype=np.uint) | packed[:, i]

    return final


def struct_unpack(struct_format, data, list_return=False):
    if list_return:
        return struct.unpack(struct_format, data)
    else:
        return struct.unpack(struct_format, data)[0]


def open_setter(path):
    """
    Some files may be zipped, opens files according to the zip status

    :param path: File path
    :type path: Path
    :return: gzip.open if the file is gzipped else open
    """
    if path.suffix == ".gz":
        return gzip.open
    else:
        return open


def decode_line(line, zip_status):
    """
    Some files may be zipped, when we open zipped files we will need to decode them

    :param line: Current line from open file, zipped or otherwise
    :param zip_status: If the file is zipped or not
    :return: decoded line from the open file
    """
    if zip_status:
        return line.decode("utf-8").split()
    else:
        return line.readline().split()


def terminal_time():
    """
    A way to remember when you initialised a cell by return the current hour and minute as a string
    """
    return f"{datetime.now().time().hour}:{datetime.now().time().minute}:{datetime.now().time().second}"


def directory_iterator(directory, file_only=True):
    """
    This takes a directory and returns a list of entries within that directory, if file_only is selected only
    files as apposed to directories will be returned

    :param directory: The directory you wish to iterate through
    :type directory: str

    :param file_only: Defaults to true where the return is just a list of files
    :type file_only: bool

    :return: List of entries from the directory
    :rtype: list
    """

    if file_only:
        return [file for file in os.listdir(directory) if os.path.isfile(f"{directory}/{file}")]
    else:
        return [file for file in os.listdir(directory)]


def flatten(list_of_lists):
    """
    Flatten a list of lists into a list
    """
    return list(chain(*list_of_lists))


def load_yaml(path_to_file):
    """
    Load the yaml file from package into scope
    """
    with open(path_to_file, "r") as f:
        try:
            return yaml.safe_load(f)
        except yaml.YAMLError:
            raise yaml.YAMLError


def line_array(line_index, numpy_line_array, type_np=None):
    """
    Construct an array of a single line index form line_array considering the type
    :rtype: np.ndarray
    """
    if type_np:
        return np.array([line[line_index] for line in numpy_line_array], dtype=type_np)
    else:
        return np.array([line[line_index] for line in numpy_line_array])


def variant_array(variant_key, variant_numpy_array):
    """
    Construct an array of a single variant_array's item using getitem via variant key
    :rtype: np.ndarray
    """
    return np.array([variant[variant_key] for variant in variant_numpy_array])


def filter_array(dict_to_filter, array_filter):
    """
    Filter out anything that is no longer required If the length of the array becomes zero, pass an error code of
    1. Otherwise return 0.
    """
    for key, value in zip(dict_to_filter.keys(), dict_to_filter.values()):
        dict_to_filter[key] = value[array_filter]

    if np.array([len(value) for value in dict_to_filter.values()])[0] == 0:
        return None
    else:
        return "OK"


def shrink_r2_matrix(distance_dp, n):
    """
    The values of R squared to be no less or greater than -1 and 1 respectively so the value calculated from the dot
    product is bounded to -1 and 1. Any value that is less than 1 divided by the number of individuals minus 1
    (to account for 0 indexing) is set to zero to reduce the complexity of the finalised r squared matrix

    :param distance_dp: The dot product tp calculate the clipped r2 matrix off
    :param n: The number of individuals in the sample
    :return: Clipped and reduced r2 matrix
    """

    clipped_r2 = np.clip(distance_dp, -1, 1)
    clipped_r2[np.absolute(clipped_r2) < (1.0 / (n - 1))] = 0
    return clipped_r2


def snps_in_window(snps, window_start, number_of_snps, window_size):
    """
    When accessing a window of snps we don't look at the snps +/- radius from the current snp, but just iterate in
    chunks of r*2 through the snps as a window. However, the last iteration may be out of range so we take the minimum
    of the number of snps as a precaution to prevent that.

    :param snps: Normalised snps
    :param window_start: Start index from the window iteration
    :param number_of_snps: total number of snps
    :param window_size: The size, radius * 2, of the window
    :return:
    """

    return snps[window_start: min(number_of_snps, (window_start + window_size))]


def posterior_mean(cd, b2, n):
    """Calculate the posterior mean from the denominator and numerator"""
    d_const_b2_exp = cd['d_const'] * np.exp(-b2 * n / 2.0)
    numerator = cd['c_const'] * np.exp(-b2 / (2.0 * cd['hdmpn']))
    if not isinstance(d_const_b2_exp, complex):
        if not isinstance(numerator, complex) and (numerator != 0.0):
            postp = numerator / (numerator + d_const_b2_exp)
            assert type(postp) != complex, "Post mean not a real number"
            return postp
        else:
            return 0.0
    else:
        return 1.0


def cleanup_dict(dict_to_clean, keys_to_remove):
    """We will sometimes no longer need certain attributes, this will clean up the dict to save memory"""
    for clean in keys_to_remove:
        dict_to_clean.pop(clean, None)

    print("Remaining Keys")
    for key in dict_to_clean:
        print(f"{key} - {len(dict_to_clean[key])}")
    print("\n")


def reshape_dict_array(dict_to_shape, dict_key, dimensions=1):
    """Reshape an array stored as dict_key within a dict_to_shape to a given dimension"""
    dict_to_shape[dict_key] = dict_to_shape[dict_key].reshape(len(dict_to_shape[dict_key]), dimensions)


def possible_combinations(com_list):
    """Of a possible list """
    return sum([list(map(list, combinations(com_list, i))) for i in range(len(com_list) + 1)], [])[1:]


def error_dict_to_terminal(error_dict):
    """Print the error dict for this chromosome then reset the initialised to default 0"""
    for index, (k, v) in enumerate(zip(error_dict.keys(), error_dict.values())):
        if index == 0:
            print(Fore.LIGHTCYAN_EX + "{:<30} {}".format(k, v))
            print(Fore.LIGHTCYAN_EX + "----------------------------------------")
        else:
            print("{:<30} {}".format(k, v))
    print("\n")

    # Reset values to 0
    for k, v in zip(error_dict.keys(), error_dict.values()):
        if isinstance(v, (int, np.int32, np.int8)):
            error_dict[k] = 0
    return time.time()
