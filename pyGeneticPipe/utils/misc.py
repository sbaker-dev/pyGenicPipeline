from datetime import datetime
from itertools import chain
from pathlib import Path
import numpy as np
import yaml
import gzip
import os


def bits_to_int(bits):
    """Converts bits to int."""
    result = 0
    for bit in bits:
        result = (result << 1) | bit

    return result


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


def line_array(line_index, line_array, type_np=None):
    """
    Construct an array of a single line index form line_array considering the type
    :rtype: np.ndarray
    """
    if type_np:
        return np.array([line[line_index] for line in line_array], dtype=type_np)
    else:
        return np.array([line[line_index] for line in line_array])


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
