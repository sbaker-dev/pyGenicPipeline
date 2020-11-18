from datetime import datetime
from itertools import chain
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
