from ..utils import errors as ec

from miscSupports import load_yaml
from itertools import combinations
from colorama import Fore
from pathlib import Path
import numpy as np
import time
import gzip


def set_args(args):
    """
    Args may be set as a dict, or as a yaml file with its path passed as the args. If the later then we need to load
    in the dict of values
    """
    if isinstance(args, dict):
        return args
    else:
        yaml_path = Path(args)
        assert (yaml_path.exists() and yaml_path.suffix == ".yaml"), ec.path_invalid(yaml_path, "_set_args")

        return load_yaml(yaml_path)


def set_current_job(operation_dict):
    """
    Set the current job to be processed

    :param operation_dict:
        May be None if the user is using the main method as an object
        May be a string if the user is using a job submission form
        May be a dict if the user is using the method for development or natively in python
    :type operation_dict: None | str | dict

    :return:
        If None, Returns None
        If str, Returns operation_dict
        If Dict, Asserts that only 1 job is selected and then returns the job name as a string

    :raises TypeError, AssertionError:
        TypeError if job is not a None, str, or dict.
        AssertionError if the job dict contains more than a single job
    """
    if not operation_dict:
        return None
    elif isinstance(operation_dict, str):
        return operation_dict
    elif isinstance(operation_dict, dict):
        job = [job_name for job_name, run in zip(operation_dict.keys(), operation_dict.values()) if run]
        assert len(job) == 1, ec.job_violation(job)
        return job[0]
    else:
        raise TypeError(ec.job_type(type(operation_dict)))


def validate_path(path, allow_none=True):
    """
    We have multiple types of files and directories, some may be allow to be None as they will not be required
    whilst others like the working directory will always be required. This method is a generalisation of individual
    setters.

    :param path: Path to a directory or file
    :type path: str

    :param allow_none: Defaults to True, if true if a path is set to none it will just return None. If False, an
        assertion will be run to validate that it is not none. In both cases, should the file not be None, then the
        path is validated via Path.exists()
    :type allow_none: Bool

    :return: Path to the current file or directory if None return is not allowed, otherwise the Path return is
        optional and the return may be none.
    """
    if allow_none and not path:
        return None
    else:
        assert path and Path(path).exists(), ec.path_invalid(path, "_validate_path")
        return Path(path)


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
        return line.split()


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


def filter_array(dict_to_filter, array_filter, error_process):
    """
    Filter out anything that is no longer required If the length of the array becomes zero, pass an error code of
    1. Otherwise return 0.
    """
    for key, value in zip(dict_to_filter.keys(), dict_to_filter.values()):
        dict_to_filter[key] = value[array_filter]

    if np.array([len(value) for value in dict_to_filter.values()])[0] == 0:
        raise TypeError(f"Filter out all elements after filtering on {error_process}")


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


def reshape_dict_array(dict_to_shape, dict_key, dimensions=1):
    """Reshape an array stored as dict_key within a dict_to_shape to a given dimension"""
    dict_to_shape[dict_key] = dict_to_shape[dict_key].reshape(len(dict_to_shape[dict_key]), dimensions)


def possible_combinations(com_list):
    """Of a possible list """
    return sum([list(map(list, combinations(com_list, i))) for i in range(len(com_list) + 1)], [])[1:]


def error_dict_to_terminal(error_dict, process, start):
    """Print the error dict for this chromosome then reset the initialised to default 0"""
    print("")
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

    print(f"Finished {process} in {round(time.time() - start, 2)} Seconds\n")
