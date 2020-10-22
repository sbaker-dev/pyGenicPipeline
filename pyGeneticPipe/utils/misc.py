import gzip


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
