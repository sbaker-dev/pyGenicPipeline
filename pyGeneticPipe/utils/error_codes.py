

# Bgen error codes
def magic_violation(file_name):
    return f"INVALID BGEN FILE for file at path: {file_name}\n" \
           f"Bgen files have a magic number, which in more modern files is 4 bytes of b'bgen' but can" \
           f" also take four zero bytes in older files.\n"


def offset_violation(file_name, offset, header_block_length):
    return f"INVALID HEADER BLOCK for file at path: {file_name}\n" \
           f"Bgen files have a header block, that must not be larger than the offset. Yet found\n" \
           f"Offset: {offset} header_block_length {header_block_length}"


def compression_violation(file_name, compression_flag):
    return f"INVALID COMPRESSION FLAG for file at path: {file_name}\n" \
           f"Bgen files have a flag, where the first two bits represent the compression of the data with 0 being " \
           f"uncompressed, 1 being compressed via zlib, and 2 being compressed via z-standard. Yet found\n" \
           f"Compression flag: {compression_flag}"


def layout_violation(file_name, layout_flag):
    return f"INVALID LAYOUT FLAG for file at path: {file_name}\n" \
           f"Bgen files have a flag, where bits 2-5 represents a Layout Flags should only take a value of 1 or 2 and " \
           f"relates to how and what genetic information is stored. Yet Found\n" \
           f"Layout flag: {layout_flag}"


def sample_identifier_violation(file_name, sample_identifier):
    return f"INVALID SAMPLE IDENTIFIER FLAG for file at path: {file_name}\n" \
           f"Bgen files have a flag, where bit 31 represents if sample identifiers as within the file at 1 or not" \
           f" at 0. Yet found\n" \
           f"Sample Identifier flag: {sample_identifier}"


def bgi_path_violation(bgi_path):
    return f"INVALID BGI TYPE for bgi_path of type: {bgi_path}\n" \
           f"Bgi_path defaults to False where you don't have an index in an external file. If you do have .bgi file" \
           f" in the same directory, name the same thing but with .bgi on the end, then it should be True. If the " \
           f"file is in another directory you may also pass the path as a str in to bgi_path. Yet Found\n" \
           f"BGI path type: {type(bgi_path)}"


# Input error codes
def path_invalid(parent_path, operation):
    return f"INVALID Path for {operation}\n" \
           f"{operation} attempt to navigate to a directory or file at {parent_path} but it does not exist"
