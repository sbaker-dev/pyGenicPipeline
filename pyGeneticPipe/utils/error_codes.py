

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
def job_violation(jobs):
    return f"TO MANY JOBS SPECIFIED\n" \
           f"pyGeneticPipe is setup to take a single process at a time, controlled via shell scripts or via separate " \
           f"calls to the respective Classes within it.\nHowever in the current operation {len(jobs)} where found to " \
           f"be set to run: {jobs}"


def process_not_run(job, project, non_process):
    return f"JOB {job} REQUIRES {non_process} SUBPROCESS BUT IT HAS NOT BEEN PERFORMED FOR {project} "


def appending_error(project, keys):
    return f"JOB HAS ALREADY BEEN UNDERTAKEN FOR PROJECT {project}\n" \
           f"The h5py file already contains {keys}, if you want to re-do the project you need to set Override to " \
           f"True\nWARNING: This will delete ALL data within the h5py file, not just this part"


def missing_arg(job, arg):
    return f"{job} requires {arg} but it was not set!\n"


def path_invalid(parent_path, operation):
    return f"INVALID Path for {operation}\n" \
           f"{operation} attempt to navigate to a directory or file at {parent_path} but it does not exist"


def invalid_effect_type(effect_types, effect_found):
    return f"INVALID EFFECT TYPE of {effect_found}\n" \
           f"Summary statistics from GWAS's can have effect types of {effect_types} yet found {effect_found}"


def z_scores_with_standard_errors():
    return f"CALCULATING Z SCORES SET YET NO STANDARD ERRORS IN SUMMARY STATS\n" \
           f"A standard error column was not specified or found in the gwas summary stats so z scores cannot be " \
           f"calculated"


def no_valid_snps(bim_file, chromosomes_status, hm3_status):
    if chromosomes_status is not None:
        chromosomes_status = True
    if hm3_status is not None:
        hm3_status = True

    return f"NO VALID SNPS FOUND WITHIN {bim_file}\n" \
           f"You currently have set custom_chromosomes to : {chromosomes_status} and HapMap3 censuring" \
           f" to {hm3_status}\n, and with these options zero snps have been found in {bim_file}"


def sample_size():
    return f"NO SAMPLE SIZE FOUND\n" \
           f"When working on constructing summary statistics the sample size of the GWAS studies summary statistics" \
           f"is required."


def validation_size_invalid(size):
    return f"INVALID SAMPLE SIZE {size}\n" \
           f"Validation_Size must be between 0 and 1 inclusively"


def mandatory_header(header, found_headers, match_headers):
    return f"MANDATORY HEADER {header} NOT SET\n" \
           f"When reading in summary statistics {header} is required yet\n" \
           f"Tried to find {header} in headers: {found_headers} by matching {match_headers}"


def ambiguous_header(header, found_headers, match_headers):
    return f"AMBIGUOUS HEADERS FOUND FOR {header}\n" \
           f"When reading the headers linknig should be unique yet\n" \
           f"Founding {found_headers} indexes for {header} by matching {match_headers}"


# bed file codes
def bed_magic_violation(file_name, magic):
    return f"INVALID BED FILE for file at path: {file_name}\n" \
           f"BED files have a magic number of 6c1b in hex for the first two bytes but found {magic}"


def bed_matrix_order(file_name, order):
    return f"INVALID BED FILE ORDER for bed file at path: {file_name}\n" \
           f"Bed files have the third byte relate to the order, a 00 relating to sample and 01 variant yet found " \
           f"BED file matrix order {order}"
