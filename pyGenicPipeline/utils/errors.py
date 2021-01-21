def job_violation(jobs):
    return f"TO MANY JOBS SPECIFIED\n" \
           f"pyGeneticPipe is setup to take a single process at a time, controlled via shell scripts or via separate " \
           f"calls to the respective Classes within it.\nHowever in the current operation {len(jobs)} where found to " \
           f"be set to run: {jobs}"


def job_type(jt):
    return f"INVALID JOB TYPE\n" \
           f"May be None if the user is using the main method as an object\n" \
           f"May be a string if the user is using a job submission form\n" \
           f"May be a dict if the user is using the method for development or natively in python\n" \
           f"Yet Found {jt}"


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


def snp_overflow(summary_length, validation_length):
    return f"CRITICAL ERROR: MORE SUMMARY SNPS SELECTED THAN EXIST IN VALIDATION\n" \
           f"Cleaning summary statistics checks if a snp is within the genetic file or not, and only then cleans it" \
           f".\nTherefore, it should be impossible to have more snps select than the number found in the validation." \
           f"\nYet Found Summary: {len(summary_length)} Validation: {validation_length}"


def all_missing(attributes_str, operation):
    return f"ERROR NO VALUES FOUND FOR {operation}\n" \
           f"{attributes_str} = 0"


def ld_radius_to_large():
    return "Warning: LD radius seems small in comparision to average LD score!\nConsider using a smaller radius or" \
           " less snps"


def herit_type(type_herit):
    return f"ERROR: PREDEFINED HERITABILITY IS NEITHER DICT NOR FLOAT\n" \
           f"If you are providing a heritability on a per chromosome baises then it needs to be a dict of chromosme: " \
           f"heritability\nIf you want to distribute heritability over the chromosomes then provide float\nIf you " \
           f"want to calculate the heritability then provide None, if float is set it will still calculate it\n" \
           f"Yet Found {type_herit}"


def heritability_to_large():
    return "Warning: Estimated Genome-wide Heritability is suspiciously large suggesting the sampple may be incorrect" \
           " or snps are enriched for heritability.\nYou may assign a dict of type chromosome: heritability to " \
           "heritability_calculated if you wish to override computed heritability"


def gibbs_convergence(variant_fraction, snp_count, sum_sq_beta, genome_herit):
    print(f"Warning: Sum Squared beta is much large than estimated heritability suggesting a lack of "
          f"convergence of Gibbs.\nCasual Variants used at Variant Fraction {variant_fraction}: "
          f"{variant_fraction * snp_count}\nSum Squared Beta > Genome Heritability: {sum_sq_beta} > {genome_herit}")


def scores_valid_headers(isolates, headers, score_files):
    assert len(isolates) >= 1, "CRITICAL ERROR - INFINITESIMAL FAILED TO LOAD"
    print(f"The following where dropped "
          f"{[f'{h}: {headers[h]}/{len(score_files)}' for h in headers.keys() if h not in isolates]}")

