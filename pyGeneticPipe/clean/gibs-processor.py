pval_derived_betas = "g[beta]"
h2 = None
n = 253288
ns = "g[ns]"  # An array of n.....
p = [1, 0.3, "etc"]  # The fraction of causal variants used in the gibs sampler
ld_radius = radius = 183
verbose = True
num_inter = 60  # Number of iterations to run the gins sampler for
burn_in_iter = 5  # number of iterations to invalidate
ld_dict = {}
start_beta = inf_reduced = 0
# boundarys = "genfile"
zero_jump_probability = 0.01  # who knows what the hell this is... its private so ???
sampl_var_shrink_factor = 1 # something to do with the bayes shrink
snp_lrld = None  # We don't need this, we filter out snps that don't meet this criteria
