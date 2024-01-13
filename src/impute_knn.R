#
# Author: Miquel Anglada Girotto
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# Impute missing values with K-Nearest Neighbors algorithm.

require(optparse)
require(tidyverse)
require(impute)

# Development
# -----------
# ROOT = here::here()
# input_file = "../data/examples/testing/paired_end/results/mapping_and_quantification/quantification/combined/cass.psi.clean.txt"
# n_neighbors = 10
# rowmax = 1
# colmax = 1
# maxp = 1500
# rng.seed = 362436069

##### FUNCTIONS #####
parseargs = function(){
    
    option_list = list( 
        make_option("--input_file", type="character"),
        make_option("--n_neighbors", type="integer", default=10),
        make_option("--rowmax", type="numeric", default=0.5),
        make_option("--colmax", type="numeric", default=0.8),
        make_option("--maxp", type="integer", default=1500),
        make_option("--rng_seed", type="integer", default=362436069),
        make_option("--output_file", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}


main = function(){
    args = parseargs()
    
    input_file = args[["input_file"]]
    n_neighbors = args[["n_neighbors"]]
    rowmax = args[["rowmax"]]
    colmax = args[["colmax"]]
    maxp = args[["maxp"]]
    rng_seed = args[["rng_seed"]]
    output_file = args[["output_file"]]
    
    # load
    X = read_tsv(input_file)
    
    # prep
    ## set first column as row names
    first_col_name = colnames(X)[1]
    X = X %>%
        as.data.frame() %>%
        column_to_rownames(first_col_name)
    
    ## drop fully missing rows
    X = X %>% 
        filter_all(any_vars(!is.na(.))) %>%
        as.matrix()
    
    # impute missing values
    result = impute.knn(
        data = X,
        k = n_neighbors,
        rowmax = rowmax,
        colmax = colmax,
        maxp = maxp,
        rng.seed = rng_seed
    )[["data"]]
    
    # save
    result = result %>% as.data.frame() %>% rownames_to_column(first_col_name)
    write_tsv(result, output_file)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}