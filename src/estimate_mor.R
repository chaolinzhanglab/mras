#
# Author: Miquel Anglada-Girotto and Daniel Moakley
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# Impute missing values with K-Nearest Neighbors algorithm.

require(optparse)
require(tidyverse)
require(viper)

# Development
# -----------
ROOT = here::here()
OUTPUT_DIR_QUANT = file.path(ROOT,"results","aracne_splicing","paired_end","quantification")
OUTPUT_DIR_ARACNEAP = file.path(ROOT,"results","aracne_splicing","paired_end","aracne")
splicing_file = file.path(OUTPUT_DIR_QUANT,"imputed","cass.psi.clean.txt")
genexpr_file = file.path(OUTPUT_DIR_QUANT,"combined","expr.clean.txt")
aracne_network_file = file.path(OUTPUT_DIR_ARACNEAP,"regulons","cass","network.txt")

splicing_file = "../src/ARACNe-AP/test/matrix.txt"
genexpr_file = "../src/ARACNe-AP/test/matrix.txt"
aracne_network_file = "../data/examples/network_reverse_engineering/network.txt"

##### FUNCTIONS #####
parseargs = function(){
    
    option_list = list( 
        make_option("--splicing_file", type="character"),
        make_option("--genexpr_file", type="integer", default=10),
        make_option("--aracne_network_file", type="numeric", default=0.5),
        make_option("--output_file", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}


main = function(){
    args = parseargs()
    
    splicing_file = args[["splicing_file"]]
    genexpr_file = args[["genexpr_file"]]
    aracne_network_file = args[["aracne_network_file"]]
    output_file = args[["output_file"]]
    
    # load
    genexpr = read_tsv(genexpr_file)
    splicing = read_tsv(splicing_file)
    
    # prep
    genexpr = genexpr %>% column_to_rownames(colnames(genexpr)[1])
    splicing = splicing %>% column_to_rownames(colnames(splicing)[1])
    X = genexpr %>% bind_rows(splicing) %>% distinct() %>% as.matrix()
    
    # make regulons
    regulons = aracne2regulon(aracne_network_file, X, format="adj")
    
    # save
    saveRDS(regulons, output_file)
}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}