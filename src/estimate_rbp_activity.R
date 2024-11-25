#
# Author: Miquel Anglada-Girotto and Daniel Moakley
# Contact: miquel [dot] anglada [at] crg [dot] eu
#
# Script purpose
# --------------
# Estimate sample-level RBP activity using the VIPER algorithm.

require(optparse)
require(tidyverse)
require(viper)

##### FUNCTIONS #####
parseargs = function(){
    
    option_list = list( 
        make_option("--splicing_file", type="character"),
        make_option("--output_folder", type="character")
    )

    args = parse_args(OptionParser(option_list=option_list))
    
    return(args)
}


main = function(){
    args = parseargs()
    
    splicing_file = args[["splicing_file"]]
    output_folder = args[["output_folder"]]
   
    # load
    splicing = read_tsv(splicing_file)
    splicing = splicing %>% column_to_rownames(colnames(splicing)[1])
   
    loc_reg = paste(output_folder,"regulon.rds",sep="/")
    regulons<-readRDS(loc_reg)

    # estimate RBP activity 
    rbp_activity = viper(splicing, regulons, method="scale", pleiotropy=F, verbose=T)
    output_act = paste(output_folder,"rbpactivity.txt",sep="/")
    write.table(rbp_activity, output_act, sep="\t")

}


##### SCRIPT #####
if (sys.nframe() == 0L) {
    main()
    print("Done!")
}
