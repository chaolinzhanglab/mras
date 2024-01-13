import os
import yaml
from pathlib import Path
import pandas as pd
import argparse

# unpack config
configfile: "mapping_and_quantification-config.yaml"
PATHS = config["PATHS"]
PARAMS = config["PARAMS"]

# inputs
INPUT_SAMPLE_INFO_FILE = PATHS["INPUT_SAMPLE_INFO_FILE"]

# outputs
OUTPUT_DIR = PATHS["OUTPUT_DIR"]
OUTPUT_DIR_MAPPING = os.path.join(OUTPUT_DIR,"mapping")
OUTPUT_DIR_QUANT = os.path.join(OUTPUT_DIR,"quantification")

# software
OLEGO_SRC_DIR = PATHS["OLEGO_SRC_DIR"]
OLEGO_INDEX_PATT = PATHS["OLEGO_INDEX_PATT"]
OLEGO_JUNCTIONS_FILE = PATHS["OLEGO_JUNCTIONS_FILE"]
QUANTAS_SRC_DIR = PATHS["QUANTAS_SRC_DIR"]
QUANTAS_ANNOTATION_DIR = PATHS["QUANTAS_ANNOTATION_DIR"]

# general
TMP_DIR = PATHS["TMP_DIR"]
COMBINE_GROUPS = PARAMS["COMBINE_GROUPS"]

# read inputs
SAMPLE_INFO = pd.read_csv(INPUT_SAMPLE_INFO_FILE)
FASTQ_R1_FILES = SAMPLE_INFO["fastq_r1_path"].tolist()
FASTQ_R2_FILES = SAMPLE_INFO["fastq_r2_path"].dropna().tolist()
SAMPLE_NAMES = SAMPLE_INFO["sample_id"].tolist() 
GROUPING = SAMPLE_INFO[["sample_id","group_label"]]

RUN_TYPE = "PAIRED_END" if len(FASTQ_R2_FILES)>0 else "SINGLE_END"

if RUN_TYPE=="PAIRED_END":
    FASTQ_R1S = {s: f for s,f in zip(SAMPLE_NAMES, FASTQ_R1_FILES)}
    FASTQ_R2S = {s: f for s,f in zip(SAMPLE_NAMES, FASTQ_R2_FILES)}
    FASTQS = {"1": FASTQ_R1S, "2": FASTQ_R2S}
else:
    FASTQS = {s: f for s,f in zip(SAMPLE_NAMES, FASTQ_R1_FILES)}


GROUPS = {g: GROUPING.loc[GROUPING["group_label"]==g,"sample_id"].tolist() for g in GROUPING["group_label"].unique()}

EVENT_TYPES = PARAMS["EVENT_TYPES"]

TEST_MODE = PARAMS["TEST_MODE"]

##### RULES #####
rule all:
    input:
        # ----- From .fastq to PSI and expression matrices -----
        # mapping
        ## align separately (if paired)
        expand(os.path.join(OUTPUT_DIR_MAPPING,"sams","{sample}_read{end}.sam"), sample=SAMPLE_NAMES, end=FASTQS.keys()) if RUN_TYPE=="PAIRED_END" else [],
        ## transform alignment to bed
        expand(os.path.join(OUTPUT_DIR_MAPPING,"beds","{sample}.bed"), sample=SAMPLE_NAMES),
        
        # quantification
        expand(os.path.join(OUTPUT_DIR_QUANT,"splicing","{sample}.{data_type}.count.txt"), sample=SAMPLE_NAMES, data_type=EVENT_TYPES),
        expand(os.path.join(OUTPUT_DIR_QUANT,"genexpr","{sample}.expr.txt"), sample=SAMPLE_NAMES),

        # grouping (optional)
        ## combine quantification
        expand(os.path.join(OUTPUT_DIR_QUANT,"splicing_grouped","{group}.{data_type}.count.txt"), group=GROUPS.keys(), data_type=EVENT_TYPES) if COMBINE_GROUPS else [],
        expand(os.path.join(OUTPUT_DIR_QUANT,"genexpr_grouped","{group}.expr.txt"), group=GROUPS.keys()) if COMBINE_GROUPS else [],

        ## create gene expression and splicing matrices
        expand(os.path.join(OUTPUT_DIR_QUANT,"combined","splicing.{data_type}.conf"), data_type=EVENT_TYPES),
        expand(os.path.join(OUTPUT_DIR_QUANT,"combined","{data_type}.psi.txt"), data_type=EVENT_TYPES),
        os.path.join(OUTPUT_DIR_QUANT,"combined","genexpr.conf"),
        os.path.join(OUTPUT_DIR_QUANT,"combined","expr.txt"),

        # cleanup output matrices
        expand(os.path.join(OUTPUT_DIR_QUANT,"combined","{data_type}.psi.clean.txt"), data_type=EVENT_TYPES),
        os.path.join(OUTPUT_DIR_QUANT,"combined","expr.clean.txt"),

        # KNN impute missing values in PSI matrix
        expand(os.path.join(OUTPUT_DIR_QUANT,"imputed","{data_type}.psi.clean.txt"), data_type=EVENT_TYPES),
        
        # testing
        expand(os.path.join("testing","TEST_PASSED-mapping_and_quantification-{run_type}-{data_type}"), run_type=RUN_TYPE.lower(), data_type=EVENT_TYPES) if TEST_MODE else []

        
# ----- Quantify splicing and gene expression for clusters of single cells -----
if RUN_TYPE=="SINGLE_END":
    rule mapping_single_end:
        input:
            fastq = lambda wildcards: FASTQS[wildcards.sample],
            junctions = OLEGO_JUNCTIONS_FILE
        output:
            sam = os.path.join(OUTPUT_DIR_MAPPING,"sams","{sample}.sam"),
            bed = os.path.join(OUTPUT_DIR_MAPPING,"beds","{sample}.bed")
        params:
            src = OLEGO_SRC_DIR,
            genome_index_patt = OLEGO_INDEX_PATT
        threads: 12
        resources:
            runtime = 3600*1, # 1h
            memory = 6 # 6GB
        shell:
            """
            set -eo pipefail

            # align
            {params.src}/olego \
                        --verbose \
                        --num-threads {threads} \
                        --regression-model {params.src}/models/mm.cfg \
                        --junction-file {input.junctions} \
                        --output-file {output.sam} \
                        {params.genome_index_patt} \
                        {input.fastq}
                        
                        
            # convert sam to bed
            perl {params.src}/sam2bed.pl \
                        --uniq \
                        --verbose \
                        {output.sam} \
                        {output.bed}


            echo Done!
            """
            

elif RUN_TYPE=="PAIRED_END":     
    rule mapping_paired_end:
        input:
            fastq = lambda wildcards: FASTQS[wildcards.end][wildcards.sample],
            junctions = OLEGO_JUNCTIONS_FILE
        output:
            sam = os.path.join(OUTPUT_DIR_MAPPING,"sams","{sample}_read{end}.sam")
        params:
            src = OLEGO_SRC_DIR,
            genome_index_patt = OLEGO_INDEX_PATT
        threads: 12
        resources:
            runtime = 3600*1, # 1h
            memory = 6 # 6GB
        shell:
            """
            set -eo pipefail

            # align
            {params.src}/olego \
                        --verbose \
                        --num-threads {threads} \
                        --regression-model {params.src}/models/mm.cfg \
                        --junction-file {input.junctions} \
                        --output-file {output.sam} \
                        {params.genome_index_patt} \
                        {input.fastq}

            echo Done!
            """
            
            
    rule merge_paired_end:
        input:
            sam_r1 = os.path.join(OUTPUT_DIR_MAPPING,"sams","{sample}_read1.sam"),
            sam_r2 = os.path.join(OUTPUT_DIR_MAPPING,"sams","{sample}_read2.sam"),
            isoform = os.path.join(QUANTAS_ANNOTATION_DIR,"mm10.exon.trio.hmr.nr.bed")
        output:
            bed = os.path.join(OUTPUT_DIR_MAPPING,"beds","{sample}.bed")
        params:
            src = QUANTAS_SRC_DIR,
            output_dir = os.path.join(OUTPUT_DIR_MAPPING,"beds","{sample}"),
            split_size = PARAMS["QUANTAS_SPLIT_SIZE"],
            E = PARAMS["QUANTAS_E"]
        threads: 1
        resources:
            runtime = 3600*1, # 1h
            memory = 6 # 6GB
        shell:
            """
            set -eo pipefail
            
            # merge
            perl {params.src}/gapless/gapless_huge_file.pl \
                        -v \
                        -sam \
                        -uniq \
                        --split-size {params.split_size} \
                        -E {params.E} \
                        -big \
                        --print-singleton \
                        -isoform {input.isoform} \
                        -o {params.output_dir} \
                        {input.sam_r1} \
                        {input.sam_r2}
                        
            # rename bed file
            mv {params.output_dir}/pair.gapless.bed {output.bed}

            echo Done!
            """
        
        
    
rule quantify_splicing:
    input:
        as_bed = os.path.join(QUANTAS_ANNOTATION_DIR,"Mm.seq.all.{data_type}.chrom.can.bed"),
        mapping_bed = os.path.join(OUTPUT_DIR_MAPPING,"beds","{sample}.bed")
    output:
        splicing = os.path.join(OUTPUT_DIR_QUANT,"splicing","{sample}.{data_type}.count.txt")
    params:
        data_type = "{data_type}",
        sample = "{sample}",
        run_type = RUN_TYPE,
        src = QUANTAS_SRC_DIR,
        tmp_dir = TMP_DIR,
    threads: 1
    resources:
        runtime = 3600*1, # 1h
        memory = 4 # 4GB
    shell:
        """
        set -eo pipefail
        
        NOW="$(date +%Y%m%d%H%M%S)"
        
        args=(
            -type {params.data_type} 
            -big 
            -v 
            -c {params.tmp_dir}/countit_{params.data_type}_{params.sample}_$NOW 
            {input.as_bed} 
            {input.mapping_bed} 
            {output.splicing}            
        )
        # add "weight" flag when paired-ended reads
        if [[ "{params.run_type}" == "PAIRED_END" ]]; then
            args+=(-weight)
        fi    
        
        perl {params.src}/countit/summarize_splicing.pl "${{args[@]}}"
        
        echo "Done!"
        """
        
        
rule quantify_genexpr:
    input:
        mapping_bed = os.path.join(OUTPUT_DIR_MAPPING,"beds","{sample}.bed"),
        exons = os.path.join(QUANTAS_ANNOTATION_DIR,"mm10.exon.uniq.core.bed"),
        e2g = os.path.join(QUANTAS_ANNOTATION_DIR,"mm10.exon.uniq.core.id2gene2symbol"),
    output:
        genexpr = os.path.join(OUTPUT_DIR_QUANT,"genexpr","{sample}.expr.txt")
    params:
        src = QUANTAS_SRC_DIR,
        run_type = RUN_TYPE,
        tmp_dir = TMP_DIR,
        sample = "{sample}"
    threads: 1
    resources:
        runtime = 3600*1, # 1h
        memory = 4 # 4GB
    shell:
        """
        set -eo pipefail

        NOW="$(date +%Y%m%d%H%M%S)"

        args=(
            -v  
            -big 
            --weight 
            --cache {params.tmp_dir}/countit_expression_{params.sample}_$NOW 
            -exon {input.exons} 
            -e2g {input.e2g} 
            {input.mapping_bed} 
            {output.genexpr}
        )
        # add "weight" flag when paired-ended reads
        if [[ "{params.run_type}" == "PAIRED_END" ]]; then
            args+=(-weight)
        fi    

        perl {params.src}/countit/summarize_expression_wrapper.pl "${{args[@]}}"

        echo "Done!"
        """
    
if COMBINE_GROUPS:
    rule combine_splicing_singlecell_clusters:
        input:
            splicing_samples = lambda wildcards: [
                os.path.join(OUTPUT_DIR_QUANT,"splicing","{sample}.{data_type}.count.txt").format(sample=s, data_type="{data_type}") 
                for s in GROUPS[wildcards.group]
            ]
        output:
            splicing_grouped = os.path.join(OUTPUT_DIR_QUANT,"splicing_grouped","{group}.{data_type}.count.txt")
        params:
            src = QUANTAS_SRC_DIR,
            data_type = "{data_type}",
            group = "{group}",
            tmp_dir = TMP_DIR
        threads: 1
        resources:
            runtime = 3600*1, # 1h
            memory = 4 # 4GB
        shell:
            """
            set -eo pipefail

            NOW="$(date +%Y%m%d%H%M%S)"

            perl {params.src}/countit/combine_replicates.pl \
                        -v \
                        -type {params.data_type} \
                        -cache {params.tmp_dir}/countit_combine_{params.data_type}_{params.group}_$NOW \
                        {input.splicing_samples} \
                        {output.splicing_grouped}

            echo "Done!"
            """


    rule combine_genexpr_singlecell_clusters:
        input:
            genexpr_samples = lambda wildcards: [
                os.path.join(OUTPUT_DIR_QUANT,"genexpr","{sample}.expr.txt").format(sample=s) 
                for s in GROUPS[wildcards.group]
            ]
        output:
            genexpr_grouped = os.path.join(OUTPUT_DIR_QUANT,"genexpr_grouped","{group}.expr.txt")
        params:
            src = QUANTAS_SRC_DIR,
            data_type = "expr",
            group = "{group}",
            tmp_dir = TMP_DIR
        threads: 1
        resources:
            runtime = 3600*1, # 1h
            memory = 4 # 6GB
        shell:
            """
            set -eo pipefail

            NOW="$(date +%Y%m%d%H%M%S)"

            perl {params.src}/countit/combine_replicates.pl \
                        -v \
                        -type {params.data_type} \
                        -cache {params.tmp_dir}/countit_combine_{params.data_type}_{params.group}_$NOW \
                        {input.genexpr_samples} \
                        {output.genexpr_grouped}

            echo "Done!"
            """


rule create_conf_splicing:
    input:
        splicing = [os.path.join(OUTPUT_DIR_QUANT,"splicing_grouped","{group}.{data_type}.count.txt").format(group=g, data_type="{data_type}") for g in GROUPS.keys()] if COMBINE_GROUPS else [os.path.join(OUTPUT_DIR_QUANT,"splicing","{sample}.{data_type}.count.txt").format(sample=s, data_type="{data_type}") for s in FASTQS.keys()]
    output:
        conf = os.path.join(OUTPUT_DIR_QUANT,"combined","splicing.{data_type}.conf")
    params:
        file_suffix = ".{data_type}.count.txt"
    threads: 1
    resources:
        runtime = int(3600*0.5), # 30 min
        memory = 1 # 1 GB
    run:
        import pandas as pd
        
        conf = []
        for f in input.splicing:
            filename = f
            file_id = os.path.basename(filename).replace(params.file_suffix,"")
            
            conf.append({
                "filename": filename,
                "file_id": file_id
            })
        conf = pd.DataFrame(conf)
        
        conf.to_csv(output.conf, sep="\t", index=False, header=False)
        
        print("Done!")
        
        
rule create_conf_genexpr:
    input:
        genexpr = [os.path.join(OUTPUT_DIR_QUANT,"genexpr_grouped","{group}.expr.txt").format(group=g) for g in GROUPS.keys()] if COMBINE_GROUPS else [os.path.join(OUTPUT_DIR_QUANT,"genexpr","{sample}.expr.txt").format(sample=s) for s in FASTQS.keys()]
    output:
        conf = os.path.join(OUTPUT_DIR_QUANT,"combined","genexpr.conf")
    params:
        file_suffix = ".expr.txt"
    threads: 1
    resources:
        runtime = int(3600*0.5), # 30 min
        memory = 1 # 1 GB
    run:
        import pandas as pd
        
        conf = []
        for f in input.genexpr:
            filename = f
            file_id = os.path.basename(filename).replace(params.file_suffix,"")
            
            conf.append({
                "filename": filename,
                "file_id": file_id
            })
        conf = pd.DataFrame(conf)
        
        conf.to_csv(output.conf, sep="\t", index=False, header=False)
        
        print("Done!")
        
        
rule make_splicing_matrix:
    input:
        splicing = [os.path.join(OUTPUT_DIR_QUANT,"splicing_grouped","{group}.{data_type}.count.txt").format(group=g, data_type="{data_type}") for g in GROUPS.keys()] if COMBINE_GROUPS else [os.path.join(OUTPUT_DIR_QUANT,"splicing","{sample}.{data_type}.count.txt").format(sample=s, data_type="{data_type}") for s in FASTQS.keys()],
        conf = os.path.join(OUTPUT_DIR_QUANT,"combined","splicing.{data_type}.conf"),
        id2gene2symbol = os.path.join(QUANTAS_ANNOTATION_DIR,"Mm.seq.all.AS.chrom.can.id2gene2symbol")
    output:
        splicing_mat = os.path.join(OUTPUT_DIR_QUANT,"combined","{data_type}.psi.txt")
    params:
        src = QUANTAS_SRC_DIR,
        data_type = "{data_type}",
        min_cov = PARAMS["QUANTAS_MIN_COV"],
        max_std = PARAMS["QUANTAS_MAX_STD"],
        na_string = "NaN"
    threads: 1
    resources:
        runtime = int(3600*0.5), # 30 min
        memory = 1 # 1 GB
    shell:
        """
        set -eo pipefail
        
        perl {params.src}/countit/gen_splicing_matrix.pl \
                    -v \
                    -type {params.data_type} \
                    --min-cov {params.min_cov} \
                    --max-std {params.max_std} \
                    --na-string {params.na_string} \
                    --print-info \
                    --id2gene2symbol {input.id2gene2symbol} \
                    {input.conf} \
                    {output.splicing_mat}
                    
        echo "Done!"
        """
        
        
rule make_genexpr_matrix:
    input:
        genexpr = [os.path.join(OUTPUT_DIR_QUANT,"genexpr_grouped","{group}.expr.txt").format(group=g) for g in GROUPS.keys()] if COMBINE_GROUPS else [os.path.join(OUTPUT_DIR_QUANT,"genexpr","{sample}.expr.txt").format(sample=s) for s in FASTQS.keys()],
        conf = os.path.join(OUTPUT_DIR_QUANT,"combined","genexpr.conf")
    output:
        genexpr_mat = os.path.join(OUTPUT_DIR_QUANT,"combined","expr.txt")
    params:
        src = QUANTAS_SRC_DIR,
        na_string = "NaN"
    threads: 1
    resources:
        runtime = int(3600*0.5), # 30 min
        memory = 1 # 1 GB
    shell:
        """
        set -eo pipefail
        
        perl {params.src}/countit/gen_expression_matrix.pl \
                    -v \
                    -pseudocount 1 \
                    -log2 \
                    {input.conf} \
                    {output.genexpr_mat}
                    
        echo "Done!"
        """
        
        
rule cleanup_splicing:
    input:
        splicing_mat = os.path.join(OUTPUT_DIR_QUANT,"combined","{data_type}.psi.txt")
    output:
        splicing_mat = os.path.join(OUTPUT_DIR_QUANT,"combined","{data_type}.psi.clean.txt")
    threads: 1
    resources:
        runtime = int(3600*0.5), # 30 min
        memory = 10 # 10 GB
    run:
        import pandas as pd
        
        # load
        splicing = pd.read_table(input.splicing_mat)
        
        # drop splicing columns
        to_drop = ["#chrom","chromStart","chromEnd","score","strand","type","isoforms","NAME"]
        splicing = splicing.drop(columns=to_drop).copy()
        
        # save
        splicing.to_csv(output.splicing_mat, sep="\t", index=False)

        print("Done!")
        

rule cleanup_genexpr:
    input:
        genexpr_mat = os.path.join(OUTPUT_DIR_QUANT,"combined","expr.txt")
    output:
        genexpr_mat = os.path.join(OUTPUT_DIR_QUANT,"combined","expr.clean.txt")
    threads: 1
    resources:
        runtime = int(3600*0.5), # 30 min
        memory = 10 # 10 GB
    run:
        import pandas as pd
        
        # load
        genexpr = pd.read_table(input.genexpr_mat)
        
        # combine genexpr columns
        genexpr["gene_id"] = genexpr["gene_id"].astype("str")+"//"+genexpr["NAME"]
        genexpr = genexpr.drop(columns=["NAME"])
        
        # save
        genexpr.to_csv(output.genexpr_mat, sep="\t", index=False)
        
        print("Done!")
        
        
rule knn_imputation_splicing:
    input:
        splicing_mat = os.path.join(OUTPUT_DIR_QUANT,"combined","{data_type}.psi.clean.txt"),
    output:
        splicing_mat = os.path.join(OUTPUT_DIR_QUANT,"imputed","{data_type}.psi.clean.txt")
    params:
        n_neighbors = PARAMS["KNN_N_NEIGHBORS"],
        rowmax = PARAMS["KNN_ROWMAX"],
        colmax = PARAMS["KNN_COLMAX"],
        maxp = PARAMS["KNN_MAXP"],
        rng_seed = PARAMS["KNN_SEED"]
    threads: 1
    resources:
        runtime = 3600*1, # 1 h
        memory = 20 # 20 GB
    shell:
        """
        set -eo pipefail
        
        Rscript src/impute_knn.R \
                    --input_file={input.splicing_mat} \
                    --n_neighbors={params.n_neighbors} \
                    --rowmax={params.rowmax} \
                    --colmax={params.colmax} \
                    --maxp={params.maxp} \
                    --rng_seed={params.rng_seed} \
                    --output_file={output.splicing_mat}
                    
        echo "Done!"
        """
    
    
if TEST_MODE:
    rule test_reproducibility:
        input:
            genexpr = os.path.join(OUTPUT_DIR_QUANT,"combined","expr.clean.txt"),
            splicing = os.path.join(OUTPUT_DIR_QUANT,"imputed","{data_type}.psi.clean.txt"),
            genexpr_ref = PATHS["TESTING_REF_GENEXPR_FILE"],
            splicing_ref = PATHS["TESTING_REF_SPLICING_FILE"]
        output:
            touch(os.path.join("testing","TEST_PASSED-mapping_and_quantification-{run_type}-{data_type}"))
        run:
            import pandas as pd
            
            # load new matrices
            genexpr = pd.read_table(input.genexpr)
            splicing = pd.read_table(input.splicing)
            
            # load testing references
            genexpr_ref = pd.read_table(input.genexpr_ref)
            splicing_ref = pd.read_table(input.splicing_ref)
            
            # test they are equal
            pd.testing.assert_frame_equal(genexpr, genexpr_ref)
            pd.testing.assert_frame_equal(splicing, splicing_ref)
            
            print("Done!")