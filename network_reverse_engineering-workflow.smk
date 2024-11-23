import os
import yaml
from pathlib import Path
import argparse

# unpack config
configfile: "network_reverse_engineering-config.yaml"
PATHS = config["PATHS"]
PARAMS = config["PARAMS"]

# inputs
INPUT_GENEXPR_FILE = PATHS["INPUT_GENEXPR_FILE"]
INPUT_SPLICING_FILES = PATHS["INPUT_SPLICING_FILES"]
INPUT_REGULATORS_OI_FILE = PATHS["INPUT_REGULATORS_OI_FILE"]

# outputs
OUTPUT_DIR = PATHS["OUTPUT_DIR"]
OUTPUT_DIR_ARACNEAP = os.path.join(OUTPUT_DIR,"aracne")
OUTPUT_DIR_RBPACT = os.path.join(OUTPUT_DIR,"rbpactivity")

# software
ARACNEAP_SRC_DIR = PATHS["ARACNEAP_SRC_DIR"]

# general
TMP_DIR = PATHS["TMP_DIR"]

EVENT_TYPES = list(INPUT_SPLICING_FILES.keys())

ESTIMATE_ACTIVITY = PARAMS["ESTIMATE_ACTIVITY"]

TEST_MODE = PARAMS["TEST_MODE"]

# ARACNe-AP
N_BOOTSTRAPS = PARAMS["ARACNE_N_BOOTSTRAPS"]
BOOTSTRAPS = list(range(N_BOOTSTRAPS))
##### RULES #####
rule all:
    input:
        # ----- ARACNe RBP-cass network reverse engineering -----
        expand(os.path.join(OUTPUT_DIR_ARACNEAP,"regulons","{data_type}",".done","done_threshold"), data_type=EVENT_TYPES),
        expand(os.path.join(OUTPUT_DIR_ARACNEAP,"regulons","{data_type}",".done","done_bootstrap_{boot_i}"), data_type=EVENT_TYPES, boot_i=BOOTSTRAPS),
        expand(os.path.join(OUTPUT_DIR_ARACNEAP,"regulons","{data_type}",".done","done_pruning"), data_type=EVENT_TYPES),
        expand(os.path.join(OUTPUT_DIR_ARACNEAP,"regulons","{data_type}","pruned","network.txt"), data_type=EVENT_TYPES),
        expand(os.path.join(OUTPUT_DIR_RBPACT,"regulons","{data_type}","regulon.rds"), data_type=EVENT_TYPES),
        expand(os.path.join(OUTPUT_DIR_RBPACT,"regulons","{data_type}","rbpactivity.txt"), data_type=EVENT_TYPES),
        
        # testing
        expand(os.path.join("testing","TEST_PASSED-network_reverse_engineering-{data_type}"), data_type=EVENT_TYPES) if TEST_MODE else []

        
# ----- ARACNe RBP-cass network reverse engineering -----
rule regulon_inference_aracne_java_threshold:
    input:
        regulators = INPUT_GENEXPR_FILE,
        targets = lambda wildcards: INPUT_SPLICING_FILES[wildcards.data_type],
        regulators_oi = INPUT_REGULATORS_OI_FILE
    output:
        touch(os.path.join(OUTPUT_DIR_ARACNEAP,"regulons","{data_type}",".done","done_threshold"))
    params:
        src = ARACNEAP_SRC_DIR,
        output_dir = os.path.join(OUTPUT_DIR_ARACNEAP,"regulons","{data_type}"),
        random_seed = PARAMS["ARACNE_MI_THRESH_SEED"],
        mi_pvalue_thresh = str(PARAMS["ARACNE_MI_THRESH_PVALUE"]).replace("e","E"),
    threads: 12
    resources:
        runtime = 3600*12, # 12h 
        memory = 5, # 5GB
    shell:
        """
        set -eo pipefail
        
        java -Xmx{resources.memory}G -jar {params.src}/dist/aracne.jar \
                --expfile_upstream {input.regulators} \
                --tfs {input.regulators_oi} \
                --expfile_downstream {input.targets} \
                --output {params.output_dir} \
                --pvalue {params.mi_pvalue_thresh} \
                --seed {params.random_seed} \
                --threads {threads} \
                --calculateThreshold
                
        echo "Done!"
        """
        
        
rule regulon_inference_aracne_java_bootstrap:
    input:
        os.path.join(OUTPUT_DIR_ARACNEAP,"regulons","{data_type}",".done","done_threshold"),
        regulators = INPUT_GENEXPR_FILE,
        targets = lambda wildcards: INPUT_SPLICING_FILES[wildcards.data_type],
        regulators_oi = INPUT_REGULATORS_OI_FILE
    output:
        touch(os.path.join(OUTPUT_DIR_ARACNEAP,"regulons","{data_type}",".done","done_bootstrap_{boot_i}"))
    params:
        src = ARACNEAP_SRC_DIR,
        output_dir = os.path.join(OUTPUT_DIR_ARACNEAP,"regulons","{data_type}"),
        random_seed = "{boot_i}",
        mi_pvalue_thresh = str(PARAMS["ARACNE_MI_THRESH_PVALUE"]).replace("e","E")
    threads: 6
    resources:
        runtime = 3600*12, # 12h 
        memory = 5, # 5GB
    shell:
        """
        set -eo pipefail

        java -Xmx{resources.memory}G -jar {params.src}/dist/aracne.jar \
                --expfile_upstream {input.regulators} \
                --tfs {input.regulators_oi} \
                --expfile_downstream {input.targets} \
                --output {params.output_dir} \
                --pvalue {params.mi_pvalue_thresh} \
                --seed {params.random_seed} \
                --threads {threads} \
                --nodpi
                
        echo "Done!"
        """

        
rule regulon_inference_aracne_prune_bootstraps:
    input:
        os.path.join(OUTPUT_DIR_ARACNEAP,"regulons","{data_type}",".done","done_threshold"),
        [os.path.join(OUTPUT_DIR_ARACNEAP,"regulons","{data_type}",".done","done_bootstrap_{boot_i}").format(data_type="{data_type}", boot_i=boot_i) for boot_i in BOOTSTRAPS]
    output:
        touch(os.path.join(OUTPUT_DIR_ARACNEAP,"regulons","{data_type}",".done","done_pruning"))
    params:
        output_dir = os.path.join(OUTPUT_DIR_ARACNEAP,"regulons","{data_type}"),
        max_targets = PARAMS["ARACNE_BOOTSTRAP_MAX_TARGETS"]
    threads: 1
    resources:
        runtime = 3600*1, # 1 h
        memory = 10 # 10 GB
    shell:
        """
        set -eo pipefail

        perl src/filter_arachne_bootstraps.pl \
                    --regulon-size {params.max_targets} \
                    {params.output_dir} \
                    {params.output_dir}/pruned

        echo "Done!"
        """


rule regulon_inference_aracne_java_consolidation:
    input:
        os.path.join(OUTPUT_DIR_ARACNEAP,"regulons","{data_type}",".done","done_threshold"),
        [os.path.join(OUTPUT_DIR_ARACNEAP,"regulons","{data_type}",".done","done_bootstrap_{boot_i}").format(data_type="{data_type}", boot_i=boot_i) for boot_i in BOOTSTRAPS],
        os.path.join(OUTPUT_DIR_ARACNEAP,"regulons","{data_type}",".done","done_pruning")
    output:
        os.path.join(OUTPUT_DIR_ARACNEAP,"regulons","{data_type}","pruned","network.txt")
    params:
        src = ARACNEAP_SRC_DIR,
        output_dir = os.path.join(OUTPUT_DIR_ARACNEAP,"regulons","{data_type}","pruned"),
        consolidate_pvalue = PARAMS["ARACNE_CONSOLIDATE_PVALUE"]
    threads: 12
    resources:
        runtime = 3600*12, # 12h 
        memory = 20, # 15GB
    shell:
        """
        set -eo pipefail
        
        java -Xmx{resources.memory}G -jar {params.src}/dist/aracne.jar \
                --output {params.output_dir} \
                --threads {threads} \
                --consolidatepvalue {params.consolidate_pvalue} \
                --consolidate
                
        echo "Done!"
        """
        

rule prepare_regulons:
    input:
        genexpr = INPUT_GENEXPR_FILE,
        splicing = lambda wildcards: INPUT_SPLICING_FILES[wildcards.data_type],
        aracne_network = os.path.join(OUTPUT_DIR_ARACNEAP,"regulons","{data_type}","pruned","network.txt"),
    output:
        os.path.join(OUTPUT_DIR_RBPACT,"regulons","{data_type}","regulon.rds"),
    params:
        output_dir = os.path.join(OUTPUT_DIR_RBPACT,"regulons","{data_type}")

    threads: 1
    resources:
        runtime = 3600*1, # 1 h
        memory = 10 # GB
    shell:
        """
        set -eo pipefail
        
        Rscript src/estimate_mor.R \
                    --splicing_file={input.splicing} \
                    --genexpr_file={input.genexpr} \
                    --aracne_network_file={input.aracne_network} \
                    --output_folder={params.output_dir}
                    
        echo "Done! estimating MOR!"
        """

if ESTIMATE_ACTIVITY:
    rule estimate_rbp_activity:
        input:
            splicing = lambda wildcards: INPUT_SPLICING_FILES[wildcards.data_type],
        output:
             os.path.join(OUTPUT_DIR_RBPACT,"regulons","{data_type}","rbpactivity.txt"),
        params:
            output_dir = os.path.join(OUTPUT_DIR_RBPACT,"regulons","{data_type}")
        threads: 1
        resources:
            runtime = 3600*1, # 1 h
            memory = 10 # GB
        shell:
            """
            set -eo pipefail
        
            Rscript src/estimate_rbp_activity.R \
                        --splicing_file={input.splicing} \
                        --output_folder={params.output_dir}
                    
            echo "Done estimating RBP activity!"
            """

if TEST_MODE:
    rule test_reproducibility:
        input:
            network = os.path.join(OUTPUT_DIR_ARACNEAP,"regulons","{data_type}","pruned","network.txt"),
            network_ref = PATHS["TESTING_REF_NETWORK_FILE"]
        output:
            touch(os.path.join("testing","TEST_PASSED-network_reverse_engineering-{data_type}"))
        run:
            import pandas as pd
            
            network = pd.read_table(input.network)
            network_ref = pd.read_table(input.network_ref)
            
            X = pd.merge(network, network_ref, on=["Target","Regulator"], how="inner", suffixes=["_new","_ref"])
            correlation = X["MI_new"].corr(X["MI_ref"], method="pearson")
            
            print("Correlation is:", correlation)
            assert correlation>0.99
            
            print("Done!")
