import os
CTAT_GENOME_LIB =  "/media/nadine/HOME/nadine/ctat-mutations_test/ctat_genome_lib_build_dir/"

print(CTAT_GENOME_LIB)
rule all:
    input:
        "allcatch_output/predictions.tsv",
        "ctat_output/variants.HC_init.wAnnot.vcf.gz"

rule install_allcatchr:
    shell:
        "Rscript -e 'devtools::install_github(\"ThomasBeder/ALLCatchR\")'"

# Rule to run ALLCatchR
rule run_allcatchr:
    input:
        r_script="scripts/run_ALLCatchR.R",
        input_file="data/counts/UKSH_counts.tsv",  # Update with the correct path
    output:
        "allcatch_output/predictions.tsv"  # TODO: Is it possible to create an output directory for that?

    shell:
        "Rscript {input.r_script} {input.input_file}  {output};"
        "mv predictions.tsv {output}"


rule clone_ctat_mutations_repo:
    shell:
        "git clone --recursive https://github.com/NCIP/ctat-mutations.git"


# Rule to run ctat-mutations
rule run_ctat_mutations:
    input:
        left="data/samples/reads_1.fastq.gz",# Update with the correct path
        right="data/samples/reads_2.fastq.gz",# Update with the correct path
    params:
        sample_id="test001",# Replace with the actual sample ID
        outputdir="ctat_output",  # Replace with the desired output directory name
        genome_lib = CTAT_GENOME_LIB

    output:
        "ctat_output/variants.HC_init.wAnnot.vcf.gz"  # Update with a proper output file


    shell:
            """
            python ctat-mutations/ctat_mutations \
            --left {input.left} \
            --right {input.right} \
            --outputdir {params.outputdir} \
            --sample_id {params.sample_id} \
            --genome_lib_dir {params.genome_lib}\
            --boosting_method none \
            --no_cravat
            """




