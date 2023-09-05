import pandas as pd
import os

#CTAT_GENOME_LIB =  "/media/nadine/HOME/nadine/ctat-mutations_test/ctat_genome_lib_build_dir"
configfile: "config.yaml"


rule all:
    input:
        "allcatch_output/predictions.tsv",
        "ctat_output_directory/test007.cancer.vcf"

rule install_allcatchr:
    shell:
        "Rscript -e 'devtools::install_github(\"ThomasBeder/ALLCatchR\")'"

# Rule to run ALLCatchR
rule run_allcatchr:
    input:
        r_script = "scripts/run_ALLCatchR.R",
        input_file = config["counts"],  # Update with the correct path
    output:
        "allcatch_output/predictions.tsv"

    shell:
        "Rscript {input.r_script} {input.input_file} {output};"
        "mv predictions.tsv {output}"


rule pull_ctat_mutations_singularity_image:
    shell:
        "singularity pull docker://trinityctat/ctat_mutations"


# Rule to run ctat-mutations
rule run_ctat_mutations:
    input:
        mount_dir = config["ctat_mount_dir"]
        # left="/home/nadine/PycharmProjects/Blast-o-Matic-Fusioninator/data/samples/reads_1.fastq.gz",  # Update with the correct path
        # right="/home/nadine/PycharmProjects/Blast-o-Matic-Fusioninator/data/samples/reads_2.fastq.gz"  # Update with the correct path
        #TODO: How can I pass here the input left and right samples?

    params:
        sample_id = "test007", # TODO Replace with the actual sample ID
        output_dir = config["ctat_output_dir"],  # Replace with the desired output directory name
        genome_lib = config["genome_lib"],
        left = config["left_samples"],
        right = config["right_samples"],

    output:
        "ctat_output_directory/{sample_id}.cancer.vcf"

    shell:
         """
         singularity exec -e -B {input.mount_dir}:/data \
         -B {params.genome_lib}:/ctat_genome_lib_dir \
         ctat_mutations.v3.2.0.simg \
         /usr/local/src/ctat-mutations/ctat_mutations \
         --left {input.mount_dir}{params.left} \
         --right {input.mount_dir}{params.right} \
         --sample_id {params.sample_id} \
         --output {params.output_dir}  \
         --cpu 10 \
         --genome_lib_dir /ctat_genome_lib_dir \
         --boosting_method none \
         --no_cravat
         """



