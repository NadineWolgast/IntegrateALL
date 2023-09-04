import os
CTAT_GENOME_LIB =  "/media/nadine/HOME/nadine/ctat-mutations_test/ctat_genome_lib_build_dir"

print(CTAT_GENOME_LIB)
rule all:
    input:
        "allcatch_output/predictions.tsv",
        "ctat_output_directory/test003.cancer.vcf"

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


#rule clone_ctat_mutations_repo:
#    shell:
#        "git clone --recursive https://github.com/NCIP/ctat-mutations.git"

rule pull_ctat_mutations_singularity_image:
    shell:
        "singularity pull docker://trinityctat/ctat_mutations"


# Rule to run ctat-mutations
rule run_ctat_mutations:
    input:
        left = "/home/nadine/PycharmProjects/Blast-o-Matic-Fusioninator/data/samples/reads_1.fastq.gz",  # Update with the correct path
        right = "/home/nadine/PycharmProjects/Blast-o-Matic-Fusioninator/data/samples/reads_2.fastq.gz"  # Update with the correct path

    params:
        sample_id="test003",# Replace with the actual sample ID
        outputdir="/home/nadine/PycharmProjects/Blast-o-Matic-Fusioninator/ctat_output_directory",  # Replace with the desired output directory name
        #genome_lib = CTAT_GENOME_LIB
    output:
        "ctat_output_directory/{sample_id}.cancer.vcf"  # Update with a proper output file

    shell:
        """
        singularity exec -e -B /home/nadine/PycharmProjects/Blast-o-Matic-Fusioninator/ctat-mutations/testing:/data \
        -B /media/nadine/HOME/nadine/ctat-mutations_test/ctat_genome_lib_build_dir:/ctat_genome_lib_dir \
        /media/nadine/HOME/nadine/ctat-mutations_test/ctat_mutations.v3.2.0.simg \
        /usr/local/src/ctat-mutations/ctat_mutations \
        --left /data/reads_1.fastq.gz \
        --right /data/reads_2.fastq.gz \
        --sample_id {params.sample_id} \
        --output /home/nadine/PycharmProjects/Blast-o-Matic-Fusioninator/ctat_output_directory \
        --cpu 10 \
        --genome_lib_dir /ctat_genome_lib_dir \
        --boosting_method none \
        --no_cravat
        """
