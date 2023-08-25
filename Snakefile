#rule all:
#    input:
#        "predictions.tsv"

#rule install_allcatchr:
#    shell:
#        "Rscript -e 'devtools::install_github(\"ThomasBeder/ALLCatchR\")'"

# Rule to run ALLCatchR
#rule run_allcatchr:
#    input:
#        r_script="scripts/run_ALLCatchR.R",
#        input_file="data/counts/UKSH_counts.tsv",  # Update with the correct path
#    output:
#        "predictions.csv"  # TODO: Is it possible to create an output directory for that?
#    shell:
#        "Rscript {input.r_script} {input.input_file}  {output.output_file}"

#rule pull_ctat_mutations_docker:
#    shell:
#        "docker pull trinityctat/ctat_mutations"


#rule run_docker_ctat_mutations:
#    input:
#        left="data/samples/reads_1.fastq.gz",  # Update with the correct path
#        right="data/samples/reads_2.fastq.gz",  # Update with the correct path
#    params:
#        sample_id="test",
#        outputdir="ctat_mutations_outdir",
#        cpu=10,
#        ctat_genome_lib="/media/nadine/HOME/nadine/ctat-mutations_test/ctat_genome_lib_build_dir/ctat_genome_lib"  # Update with the correct path
#    output:
#        touch("ctat_mutations_completed.txt")
#    shell:
        """
        DATA_FOLDER=`pwd`
        docker run --rm -v ${DATA_FOLDER}:${DATA_FOLDER} -v /tmp:/tmp \
            -v {params.ctat_genome_lib}:/ctat_genome_lib_dir:ro \
            trinityctat/ctat_mutations \
            /usr/local/src/ctat-mutations/ctat_mutations \
            --left {input.left} \
            --right {input.right} \
            --sample_id {params.sample_id} \
            --output {params.outputdir} \
            --cpu {params.cpu} \
            --genome_lib_dir /ctat_genome_lib_dir
        touch {output}
        """



''' funktioniert noch nicht
# Rule to install ctat-mutations using Conda
rule install_ctat_mutations:
    conda:
        "envs/ctat.yaml"  # Conda environment file that specifies dependencies
    shell:
        "conda env create -n ctat_mutations_env --file envs/ctat.yaml"


# Rule to run ctat-mutations
rule run_ctat_mutations:
    input:
        left="path/to/left.fastq",# Update with the correct path
        right="path/to/right.fastq",# Update with the correct path
    params:
        sample_id=config["samples"],# Replace with the actual sample ID
        outputdir="output_directory"  # Replace with the desired output directory name

    output:
        touch("path/to/ctat_mutations_completed.txt")  # Update with a proper output file
    #params:
    #    conda_env="ctat_mutations_env"  # Name of the Conda environment
    #conda:
    #    "envs/ctat.yaml"  # Conda environment file that specifies dependencies
    shell:
        "          python /path/to/ctat_mutations \\n"
        "          --left {input.left} \\n"
        "          --right {input.right} \\n"
        "          --sample_id {params.sample_id} \\n"
        "          --outputdir {params.outputdir}\n"
        "          touch {output}\n"

        #"source activate {params.conda_env} && ctat-mutations {input.input_bam} {input.reference} {output.output_vcf}"
        # Activate the Conda environment and run ctat-mutations

 '''
