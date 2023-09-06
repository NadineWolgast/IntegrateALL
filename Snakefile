configfile: "config.yaml"


rule all:
    input:
        "allcatch_output/predictions.tsv",
        "ctat_output_directory/test007.cancer.vcf",
        "data/config.txt",
        "data/meta.txt"


#rule star_pe_multi:
#    input:
#        # use a list for multiple fastq files for one sample
#        # usually technical replicates across lanes/flowcells
#        fq1=["data/samples/{sample}_1.fastq", "data/samples/{sample}_2.fastq"],
#        # paired end reads needs to be ordered so each item in the two lists match
#        #fq2=["data/samples/{sample}_R2.1.fastq", "data/samples/{sample}_R2.2.fastq"],  #optional
#        # path to STAR reference genome index
#        idx="index",
#    output:
#        # see STAR manual for additional output files
#        aln="star/pe/{sample}/pe_aligned.sam",
#        log="logs/pe/{sample}/Log.out",
#        sj="star/pe/{sample}/SJ.out.tab",
#        unmapped=["star/pe/{sample}/unmapped.1.fastq.gz","star/pe/{sample}/unmapped.2.fastq.gz"],
#    log:
 #       "logs/pe/{sample}.log",
#    params:
#        # optional parameters
#        extra="",
#    threads: 8
#    wrapper:
#        "v2.6.0/bio/star/align"


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
        right = config["right_samples"]

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


rule install_rnaseq_cnv:
    shell:
        "Rscript -e 'devtools::install_github(\"honzee/RNAseqCNV\")'"


rule write_rnaseq_cnv_config_file:
    input:
        r_script = "scripts/write_rnaseq_cnv_config.R",
        count_dir = config["counts"],
        snv_dir = "ctat_output_directory/",
        out_dir= "rnaseq_cnv_output_directory"

    output:
        "data/config.txt"

    shell:
        "Rscript {input.r_script} {input.out_dir} {input.count_dir} {input.snv_dir} {output};"
        "mv config.txt {output}"


rule write_rnaseq_cnv_meta_file:
    input:
        r_script = "scripts/write_rnaseq_cnv_meta.R",
        count_dir = config["count_directory"],
        vcf_dir = "ctat_output_directory/"

    params:
        sample_id = "test006" # TODO: Change to sample id of config

    output:
        "data/meta.txt"

    shell:
        "Rscript {input.r_script} {params.sample_id} {input.count_dir} {input.vcf_dir} {output};"
        "mv meta.txt {output}"


#rule run_rnaseq_cnv:
#    input:
#        r_script = "scripts/run_rnaseq_cnv.R",
#        config_file = "data/config.txt",
#        metadata_file = config["meta"]
#    output:
#        "rnaseq_cnv_output_directory/"

#    shell:
#        "Rscript {input.r_script} {input.config_file} {input.metadata_file}  {output};"
#        "mv predictions.tsv {output}"


