configfile: "config.yaml"


rule all:
    input:
        #"allcatch_output/predictions.tsv",
        #"data/config.txt",
        #"data/meta.txt",
        #"rnaseq_cnv_output_directory/",
        #"data/vcf_files/21Ord12062.vcf",
        #"data/vcf_files/21Ord12062.tsv",
        #"rnaseq_cnv_output_directory/21Ord12062",
        #"data/single_counts/21Ord12062.txt"
        #'Star_output/'
        '/media/nadine/INTENSO/STAR/hg38_index'


#rule star_pe_multi:
#    input:
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


#rule bam_wta_index:
#    input:
#        "star/Line{index}/WTA_Aligned.sortedByCoord.out.bam"
#    output:
#        "star/Line{index}/WTA_Aligned.sortedByCoord.out.bam.bai"
#    shell:
#        "samtools index {input}"

# TODO: rule index not tested yet
rule index:
        input:
            fa = config['star_ref'], # provide your reference FASTA file
            gtf = config['star_gtf'] # provide your GTF file
        output:
            directory('/media/nadine/INTENSO/STAR/hg38_index') # you can rename the index folder
        threads: 20 # set the maximum number of available cores
        shell:
            'mkdir {output} && '
            'STAR --runThreadN {threads} '
            '--runMode genomeGenerate '
            '--genomeDir {output} '
            '--genomeFastaFiles {input.fa} '
            '--sjdbGTFfile {input.gtf} '
            '--sjdbOverhang 100'


rule run_star_aligner:
    input:
        fastq1= config["STAR_left_samples"],
        fastq2= config["STAR_right_samples"],

    output:
        directory('data/Star_output')

    resources:
        threads=config['threads'],
        mem=config['star_mem']
    shell:
        'mkdir {output} && '
        'STAR --runThreadN {config[threads]} '
        '--runMode alignReads '
        '--genomeDir {config[genome_index]} '
        '--readFilesIn {input.fastq1} {input.fastq2} '
        '--outFileNamePrefix {output}/{input.fastq1} '
        '--quantMode GeneCounts '
        '--sjdbOverhang 100 '
        '--twopassMode Basic '
        '--outSAMtype BAM SortedByCoordinate '
        '--genomeLoad NoSharedMemory '
        '--outFilterMultimapNmax 10 '
        '--outTmpDir {config[star_tmp_directory]}'
        '--chimOutType WithinBAM '
        '--chimSegmentMin 10 '
        '--readFilesCommand zcat'



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
        #TODO: How can I pass here the input left and right samples?

    params:
        sample_id = "test003", # TODO Replace with the actual sample ID
        output_dir = config["ctat_output_dir"],
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
# TODO: Test with boosting method none & ctat_mutations_latest.sif

rule write_rnaseq_cnv_config_file:
    input:
        r_script = "scripts/write_rnaseq_cnv_config.R",
        count_dir = config["single_counts_directory"],
        snv_dir = "data/vcf_files",
        out_dir= "rnaseq_cnv_output_directory"

    output:
        "data/config.txt"

    shell:
        "Rscript {input.r_script} {input.out_dir} {input.count_dir} {input.snv_dir} {output};"
        "mv config.txt {output}"


rule write_rnaseq_cnv_meta_file:
    input:
        r_script = "scripts/write_rnaseq_cnv_meta.R",
        count_dir = config["single_counts_directory"],
        vcf_dir = "data/vcf_files"

    params:
        sample_id = config["sample_id"] # TODO: Change to sample id of config

    output:
        "data/meta.txt"

    shell:
        "Rscript {input.r_script} {params.sample_id} {input.count_dir} {input.vcf_dir} {output};"
        "mv meta.txt {output}"


rule install_rnaseq_cnv:
    shell:
        "Rscript -e 'devtools::install_github(\"honzee/RNAseqCNV\")'"


rule prepare_count_files:
    input:
        reads_dir = config["reads_per_gene_directory"],
        r_script = "scripts/prepare_count_files.R"

#21Ord12062
    params:
        sample_id=config["sample_id"],  # TODO: Change to sample id of config
        count_directory= config["single_counts_directory"]

    output:
        "{count_directory}{sample_id}.txt"

    shell:
        "Rscript {input.r_script} {params.sample_id} {input.reads_dir} {output};"
        "mv {params.sample_id}.txt {output}"



rule prepare_vcf_files:
    input:
        vcf_dir = config["ctat_output_dir"],
        r_script = "scripts/prepare_vcf_files.R"

    params:
        sample_id = config["sample_id"],
        vcf_identifier = "-S1-L1_S69_L003_R1"

    output:
        "data/vcf_files/{sample_id}.tsv"

    shell:
        "Rscript {input.r_script} {params.sample_id} {params.vcf_identifier} {input.vcf_dir} {output};"
        "mv {input.vcf_dir}/{params.sample_id}.tsv {output}"


rule run_rnaseq_cnv:
    input:
        r_script = "scripts/run_rnaseq_cnv.R",
        config_file = "data/config.txt",
        metadata_file = "data/meta.txt"
    params:
        sample_id = config["sample_id"]
    output:
        directory("rnaseq_cnv_output_directory/{sample_id}")

    shell:
        "Rscript {input.r_script} {input.config_file} {input.metadata_file} {output};"
        "sleep 0.10;"
        "mv rnaseq_cnv_output_directory/estimation_table.tsv {output};"
        "mv rnaseq_cnv_output_directory/manual_an_table.tsv {output};"
        "mv rnaseq_cnv_output_directory/log2_fold_change_per_arm.tsv {output};"
        "mv rnaseq_cnv_output_directory/alteration_matrix.tsv {output};"


