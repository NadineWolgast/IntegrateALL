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
        #'STAR_output/'
        "fusioncatcher_output/"
        #'fusioncatcher_output/'
        #"multiqc/multiqc_report_old.html",
        #"fusions/reads.tsv",
        #'/media/nadine/INTENSO/STAR/hg38_index'
        'fastqc/'

def get_input_fastqs(wildcards):
    return config["samples"][wildcards.sample]

rule fastqc:
    input:
        get_input_fastqs
    output:
        html="fastqc/{sample}.html",
        zip="fastqc/{sample}.zip"
    log:
        "logs/fastqc/{sample}/fastqc.log"
    wrapper:
        "0.31.1/bio/fastqc"

rule unzip:
    input:
        "fastqc/{sample}.zip"
    output:
        "{sample}"
    shell:
        "tar -xf {input}"

rule multiqc_dir:
    input:
        "fastqc/"
        #expand("fastqc/{sample}.html", sample=["I29799-L1_S33_L001_R1_001", "I29799-L1_S33_L001_R2_001"])
    output:
        "multiqc/multiqc_report.html"
    log:
        "logs/multiqc/multiqc.log",
    wrapper:
        "v2.6.0/bio/multiqc"


rule multiqc_file:
    input:
        expand("fastqc/{sample}", sample=["I29799-L1_S33_L001_R1_001"])
    output:
        "multiqc/multiqc_I29799.html"
    params:
        extra="",  # Optional: extra parameters for multiqc.
        use_input_files_only=True, # Optional, use only a .txt and don't search folder for files
    log:
        "logs/multiqc/multiqc.log"
    wrapper:
        "v2.6.0/bio/multiqc"


#conda install star=2.7.1a

# TODO: rule index not tested yet
rule index:
        input:
            fa = config['star_ref'], # provide your reference FASTA file
            gtf = config['star_gtf'] # provide your GTF file
        output:
            directory('/media/nadine/HOME/nadine/STAR/ensembl_94_100') # TODO: Change to config for genome index
        threads: 20 # set the maximum number of available cores
        shell:
            'STAR --runThreadN {threads} '
            '--runMode genomeGenerate '
            '--genomeDir {output} '
            '--genomeFastaFiles {input.fa} '
            '--sjdbGTFfile {input.gtf} '
            '--sjdbOverhang 100'


rule run_star_aligner:
    input:
        fastq1= config["left_samples"],
        fastq2= config["right_samples"]

    output:
        directory('STAR_output/')

    resources:
        threads=config['threads'],
        mem=config['star_mem']
    shell:
        'mkdir {output} && '
        'cd {output} &&'
        'STAR --runThreadN {config[threads]} '
        '--runMode alignReads '
        '--genomeDir {config[genome_index]} '
        '--readFilesIn ../{input.fastq1} ../{input.fastq2} '
        '--outFileNamePrefix {output} '
        '--quantMode GeneCounts '
        '--sjdbOverhang 100 '
        '--twopassMode Basic '
        '--outSAMtype BAM SortedByCoordinate '
        '--genomeLoad NoSharedMemory '
        '--outFilterMultimapNmax 10 '
        #'--outTmpDir {config[star_tmp_directory]}'
        '--chimOutType WithinBAM '
        '--chimSegmentMin 10 '
        '--readFilesCommand zcat'


rule arriba:
    input:
        # STAR bam containing chimeric alignments
        bam="STAR_output/STAR_outputAligned.sortedByCoord.out.bam",
        #bam = "/media/nadine/HOME/nadine/BAM/I29799-L1Aligned.sortedByCoord.out.bam",
        # path to reference genome
        genome=config["star_ref"],
        # path to annotation gtf
        annotation=config["star_gtf"],
        # optional arriba blacklist file
        custom_blacklist=[],
    output:
        # approved gene fusions
        #fusions="fusions/{sample}.tsv",
        fusions="fusions/reads.tsv",
        # discarded gene fusions
        #discarded="fusions/{sample}.discarded.tsv",  # optional
        discarded="fusions/reads.discarded.tsv",  # optional
    log:
        #"logs/arriba/{sample}.log",
        "logs/arriba/reads.log",
    params:
        # required when blacklist or known_fusions is set
        genome_build="GRCh38",
        # strongly recommended, see https://arriba.readthedocs.io/en/latest/input-files/#blacklist
        # only set blacklist input-file or blacklist-param
        default_blacklist=False,  # optional
        default_known_fusions=True,  # optional
        # file containing information from structural variant analysis
        sv_file="",
        # optional parameters
        extra="-i 1,2",
    threads: 1
    wrapper:
        "v2.6.0/bio/arriba"


rule run_fusioncatcher:
    input:
        fastq_directory = "data/samples/",
        data_directory = "/media/nadine/HOME/nadine/fusioncatcher/data/human_v102",
        mount_dir= config["ctat_mount_dir"]

    output:
        directory("fusioncatcher_output/")

    #singularity: "fusioncatcher-1.33.sif"

    conda:
        "envs/fusioncatcher.yaml"

    shell:
        '''mkdir {output} &&
        fusioncatcher \
        -d {input.data_directory} \
        -i {input.fastq_directory} \
        -o {output}'''


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

    log:
        "logs/write_config_file/meta.log"

    shell:
        "Rscript {input.r_script} {input.out_dir} {input.count_dir} {input.snv_dir} {output} > {log};"
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
    log:
        "logs/write_meta_file/meta.log"

    shell:
        "Rscript {input.r_script} {params.sample_id} {input.count_dir} {input.vcf_dir} {output} > {log};"
        "mv meta.txt {output}"


rule install_rnaseq_cnv:
    shell:
        "Rscript -e 'devtools::install_github(\"honzee/RNAseqCNV\")'"


rule prepare_count_files:
    input:
        reads_dir = config["reads_per_gene_directory"],
        r_script = "scripts/prepare_count_files.R"

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


