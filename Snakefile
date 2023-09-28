import os
import csv
import pandas as pd


configfile: "config.yaml"


# Import your custom Python scripts
from scripts.extract_filenames_from_csv import extract_filenames_from_csv
from scripts.create_sample_dataframe import create_sample_dataframe
from scripts.get_unique_paths_without_extension import get_unique_paths_without_extension
from scripts.merge_reads_per_gene_files import merge_reads_per_gene_files


samples = {}
with open(config["sample_file"], "r") as f:
    next(f)
    for line in f:
        sample_id, left, right = line.strip().split(",")
        samples[sample_id] = (left, right)

# Get data from input sample sheet for the rules:
sample_file = config["sample_file"]
left_files, right_files = extract_filenames_from_csv(sample_file)
fastq_dataframe = create_sample_dataframe(sample_file)
fastq_directory = get_unique_paths_without_extension(fastq_dataframe)


#sample_id,left,right
#testReads,data/samples/test-reads-A01_R1_001.fastq.gz,data/samples/test-reads-A01_R2_001.fastq.gz
#reads,data/samples/reads_1.fq.gz,data/samples/reads_2.fq.gz

df = pd.read_csv(config['sample_file'],sep=',')

# Initialize an empty dictionary to store the results
samples_test = {}

# Iterate over the rows of the DataFrame
for index, row in df.iterrows():
    sample_id = row['sample_id']
    left = os.path.basename(row['left'])
    right = os.path.basename(row['right'])

    # Add a backslash in front of the filenames
    left = '/' + left
    right = '/' + right

    # Create a nested dictionary for each sample
    sample_info = {
        "left": left,
        "right": right
    }

    # Add the sample to the samples_test dictionary
    samples_test[sample_id] = sample_info

# Print the resulting dictionary
print(samples_test)

def extract_sample_ids_from_meta(meta_file):
    df = pd.read_csv(meta_file, sep='\t', header=None)  # Avoid header inference

    # Extract and return the values from the first column as a list
    sample_ids = df.iloc[:, 0].tolist()

    return sample_ids


rnaseqcnv_sample_ids = extract_sample_ids_from_meta("data/meta.txt")


rule all:
    input:
        #expand("STAR_output/{sample_id}/Aligned.sortedByCoord.out.bam",sample_id=list(samples.keys())),
        #expand("data/single_counts/{sample_id}.txt", sample_id=list(samples.keys())),
        #expand("fusions/{sample_id}.tsv",sample_id=list(samples.keys())),
        #expand("fastqc/{sample}", sample=fastq_dataframe['sample_id']),
        expand("multiqc/{sample}/multiqc_report.html",sample=fastq_dataframe['sample_id']),
        "multiqc/multiqc_report.html",
        #"allcatch_output/predictions.tsv",
        #expand("ctat_output_directory/{sample_id}/{sample_id}.cancer.vcf",sample_id=samples_test.keys()), # Funktioniert
        #expand("ctat_output_directory/{sample_id}/{sample_id}.cancer.vcf",sample_id=samples_test.keys()),
        #expand("fusions/{sample_id}.tsv",sample_id=samples_test.keys()),
        #expand("data/single_counts/{sample_id}.txt", sample_id=samples.keys()),
        #expand("data/vcf_files/{sample_id}.tsv", sample_id=list(samples.keys())),
        #"data/config.txt",
        #"data/meta.txt",
        #expand("RNAseqCNV_output/{sample_id}", sample_id=rnaseqcnv_sample_ids),
        #"data/single_counts",
        #"fusioncatcher_output/",
        #expand("fusioncatcher_output/{sample_id}",sample_id=rnaseqcnv_sample_ids),
        #'/media/nadine/HOME/nadine/STAR/ensembl_94_100'
        #expand("STAR_output/{sample}/Aligned.sortedByCoord.out.bam.bai",sample=list(samples.keys())),
        #expand("pysamstats_output_dir/{sample_id}.coverage.txt",sample_id=list(samples.keys()))



# Function to get input FASTQ files based on sample_id
def get_input_fastqs(wildcards):
    sample_id = wildcards.sample
    fastq_file = fastq_dataframe[fastq_dataframe['sample_id'] == sample_id]['FASTQ'].values[0]
    return fastq_file


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
        sample="fastqc/{sample}.zip"
    output:
        directory("fastqc/{sample}")
    shell:
        "unzip -q {input.sample} -d {output}"


rule multiqc_dir:
    input:
        "fastqc/"
    output:
        "multiqc/multiqc_report.html"
    params:
        output_dir= "multiqc/"
    shell:
        "multiqc {input} -o {params.output_dir}"


rule multiqc_file:
    input:
        "fastqc/{sample}/"
    output:
        "multiqc/{sample}/multiqc_report.html"
    params:
        output_dir= "multiqc/{sample}/"
    shell:
        "multiqc {input} -o {params.output_dir}"



#conda install star=2.7.1a

# TODO: rule index not tested yet
rule index:
        input:
            fa = config['star_ref'], # provide your reference FASTA file
            gtf = config['star_gtf'] # provide your GTF file
        output:
            directory(config["genome_index"]) # TODO: Change to config for genome index
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
        fastq1 = lambda wildcards: samples[wildcards.sample_id][0],  # Path to left FASTQ from sample sheet
        fastq2 = lambda wildcards: samples[wildcards.sample_id][1]  # Path to right FASTQ from sample sheet

    output:
        directory = directory("STAR_output/{sample_id}"),
        bam = "STAR_output/{sample_id}/Aligned.sortedByCoord.out.bam",
        reads = "STAR_output/{sample_id}/ReadsPerGene.out.tab"

    resources:
        threads=config['threads'],
        mem=config['star_mem']
    shell:
        'mkdir {output.directory} && '
        'sleep .10;'
        'STAR --runThreadN {config[threads]} '
        '--runMode alignReads '
        '--genomeDir {config[genome_index]} '
        '--readFilesIn {input.fastq1} {input.fastq2} '
        '--outFileNamePrefix  {output.directory}/ '
        '--quantMode GeneCounts '
        '--sjdbOverhang 100 '
        '--twopassMode Basic '
        '--outSAMtype BAM SortedByCoordinate '
        '--genomeLoad NoSharedMemory '
        '--outFilterMultimapNmax 10 '
        '--chimOutType WithinBAM '
        '--chimSegmentMin 10 '
        '--readFilesCommand zcat '


rule samtools_index:
    input:
        "STAR_output/{sample}/Aligned.sortedByCoord.out.bam"
    output:
        "STAR_output/{sample}/Aligned.sortedByCoord.out.bam.bai",
    log:
        "logs/samtools_index/{sample}.log",
    params:
        extra="",  # optional params string
    threads: 4  # This value - 1 will be sent to -@
    wrapper:
        "v2.6.0/bio/samtools/index"


rule pysamstat:
    input:
        bam="STAR_output/{sample_id}/Aligned.sortedByCoord.out.bam"

    output:
        "pysamstats_output_dir/{sample_id}.coverage.txt"

    shell:
        "pysamstats --type coverage {input.bam} > pysamstats_output_dir/{sample_id}.coverage.txt"


rule run_arriba:
    input:
        # STAR BAM containing chimeric alignments from 'run_star_aligner'
        bam="STAR_output/{sample_id}/Aligned.sortedByCoord.out.bam",
        genome=config["star_ref"],       # Path to reference genome
        annotation=config["star_gtf"],   # Path to annotation GTF
        custom_blacklist=[]
    output:
        # Approved gene fusions
        fusions="fusions/{sample_id}.tsv",
        # Discarded gene fusions (optional)
        discarded="fusions/{sample_id}.discarded.tsv"
    log:
        "logs/arriba/{sample_id}.log"
    params:
        genome_build="GRCh38",           # Required when blacklist or known_fusions is set
        default_blacklist=False,         # Optional
        default_known_fusions=True,      # Optional
        sv_file="",                      # File containing information from structural variant analysis
        extra="-i 1,2"                  # Optional parameters
    threads: 1
    wrapper:
        "v2.6.0/bio/arriba"





rule run_fusioncatcher:
    input:
        fastq_directory = get_unique_paths_without_extension(fastq_dataframe),
        data_directory = "/media/nadine/HOME/nadine/fusioncatcher/data/human_v102",
        mount_dir= config["ctat_mount_dir"]

    output:
        directory("fusioncatcher_output/{sample_id}")

    conda:
        "envs/fusioncatcher.yaml"

    params:
        sample_id=lambda wildcards: wildcards.sample_id  # Extract sample_id from wildcards

    shell:
        '''mkdir {output} &&
        fusioncatcher \
        -d {input.data_directory} \
        -i {input.fastq_directory} \
        -o {output}'''


"""
rule run_fusioncatcher:
    input:
        fastq_directory = get_unique_paths_without_extension(fastq_dataframe),
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

"""

input_directory = 'STAR_output'
output_file = 'data/combined_counts/merged_reads_per_gene.tsv'
merge_reads_per_gene_files(input_directory,output_file)


rule install_allcatchr:
    shell:
        "Rscript -e 'devtools::install_github(\"ThomasBeder/ALLCatchR\")'"

# TODO: Check if rule should run for each individual ReadsPerGene.out.tab file of STAR or only for aggregated file
# TODO: Run ALLCatchR for each individual ReadsPerGene.out.tab file, merge final outputs in one table.
#  Need to rename the output files predictions.tsv from running the command for each individual file.
# Rule to run ALLCatchR
rule run_allcatchr:
    input:
        r_script = "scripts/run_ALLCatchR.R",
        #input_file = config["counts"],  # Update with the correct path
        input_file = 'data/combined_counts/merged_reads_per_gene.tsv'
    output:
        "allcatch_output/predictions.tsv"

    shell:
        "Rscript {input.r_script} {input.input_file} {output};"
        "mv predictions.tsv {output}"



rule pull_ctat_mutations_singularity_image:
    shell:
        "singularity pull docker://trinityctat/ctat_mutations"


# Rule to run ctat-mutations
# TODO: Test with alternative mount directory and sample directories
# TODO: Test with boosting method none & ctat_mutations_latest.sif
rule run_ctat_mutations:
    input:
        mount_dir=config["ctat_mount_dir"],
        fastq1= lambda wildcards: samples[wildcards.sample_id][0],  # Path to left FASTQ from sample shee
        fastq2= lambda wildcards: samples[wildcards.sample_id][1],  # Path to right FASTQ from sample sheet
        input_directory= config["ctat_input_directory"]

    params:
        sample_id= lambda wildcards: wildcards.sample_id,
        genome_lib=config["genome_lib"],
        left= lambda wildcards: samples_test[wildcards.sample_id]['left'],
        right= lambda wildcards: samples_test[wildcards.sample_id]['right'],
        #input_directory= config["ctat_input_directory"]

    wildcard_constraints:
        sample_id="|".join(samples_test.keys())

    output:
        vcf ="ctat_output_directory/{sample_id}/{sample_id}.cancer.vcf",
        directory= directory("ctat_output_directory/{sample_id}/")


    shell:
        '''
        singularity exec -e -B {input.mount_dir}:/data \
        -B {params.genome_lib}:/ctat_genome_lib_dir \
        -B {input.input_directory}:/ctat_input \
        ctat_mutations.v3.2.0.simg \
        /usr/local/src/ctat-mutations/ctat_mutations \
        --left /ctat_input/{params.left} \
        --right /ctat_input/{params.right} \
        --sample_id {params.sample_id} \
        --output {output.directory}  \
        --cpu 10 \
        --genome_lib_dir /ctat_genome_lib_dir \
        --boosting_method none \
        --no_cravat
        '''

# Rule to process ReadsPerGene.out.tab files for RNASeq-CNV
rule process_reads_per_gene:
    input:
        reads_per_gene="STAR_output/{sample_id}/ReadsPerGene.out.tab"  # Input pattern for each sample
    output:
        counts="data/single_counts/{sample_id}.txt"  # Output file for each sample
    params:
        skip_rows=4  # Number of rows to skip in the input file
    shell:
        """
        awk 'NR > {params.skip_rows} {{print $1 "\t" $2}}' {input.reads_per_gene} > {output.counts}
        """


rule write_rnaseq_cnv_files:
    input:
        r_script_config = "scripts/write_rnaseq_cnv_config.R",
        r_script_meta = "scripts/write_rnaseq_cnv_meta.R",
        count_dir = "data/single_counts/",
        vcf_dir = "data/vcf_files/",
        out_dir= "rnaseq_cnv_output_directory/"


    output:
        config="data/config.txt",
        meta="data/meta.txt"

    log:
        "logs/write_config_file/meta.log"

    shell:
        """
        # Generate config.txt
        Rscript {input.r_script_config} {input.out_dir} {input.count_dir} {input.vcf_dir} {output.config} > {log};

        # Generate meta.txt
        Rscript {input.r_script_meta}  {input.count_dir} {input.vcf_dir} {output.meta} > {log};
        """

rule install_rnaseq_cnv:
    shell:
        "Rscript -e 'devtools::install_github(\"honzee/RNAseqCNV\")'"

# TODO: Check if all the cancer.vcf files need to skip 263 lines. Needs maybe a modified skript
rule prepare_vcf_files:
    input:
        vcf_dir="ctat_output_directory/{sample_id}.cancer.vcf",
        r_script= "scripts/prepare_vcf_files.R"
    output:
        counts="data/vcf_files/{sample_id}.tsv"

    shell:
        # Execute the R script using Rscript command
        "Rscript {input.r_script} {input.vcf_dir} {output.counts};"




rule run_rnaseq_cnv:
    input:
        r_script = "scripts/run_rnaseq_cnv.R",
        config_file = "data/config.txt",
        metadata_file = "data/meta.txt"
    params:
        sample_id = config["sample_id"]
    output:
        directory("RNAseqCNV_output/{sample_id}/")

    shell:
        "Rscript {input.r_script} {input.config_file} {input.metadata_file} {output};"
        "sleep 0.10;"
        "mkdir {output};"
        "mv estimation_table.tsv {output};"
        "mv manual_an_table.tsv {output};"
        "mv log2_fold_change_per_arm.tsv {output};"
        "mv alteration_matrix.tsv {output};"
        "mv {params.sample_id}/{params.sample_id}_CNV_main_fig.png  {output}"


