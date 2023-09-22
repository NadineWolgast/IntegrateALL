import os
import csv
import pandas as pd
import shutil

configfile: "config.yaml"


samples = {}
with open(config["sample_file"], "r") as f:
    next(f)
    for line in f:
        sample_id, left, right = line.strip().split(",")
        samples[sample_id] = (left, right)


def extract_filenames_from_csv(csv_file):
    left_filenames = []
    right_filenames = []

    with open(csv_file, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            left_filenames.append(row['left'].split('/')[-1])  # Extract only the filename
            right_filenames.append(row['right'].split('/')[-1])  # Extract only the filename

    return left_filenames, right_filenames

# Usage example:
csv_file = 'samples.csv'  # Replace with the actual path to your samples.csv file
left_files, right_files = extract_filenames_from_csv(csv_file)
print("Left Filenames:", left_files)
print("Right Filenames:", right_files)

sample_ids_list = ['reads', 'testReads', 'reads2']

def create_sample_dataframe(sample_sheet_path):
    sample_df = pd.read_csv(sample_sheet_path)
    result_data = []
    for index, row in sample_df.iterrows():
        sample_id = row['sample_id']
        left_fastq = row['left']
        right_fastq = row['right']
        result_data.append({'sample_id': f"{sample_id}_left", 'FASTQ': left_fastq})
        result_data.append({'sample_id': f"{sample_id}_right", 'FASTQ': right_fastq})
    result_df = pd.DataFrame(result_data)
    return result_df

# Example usage:
sample_sheet_path = config["sample_file"]
fastq_dataframe = create_sample_dataframe(sample_sheet_path)
print(fastq_dataframe)

def get_unique_paths_without_extension(dataframe):
    unique_paths = set()

    for fastq_path in dataframe['FASTQ']:
        # Extract the directory portion of the path and remove the file extension
        path_without_extension = os.path.splitext(os.path.dirname(fastq_path))[0]
        unique_paths.add(path_without_extension)

    return list(unique_paths)

fastq_directory = get_unique_paths_without_extension(fastq_dataframe)
print(fastq_directory)


# Define the samples dictionary with sample IDs and corresponding left/right files
samples_test = {
    "reads": {"left": "/reads_1.fq.gz", "right": "/reads_2.fq.gz"},
    "testReads": {"left": "/test-reads-A01_R1_001.fastq.gz", "right": "/test-reads-A01_R2_001.fastq.gz"},
}



rule all:
    input:
        #expand("STAR_output/{sample_id}/Aligned.sortedByCoord.out.bam",sample_id=list(samples.keys())),
        #expand("data/single_counts/{sample_id}.txt", sample_id=list(samples.keys())),
        #expand("fusions/{sample_id}.tsv",sample_id=list(samples.keys())),
        #expand("fastqc/{sample}", sample=fastq_dataframe['sample_id']),
        #"multiqc/multiqc_report.html",
        #"allcatch_output/predictions.tsv",
        #expand("ctat_output_directory/{sample_id}/{sample_id}.cancer.vcf",sample_id=samples_test.keys()),
        #expand("fusions/{sample_id}.tsv",sample_id=samples_test.keys()),
        #expand("data/single_counts/{sample_id}.txt", sample_id=samples.keys()),
        #expand("data/vcf_files/{sample_id}.tsv", sample_id=list(samples.keys())),
        #"data/config.txt",
        #"data/meta.txt",
        #"rnaseq_cnv_output_directory/",
        #"data/vcf_files/21Ord12062.vcf",
        #"data/vcf_files/21Ord12062.tsv",
        "rnaseq_cnv_output_directory/reads",
        #"data/single_counts"
        #"fusioncatcher_output/",
        #"multiqc/multiqc_report.html"
        #"fusions/reads.tsv",
        #'/media/nadine/INTENSO/STAR/hg38_index'




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
        #expand("fastqc/{Sample}.zip",Sample=inputdf['Sample_ID'].unique())
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

'''
rule multiqc_file:
    input:
        expand("fastqc/{sample}", sample=config["samples"])

    output:
        "multiqc/{sample}/multiqc_report.html"

    params:
        extra = "-o multiqc/{sample}/",
        use_input_files_only=True,

    log:
        "logs/multiqc/multiqc.log"
    wrapper:
        "v2.6.0/bio/multiqc"
'''

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
        fastq1 = lambda wildcards: samples[wildcards.sample_id][0],  # Path to left FASTQ from sample sheet
        fastq2 = lambda wildcards: samples[wildcards.sample_id][1],  # Path to right FASTQ from sample sheet

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



rule run_arriba:
    input:
        # STAR BAM containing chimeric alignments from 'run_star_aligner'
        bam="STAR_output/{sample_id}/Aligned.sortedByCoord.out.bam",
        genome=config["star_ref"],       # Path to reference genome
        annotation=config["star_gtf"],   # Path to annotation GTF
        custom_blacklist=[],
    output:
        # Approved gene fusions
        fusions="fusions/{sample_id}.tsv",
        # Discarded gene fusions (optional)
        discarded="fusions/{sample_id}.discarded.tsv",
    log:
        "logs/arriba/{sample_id}.log",
    params:
        genome_build="GRCh38",           # Required when blacklist or known_fusions is set
        default_blacklist=False,         # Optional
        default_known_fusions=True,      # Optional
        sv_file="",                      # File containing information from structural variant analysis
        extra="-i 1,2",                  # Optional parameters
    threads: 1
    wrapper:
        "v2.6.0/bio/arriba"



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


#TODO: Should this code be stored as an python script and then executed in a rule?
def merge_reads_per_gene_files(input_dir, output_file):
    # Initialize an empty DataFrame to store the merged data
    merged_df = pd.DataFrame()

    # Iterate through subdirectories (each corresponds to a sample)
    for sample_dir in os.listdir(input_dir):
        sample_path = os.path.join(input_dir,sample_dir)

        # Check if the path is a directory
        if os.path.isdir(sample_path):
            sample_name = os.path.basename(sample_path)
            sample_file = os.path.join(sample_path,'ReadsPerGene.out.tab')

            # Check if the ReadsPerGene.out.tab file exists
            if os.path.exists(sample_file):
                # Read the file, skip the first 4 rows, and select the first and fourth columns
                data = pd.read_csv(sample_file,sep='\t',header=None,skiprows=4,usecols=[0, 1],names=['Gene',
                                                                                                     sample_name])

                # Set the gene name column as the index
                data.set_index('Gene',inplace=True)

                # If merged_df is empty, initialize it with the first data
                if merged_df.empty:
                    merged_df = data
                else:
                    # Merge data with the existing merged_df based on the index (gene names)
                    merged_df = pd.merge(merged_df,data,left_index=True,right_index=True,how='outer')

    # Reset the index to have the gene names as a separate column
    merged_df.reset_index(inplace=True)

    # Save the merged DataFrame to the output file
    merged_df.to_csv(output_file,sep='\t',index=False)
input_directory = 'STAR_output'
output_file = 'merged_reads_per_gene.tsv'
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
        input_file = 'merged_reads_per_gene.tsv'
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

    params:
        sample_id= lambda wildcards: wildcards.sample_id,
        #output_dir=config["ctat_output_dir"],
        genome_lib=config["genome_lib"],
        left= lambda wildcards: samples_test[wildcards.sample_id]['left'],
        right= lambda wildcards: samples_test[wildcards.sample_id]['right']

    wildcard_constraints:
        sample_id="|".join(samples_test.keys())

    output:
        vcf ="ctat_output_directory/{sample_id}/{sample_id}.cancer.vcf",
        directory= directory("ctat_output_directory/{sample_id}/"),


    shell:
        '''
        singularity exec -e -B {input.mount_dir}:/data \
        -B {params.genome_lib}:/ctat_genome_lib_dir \
        ctat_mutations.v3.2.0.simg \
        /usr/local/src/ctat-mutations/ctat_mutations \
        --left {input.mount_dir}{params.left} \
        --right {input.mount_dir}{params.right} \
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
        directory("rnaseq_cnv_output_directory/{sample_id}")

    shell:
        "Rscript {input.r_script} {input.config_file} {input.metadata_file} {output};"
        "sleep 0.10;"
        "mv rnaseq_cnv_output_directory/estimation_table.tsv {output};"
        "mv rnaseq_cnv_output_directory/manual_an_table.tsv {output};"
        "mv rnaseq_cnv_output_directory/log2_fold_change_per_arm.tsv {output};"
        "mv rnaseq_cnv_output_directory/alteration_matrix.tsv {output};"


