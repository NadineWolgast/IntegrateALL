import glob
import os
import pandas as pd


configfile: "config.yaml"

# TODO: snakemake --list-conda-envs to get all installed packages

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
print("fastq dataframe: ", fastq_dataframe)

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
print("samples_test:",samples_test)

def extract_sample_ids_from_meta(meta_file):
    df = pd.read_csv(meta_file, sep='\t', header=None)  # Avoid header inference

    # Extract and return the values from the first column as a list
    sample_ids = df.iloc[:, 0].tolist()

    return sample_ids


rnaseqcnv_sample_ids = extract_sample_ids_from_meta("data/meta.txt")
print("sample_keys: ", list(samples.keys()))
print("rnaseqcnv sample_keys: ", list(rnaseqcnv_sample_ids))


rule all:
    input:
        #expand("STAR_output/{sample_id}/Aligned.sortedByCoord.out.bam",sample_id=list(samples.keys())),
        #expand("data/single_counts/{sample_id}.txt", sample_id=list(samples.keys())),
        #expand("fusions/{sample_id}.tsv",sample_id=list(samples.keys())),
        #expand("fastqc/{sample}", sample=fastq_dataframe['sample_id']),
        #expand("multiqc/{sample}/multiqc_report.html",sample=fastq_dataframe['sample_id']),
        #expand("multiqc/{sample}/multiqc_data/multiqc_fastqc.txt", sample=fastq_dataframe['sample_id']),
        #expand("data/counts/{sample_id}.tsv",sample_id= samples.keys()),
        #expand("allcatch_output/{sample_id}/predictions.tsv",sample_id= samples.keys()),
        expand("ctat_output_directory/{sample_id}/{sample_id}.cancer.vcf",sample_id=samples_test.keys()), # Funktioniert
        #expand("fusions/{sample_id}.tsv",sample_id=samples_test.keys()),
        #expand("data/single_counts/{sample_id}.txt", sample_id=samples.keys()),
        #expand("data/vcf_files/{sample_id}.tsv", sample_id=list(samples.keys())),
        #"data/config.txt",
        #"data/meta.txt",
        #expand("/media/nadine/Alina/Blast-o-Matic-Fusioninator/RNAseqCNV_output/{sample_id}", sample_id=rnaseqcnv_sample_ids),
        #"data/single_counts",
        #"fusioncatcher_output/",
        #expand("fusioncatcher_output/{sample_id}",sample_id= samples.keys()),
        #expand("data/total_mapped_reads/{sample_id}.txt", sample_id= samples.keys()), # doesn't work yet
        #expand("data/tpm/{sample_id}.tsv",sample_id= samples.keys())  # doesn't work yet
        #'/media/nadine/HOME/nadine/STAR/ensembl_94_100'
        #expand("STAR_output/{sample}/Aligned.sortedByCoord.out.bam.bai",sample=list(samples.keys())),
        #expand("pysamstats_output_dir/{sample_id}.coverage.txt",sample_id=list(samples.keys())),
        #expand("aggregated_output/{sample}.csv", sample=list(samples.keys()))



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
        multiqc_report="multiqc/{sample}/multiqc_report.html",
        multiqc_fqc= "multiqc/{sample}/multiqc_data/multiqc_fastqc.txt",

    params:
        output_dir= "multiqc/{sample}/",

    shell:
        "multiqc {input} -o {params.output_dir}"



# TODO: need to add file prefix left and right or it overwrites itself...
'''
rule extract_and_rename_QC_files:
    input:
        multiqc_dir="multiqc/",
        fastqc_dir="fastqc/",

    output:
        directory("QC/{sample}"),
        multiqc_fastqc="QC/{sample}/multiqc_fastqc.txt",
        multiqc_general_stats="QC/{sample}/multiqc_general_stats.txt",
        fastqc_data="QC/{sample}/fastqc_data.txt",
        summary="QC/{sample}/summary.txt"


    shell:
        """
            mkdir -p QC/{wildcards.sample}
            echo Processing wildcards sample: {wildcards.sample}
            # Copy the MultiQC files\n"
            cp {input.multiqc_dir}/{wildcards.sample}_left/multiqc_data/multiqc_fastqc.txt {output.multiqc_fastqc}
            cp {input.multiqc_dir}/{wildcards.sample}_right/multiqc_data/multiqc_fastqc.txt {output.multiqc_fastqc}
            cp {input.multiqc_dir}/{wildcards.sample}_left/multiqc_data/multiqc_general_stats.txt {output.multiqc_general_stats}
            cp {input.multiqc_dir}/{wildcards.sample}_right/multiqc_data/multiqc_general_stats.txt {output.multiqc_general_stats}
            find {input.fastqc_dir}/{wildcards.sample}_left/ -type f -name 'fastqc_data.txt' -exec cp {} QC/{wildcards.sample}/ \ 
            find {input.fastqc_dir}/{wildcards.sample}_right/ -type f -name 'fastqc_data.txt' -exec cp {} QC/{wildcards.sample}/ \; 
        """
#find fastqc/reads_left/ -type f -name 'fastqc_data.txt' -exec cp {} QC/reads/ \;
'''

#conda install star=2.7.1a

# TODO: rule index not tested yet, too large for my machine
rule index:
        input:
            fa = config['star_ref'], # provide your reference FASTA file
            gtf = config['star_gtf'] # provide your GTF file
        output:
            directory(config["genome_index"])
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
        reads = "STAR_output/{sample_id}/ReadsPerGene.out.tab",
        log_out="STAR_output/{sample_id}/Log.final.out"

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
        "STAR_output/{sample}/Aligned.sortedByCoord.out.bam.bai"
    log:
        "logs/samtools_index/{sample}.log"
    params:
        extra=""  # optional params string
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

# TODO:
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
    threads: 10
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



input_directory = 'STAR_output'
output_file = 'data/combined_counts/merged_reads_per_gene.tsv'
merge_reads_per_gene_files(input_directory,output_file)


rule install_allcatchr:
    shell:
        "Rscript -e 'devtools::install_github(\"ThomasBeder/ALLCatchR\")'"


# Rule to run ALLCatchR
rule run_allcatchr:
    input:
        r_script = "scripts/run_ALLCatchR.R",
        input_file = 'data/combined_counts/merged_reads_per_gene.tsv'
    output:
        "allcatch_output/predictions.tsv"

    shell:
        "Rscript {input.r_script} {input.input_file} {output};"
        "mv predictions.tsv {output}"



# Rule to process ReadsPerGene.out.tab files for ALLCatchR
rule process_reads_per_gene_to_counts:
    input:
        reads_per_gene="STAR_output/{sample_id}/ReadsPerGene.out.tab"  # Input pattern for each sample
    output:
        counts="data/counts/{sample_id}.tsv"  # Output file for each sample
    params:
        skip_rows=4  # Number of rows to skip in the input file
    shell:
        """
        echo -e "Gene\t{wildcards.sample_id}" > {output.counts}
        awk 'NR > {params.skip_rows} {{print $1 "\t" $2}}' {input.reads_per_gene} >> {output.counts}
        """



rule run_allcatchr_on_single_count_files:
    input:
        r_script = "scripts/run_ALLCatchR.R",
        input_file = 'data/counts/{sample}.tsv'
    output:
        "allcatch_output/{sample}/predictions.tsv"

    shell:
        "Rscript {input.r_script} {input.input_file} {output};"
        "mv predictions.tsv {output}"



rule pull_ctat_mutations_singularity_image:
    shell:
        "singularity pull docker://trinityctat/ctat_mutations"


# Rule to run ctat-mutations
rule run_ctat_mutations:
    input:
        input_directory= config["ctat_input_directory"]


    params:
        sample_id= lambda wildcards: wildcards.sample_id,
        genome_lib=config["genome_lib"],
        left= lambda wildcards: samples_test[wildcards.sample_id]['left'],
        right= lambda wildcards: samples_test[wildcards.sample_id]['right'],
        input_directory= config["ctat_input_directory"],
        output_dir= config["ctat_output_dir"]

    wildcard_constraints:
        sample_id="|".join(samples_test.keys())

    output:
        vcf ="ctat_output_directory/{sample_id}/{sample_id}.cancer.vcf",
        directory= directory("ctat_output_directory/{sample_id}/")


    shell:
        '''
        singularity exec -e -B {input.input_directory}:/ctat_input \
        -B {params.genome_lib}:/ctat_genome_lib_dir:ro \
        -B {params.output_dir}/{output.directory}:/outdir \
        ctat_mutations.v3.2.0.simg \
        /usr/local/src/ctat-mutations/ctat_mutations \
        --left /ctat_input/{params.left} \
        --right /ctat_input/{params.right} \
        --sample_id {params.sample_id} \
        --output /outdir \
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


'''
rule calculate_total_mapped_reads:
    input:
        bam="STAR_output/{sample_id}/Aligned.sortedByCoord.out.bam"
    output:
        total_mapped_reads="data/total_mapped_reads/{sample_id}.txt"
    script:
        "scripts/calculate_total_mapped_reads.sh"



rule calculate_tpm:
    input:
        reads_per_gene="STAR_output/{sample_id}/ReadsPerGene.out.tab",
        total_mapped_reads="data/total_mapped_reads/{sample_id}.txt"
    output:
        tpm="data/tpm/{sample_id}.tsv"
    params:
        skip_rows=4
    script:
        "scripts/calculate_tpm.py"
'''
#TODO: Change count_dir, vcf_dir and outdir to current path
rule write_rnaseq_cnv_files:
    input:
        r_script_config = "scripts/write_rnaseq_cnv_config.R",
        r_script_meta = "scripts/write_rnaseq_cnv_meta.R",
        count_dir = "/media/nadine/Alina/Blast-o-Matic-Fusioninator/data/single_counts/",
        vcf_dir = "/media/nadine/Alina/Blast-o-Matic-Fusioninator/data/vcf_files/",
        out_dir= "/media/nadine/Alina/Blast-o-Matic-Fusioninator/RNAseqCNV_output"


    output:
        config="/media/nadine/Alina/Blast-o-Matic-Fusioninator/data/config.txt",
        meta="/media/nadine/Alina/Blast-o-Matic-Fusioninator/data/meta.txt"

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
        vcf_dir="ctat_output_directory/{sample_id}/{sample_id}.cancer.vcf",
        r_script= "scripts/prepare_vcf_files.R"
    output:
        counts="data/vcf_files/{sample_id}.tsv"

    shell:
        # Execute the R script using Rscript command
        "Rscript {input.r_script} {input.vcf_dir} {output.counts};"




rule run_rnaseq_cnv:
    input:
        r_script = "scripts/run_rnaseq_cnv.R",
        config_file = "/media/nadine/Alina/Blast-o-Matic-Fusioninator/data/config.txt",
        metadata_file = "/media/nadine/Alina/Blast-o-Matic-Fusioninator/data/meta.txt"

    output:
        directory("/media/nadine/Alina/Blast-o-Matic-Fusioninator/RNAseqCNV_output/{sample_id}/")

    shell:
        "Rscript {input.r_script} {input.config_file} {input.metadata_file} {output};"
        "sleep 0.10;"
        "mkdir {output};"
        "mv estimation_table.tsv {output};"
        "mv manual_an_table.tsv {output};"
        "mv log2_fold_change_per_arm.tsv {output};"
        "mv alteration_matrix.tsv {output};"
        "mv {sample_id}/{sample_id}_CNV_main_fig.png  {output}"

#sample_id,left,right
#reads,data/samples/reads_1.fq.gz,data/samples/reads_2.fq.gz
#testReads,data/samples/test-reads-A01_R1_001.fastq.gz,data/samples/test-reads-A01_R2_001.fastq.gz

rule aggregate_output:
    input:
        prediction_file = "/media/nadine/Alina/Blast-o-Matic-Fusioninator/allcatch_output/{sample}/predictions.tsv",
        fusioncatcher_file = "/media/nadine/Alina/Blast-o-Matic-Fusioninator/fusioncatcher_output/{sample}/final-list_candidate-fusion-genes.txt",
        arriba_file = "/media/nadine/Alina/Blast-o-Matic-Fusioninator/fusions/{sample}.tsv",
        rna_seq_cnv_alteration_file= "/media/nadine/Alina/Blast-o-Matic-Fusioninator/RNAseqCNV_output/{sample}/estimation_table.tsv",
        rna_seq_cnv_log2foldchange_file= "/media/nadine/Alina/Blast-o-Matic-Fusioninator/RNAseqCNV_output/{sample}/log2_fold_change_per_arm.tsv",
        rna_seq_cnv_manual_an_table_file= "/media/nadine/Alina/Blast-o-Matic-Fusioninator/RNAseqCNV_output/{sample}/manual_an_table.tsv",
        star_log_final_out_file= "/media/nadine/Alina/Blast-o-Matic-Fusioninator/STAR_output/{sample}/Log.final.out",
        multiqc_fqc_right= "/media/nadine/Alina/Blast-o-Matic-Fusioninator/multiqc/{sample}_right/multiqc_data/multiqc_fastqc.txt",
        multiqc_fqc_left= "/media/nadine/Alina/Blast-o-Matic-Fusioninator/multiqc/{sample}_left/multiqc_data/multiqc_fastqc.txt",


    output:
        csv= "aggregated_output/{sample}.csv"


    shell:
        """
            uniquely_mapped_reads=$(awk -F'\t' '/Uniquely mapped reads number/ {{print $2}}' {input.star_log_final_out_file})
            echo "The transcriptome sequencing of {wildcards.sample} produced $uniquely_mapped_reads uniquely aligned sequencing reads, enabling quantification of protein coding genes." > {output.csv}
            echo -e "Quality metrics (fastQC / MultiQC) indicated:" >> {output.csv}
            echo -e "Filename\tTotal Sequences\tSequences flagged as poor quality\tSequence length\t%GC\ttotal_deduplicated_percentage\tavg_sequence_length\tmedian_sequence_length\tbasic_statistics\tper_base_sequence_quality\tper_sequence_quality_scores\tper_base_sequence_content\tper_sequence_gc_content\tper_base_n_content\tsequence_length_distribution\tsequence_duplication_levels\toverrepresented_sequences\tadapter_content" >> {output.csv}
            awk -F'\t' '{{if (NR == 1) next; printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n", $2, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20}}' {input.multiqc_fqc_left} >> {output.csv}
            awk -F'\t' '{{if (NR == 1) next; printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n", $2, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20}}' {input.multiqc_fqc_right} >> {output.csv}
            echo -e "ALLCatchR allocated for sample {wildcards.sample} the following molecular subtype:" >> {output.csv}
            echo -e "subtype prediction\tscore\tconfidence" >> {output.csv}
            awk -F'\t' '{{if (NR == 1) next; printf "%s\\t%s\\t%s\\n", $3, $2, $4}}' {input.prediction_file} >> {output.csv}
            echo -e "fusioncatcher / ARRIBA identified the following driver fusion candidates:" >> {output.csv}
            echo -e "Fusioncatcher:" >> {output.csv}
            echo -e "5’ gene name\t5’ chr.position\t3’ gene name\t3’chr. position\tfusion unique spanning reads"  >> {output.csv}
            awk -F'\t' '{{if (NR == 1) next; printf "%s\\t%s\\t%s\\t%s\\t%s\\n", $1, $9, $2, $10, $6}}' {input.fusioncatcher_file} >> {output.csv}
            echo -e "ARRIBA:" >> {output.csv}
            echo -e "5’ gene name\t5’ chr.position\t3’ gene name\t3’chr. position\tconfidence"  >> {output.csv}
            awk -F'\t' '{{if (NR == 1) next; printf "%s\\t%s\\t%s\\t%s\\t%s\\n", $1, $5, $2, $6, $13}}' {input.arriba_file} >> {output.csv}
            echo -e "RNASeqCNV identified the following karyotype:" >> {output.csv}
            echo -e "gender	chrom_n	alterations	status	comments"  >> {output.csv}
            awk -F'\t' '{{if (NR == 1) next; printf "%s\\t%s\\t%s\\t%s\\t%s\\n", $2, $3, $4, $5, $6}}' {input.rna_seq_cnv_manual_an_table_file} >> {output.csv}
            echo -e "Chromosome arm calls" >> {output.csv}
            echo -e "chr\tarm\tlog2FC per arm"  >> {output.csv}
            awk -F'\t' '{{if (NR == 1) next; printf "%s\\t%s\\t%s\\n", $1, $2, $3}}' {input.rna_seq_cnv_log2foldchange_file} >> {output.csv}
        """




'''

The transcriptome sequencing of (sample ID) produced xxx (STAR) uniquely aligend sequencing reads(aus Log.final.out), enabling quantification XXX (Alinas script) protein coding genes. 

Quality metrics (fastQC / MultiQC) indicated … (in progress -> Ampel-System)

ALLCatchR allocated this sample to the following molecular subtype:

subtype prediction		score		confidence
Xxx				xxx		xxx

fusioncatcher / ARRIBA identified the following driver fusion candidates

fusioncatchr
5’ gene name	5’ chr.position	3’ gene name	3’chr. position	fusion	unique spanning reads

ARRIBA
5’ gene name	5’ chr.position	3’ gene name	3’chr. position	fusion	unique spanning reads

 alle Fusionen eine BCP-ALL whitelist

Add on: Abfrage auf welches Exon im jeweiligen MANE transkript passt die Bruchpunktkoordinate mit einer Toleranz von 10 bp

RNASeqCNV identified the following karyotype

Graphical output with boxplots

Table with chromosome arm calls

Number of choromosomes in total

The gene expression-based subtype definition of xxx, the fusion calls of xxx and a karyotype profile of xxx fit the diagnosis yyyy / do not fit a single definition. 


Inhaltlich offene Punkte:
Integration von karytype und gene expressionsprofil -> machine learning classifier
Definition einer Mutations-Positiv-Liste: CITAT vs. Pysamstat

'''