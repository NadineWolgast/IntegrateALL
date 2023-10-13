import os
import csv
import pandas as pd


configfile: "config.yaml"

# Import custom Python scripts
#from scripts.extract_filenames_from_csv import extract_filenames_from_csv
from scripts.create_sample_dataframe import create_sample_dataframe
from scripts.get_unique_paths_without_extension import get_unique_paths_without_extension
from scripts.merge_reads_per_gene_files import merge_reads_per_gene_files


def validate_input(samples_csv, left_col, right_col):
    import os
    import csv

    # Initialize a list to store error messages
    errors = []

    # Check if the samples.csv file exists
    if not os.path.exists(samples_csv):
        errors.append(f"{samples_csv} does not exist")

    # Check if the samples.csv file has the correct format
    try:
        with open(samples_csv, "r") as csv_file:
            csv_reader = csv.reader(csv_file)
            header = next(csv_reader)  # Read the header row

            if len(header) != 3 or header[0] != "sample_id" or header[1] != "left" or header[2] != "right":
                errors.append(f"Invalid header in {samples_csv}. Expected: 'sample_id, left, right'")

            for row in csv_reader:
                sample, left, right = row
                if not os.path.exists(left):
                    errors.append(f"File referenced in 'left' column does not exist: {left}")
                if not os.path.exists(right):
                    errors.append(f"File referenced in 'right' column does not exist: {right}")
    except Exception as e:
        errors.append(f"Error reading {samples_csv}: {str(e)}")
    return errors



samples = {}
with open(config["sample_file"], "r") as f:
    next(f)
    for line in f:
        sample_id, left, right = line.strip().split(",")
        samples[sample_id] = (left, right)

# Get data from input sample sheet for the rules:
sample_file = config["sample_file"]
#left_files, right_files = extract_filenames_from_csv(sample_file)
fastq_dataframe = create_sample_dataframe(sample_file)
fastq_directory = get_unique_paths_without_extension(fastq_dataframe)
print("fastq dataframe: ", fastq_dataframe)


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

print("sample_keys: ", list(samples.keys()))



def generate_files(input_csv, output_dir):
    # Create the output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Initialize list for meta.txt data
    meta_data = []
    config_data = [
        "out_dir = 'RNAseqCNV_output'",
        "count_dir = 'data/single_counts'",
        "snv_dir = 'data/vcf_files'"
    ]

    # Read the sample.csv file
    with open(input_csv, 'r') as csvfile:
        csvreader = csv.DictReader(csvfile)
        for row in csvreader:
            sample_id = row['sample_id']
            # Add data to meta.txt
            meta_data.append(f"{sample_id}\t{sample_id}.txt\t{sample_id}.tsv")

    # Write data to meta.txt
    meta_file = os.path.join(output_dir, 'meta.txt')
    with open(meta_file, 'w') as meta:
        meta.write("\n".join(meta_data))

    # Write data to config.txt
    config_file = os.path.join(output_dir, 'config.txt')
    with open(config_file, 'w') as config:
        config.write("\n".join(config_data))

    print(f"Meta.txt and config.txt files generated in {output_dir}")

output_directory = 'data/'
generate_files(sample_file, output_directory)




rule all:
    input:
        "check_samples.txt",
        expand("STAR_output/{sample_id}/Aligned.sortedByCoord.out.bam",sample_id=list(samples.keys())),
        #expand("data/single_counts/{sample_id}.txt", sample_id=list(samples.keys())),
        #expand("fusions/{sample_id}.tsv",sample_id=list(samples.keys())),
        #expand("fastqc/{sample}", sample=fastq_dataframe['sample_id']),
        #expand("multiqc/{sample}/multiqc_report.html",sample=fastq_dataframe['sample_id']),
        expand("multiqc/{sample}/multiqc_data/multiqc_fastqc.txt", sample=fastq_dataframe['sample_id']),
        #expand("data/counts/{sample_id}.tsv",sample_id= samples.keys()),
        #expand("allcatch_output/{sample_id}/predictions.tsv",sample_id= samples.keys()),
        expand("ctat_output_directory/{sample_id}/{sample_id}.cancer.vcf",sample_id=samples_test.keys()), # Funktioniert
        #expand("fusions/{sample_id}.tsv",sample_id=samples_test.keys()),
        #expand("data/single_counts/{sample_id}.txt", sample_id=samples.keys()),
        #expand("data/vcf_files/{sample_id}.tsv", sample_id=list(samples.keys())),
        expand("RNAseqCNV_output/{sample_id}",sample_id=samples.keys()),
        #"data/single_counts",
        #"fusioncatcher_output/",
        #expand("fusioncatcher_output/{sample_id}",sample_id= samples.keys()),
        expand("data/total_mapped_reads/{sample_id}.txt", sample_id= samples.keys()),
        expand("data/tpm/{sample_id}.tsv", sample_id=list(samples.keys())),
        #expand("STAR_output/{sample}/Aligned.sortedByCoord.out.bam.bai",sample=list(samples.keys())),
        #expand("pysamstats_output_dir/{sample_id}.coverage.txt",sample_id=list(samples.keys())),
        expand("comparison/{sample_id}.csv", sample_id= samples.keys()),
        expand("aggregated_output/{sample}.csv", sample=list(samples.keys()))



rule check_samples:
    input:
        samples_csv= "samples.csv"
    output:
        "check_samples.txt"
    params:
        left_col=1,
        right_col=2
    run:
        errors = validate_input(input.samples_csv, params.left_col, params.right_col)
        with open(output[0], "w") as output_file:
            if errors:
                for error in errors:
                    output_file.write(f"Error: {error}\n")
            else:
                output_file.write("Sample format and file existence checks passed")




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
        multiqc_fqc= "multiqc/{sample}/multiqc_data/multiqc_fastqc.txt"

    params:
        output_dir= "multiqc/{sample}/"

    shell:
        "multiqc {input} -o {params.output_dir}"



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
        data_directory = config["rna_fusion_data_directory"]

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


#Create merged reads for running final ALLCatchR on all counts
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

#cp -L files tmp/ && rm files && cp tmp/files .
# Rule to run ctat-mutations
rule run_ctat_mutations:
    input:
        input_directory= config["ctat_input_directory"]


    params:
        sample_id= lambda wildcards: wildcards.sample_id,
        genome_lib=config["genome_lib"],
        left= lambda wildcards: samples_test[wildcards.sample_id]['left'],
        right= lambda wildcards: samples_test[wildcards.sample_id]['right'],
        input_directory= config["ctat_input_directory"]


    wildcard_constraints:
        sample_id="|".join(samples_test.keys())

    output:
        vcf ="ctat_output_directory/{sample_id}/{sample_id}.cancer.vcf",
        directory= directory("ctat_output_directory/{sample_id}/")


    shell:
        '''
        singularity exec -e -B {input.input_directory}:/ctat_input \
        -B {params.genome_lib}:/ctat_genome_lib_dir:ro \
        ctat_mutations.v3.2.0.simg \
        /usr/local/src/ctat-mutations/ctat_mutations \
        --left /ctat_input/{params.left} \
        --right /ctat_input/{params.right} \
        --sample_id {params.sample_id} \
        --output {output.directory} \
        --cpu 10 \
        --genome_lib_dir /ctat_genome_lib_dir \
        --boosting_method none \
        --no_cravat &&
        sleep .10 &&
        cp -rL /home/nadine/ctat_output_directory/{params.sample_id}/{params.sample_id}.cancer.vcf /media/nadine/Alina/Blast-o-Matic-Fusioninator/ctat_output_directory/{params.sample_id} &&
        cp -rL /home/nadine/ctat_output_directory/{params.sample_id}/{params.sample_id}.cancer.tsv /media/nadine/Alina/Blast-o-Matic-Fusioninator/ctat_output_directory/{params.sample_id}
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



rule calculate_total_mapped_reads:
    input:
        bam="STAR_output/{sample_id}/Aligned.sortedByCoord.out.bam"
    output:
        total_mapped_reads="data/total_mapped_reads/{sample_id}.txt"
    shell:
        """
        sample_id="{wildcards.sample_id}"
        total_mapped_reads=$(samtools view -F 4 -c "STAR_output/${{sample_id}}/Aligned.sortedByCoord.out.bam")
        echo "$total_mapped_reads" > "data/total_mapped_reads/${{sample_id}}.txt"
        """



rule calculate_tpm:
    input:
        reads_per_gene="STAR_output/{sample_id}/ReadsPerGene.out.tab",
        total_mapped_reads="data/total_mapped_reads/{sample_id}.txt",
        py_script="scripts/calculate_tpm.py"
        
    output:
        tpm="data/tpm/{sample_id}.tsv"

    shell:
        "python {input.py_script} {input.reads_per_gene} {input.total_mapped_reads} {output.tpm}"
   
        


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
        config_file = "data/config.txt",
        metadata_file = "data/meta.txt"

    output:
        directory="RNAseqCNV_output/{sample_id}/",
        rna_seq_cnv_alteration_file="RNAseqCNV_output/{sample_id}/estimation_table.tsv",
        rna_seq_cnv_log2foldchange_file="RNAseqCNV_output/{sample_id}/log2_fold_change_per_arm.tsv",
        rna_seq_cnv_manual_an_table_file= "RNAseqCNV_output/{sample_id}/manual_an_table.tsv"

    shell:
        "Rscript {input.r_script} {input.config_file} {input.metadata_file} {output.directory};"
        "sleep 0.10;"
        "mkdir {output.directory};"
        "mv estimation_table.tsv {output.directory};"
        "mv manual_an_table.tsv {output.directory};"
        "mv log2_fold_change_per_arm.tsv {output.directory};"
        "mv alteration_matrix.tsv {output.directory};"
        "mv {sample_id}/{sample_id}_CNV_main_fig.png  {output.directory}"



rule check_subtype_and_karyotype:
    input:
        prediction_file = "allcatch_output/{sample}/predictions.tsv",
        rna_seq_cnv_estimation_file="RNAseqCNV_output/{sample}/estimation_table.tsv",
        fusioncatcher_file= "fusioncatcher_output/{sample}/final-list_candidate-fusion-genes.hg19.txt",
        arriba_file= "fusions/{sample}.tsv",
        chromosome_counts_karyotype_file= "annotation/chromosome_counts_vs_subtype.txt",
        anno_gene_fusions_file= "annotation/anno_Gene_fusions_vs_Subtypes.txt",
        r_script= "scripts/define_subtype_and_caryotype.R"

    output:
        csv= "comparison/{sample}.csv"

    shell:
        """
            Rscript {input.r_script} {input.prediction_file} {input.rna_seq_cnv_estimation_file} {input.fusioncatcher_file} {input.arriba_file} {input.chromosome_counts_karyotype_file} {input.anno_gene_fusions_file} {output.csv}
        """


rule aggregate_output:
    input:
        prediction_file = "allcatch_output/{sample}/predictions.tsv",
        fusioncatcher_file = "fusioncatcher_output/{sample}/final-list_candidate-fusion-genes.txt",
        arriba_file = "fusions/{sample}.tsv",
        rna_seq_cnv_log2foldchange_file= "RNAseqCNV_output/{sample}/log2_fold_change_per_arm.tsv",
        rna_seq_cnv_manual_an_table_file= "RNAseqCNV_output/{sample}/manual_an_table.tsv",
        star_log_final_out_file= "STAR_output/{sample}/Log.final.out",
        multiqc_fqc_right= "multiqc/{sample}_right/multiqc_data/multiqc_fastqc.txt",
        multiqc_fqc_left= "multiqc/{sample}_left/multiqc_data/multiqc_fastqc.txt",
        comparison_file= "comparison/{sample}.csv"


    output:
        csv= "aggregated_output/{sample}.csv"


    shell:
        """
            cat {input.comparison_file} >> {output.csv}
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
            echo -e "5’ gene name\t5’ chr.position\t3’ gene name\t3’chr. position\tdiscordant_mates"  >> {output.csv}
            awk -F'\t' '{{if (NR == 1) next; printf "%s\\t%s\\t%s\\t%s\\t%s\\n", $1, $5, $2, $6, $12}}' {input.arriba_file} >> {output.csv}
            echo -e "RNASeqCNV identified the following karyotype:" >> {output.csv}
            echo -e "gender	chrom_n	alterations	status	comments"  >> {output.csv}
            awk -F'\t' '{{if (NR == 1) next; printf "%s\\t%s\\t%s\\t%s\\t%s\\n", $2, $3, $4, $5, $6}}' {input.rna_seq_cnv_manual_an_table_file} >> {output.csv}
            echo -e "Chromosome arm calls" >> {output.csv}
            echo -e "chr\tarm\tlog2FC per arm"  >> {output.csv}
            awk -F'\t' '{{if (NR == 1) next; printf "%s\\t%s\\t%s\\n", $1, $2, $3}}' {input.rna_seq_cnv_log2foldchange_file} >> {output.csv}
        """




'''
Inhaltlich offene Punkte:
Integration von karytype und gene expressionsprofil -> machine learning classifier
Definition einer Mutations-Positiv-Liste: CITAT vs. Pysamstat
'''