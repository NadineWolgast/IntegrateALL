configfile: "config.yaml"

# Import custom Python scripts
from scripts.create_sample_dataframe import create_sample_dataframe
from scripts.merge_reads_per_gene_files import merge_reads_per_gene_files
from scripts.generate_files import generate_files
from scripts.validate_input import validate_input
from scripts.get_ctat_input_files import get_ctat_input_files


def creatingfolders(specificfolder: str) -> str:
    """
    As the name suggest it will create a folder if the folder do not exist. it will also check if the end is '/' and add
     it if not there. it will return the folder it created

    :param specificfolder: The folder needs to be created
    :return: will not return anything. Either it will create if the folder do not exist or not return anything
    """
    import os
    if specificfolder != '':
        if specificfolder[-1] != '/':
            specificfolder = specificfolder + '/'

        specificfolder = os.path.expanduser(specificfolder)
        if not os.path.exists(specificfolder):
            os.makedirs(specificfolder)
    return specificfolder

creatingfolders('STAR_output')
creatingfolders('data/combined_counts')

# Get data from input sample sheet for the rules:
sample_file = config["sample_file"]
samples = {}
with open(sample_file, "r") as f:
    next(f)
    for line in f:
        sample_id, left, right = line.strip().split(",")
        samples[sample_id] = (left, right)


# Get FASTQs for QC
fastq_dataframe = create_sample_dataframe(sample_file)

#Get FASTQ's without path for CTAT
samples_test = get_ctat_input_files(sample_file)

#Create config.txt and meta.txt for RNASeqCNV
generate_files(sample_file, 'data/')

rule all:
    input:
        "check_samples.txt",
        expand("STAR_output/{sample_id}/Aligned.sortedByCoord.out.bam",sample_id=list(samples.keys())),
        expand("STAR_output/{sample_id}/Aligned.sortedByCoord.out.bam.bai", sample_id=list(samples.keys())),
        expand("fusions/{sample_id}.pdf",sample_id=samples.keys()),
        #expand("multiqc/{sample}/multiqc_data/multiqc_fastqc.txt", sample=fastq_dataframe['sample_id']),
        #expand("fusioncatcher_output/{sample_id}/final-list_candidate-fusion-genes.txt",sample_id=list(samples.keys())),
        #expand("ctat_output_directory/{sample_id}/{sample_id}.filtered.vcf.gz",sample_id=samples_test.keys()),
        #expand("RNAseqCNV_output/{sample_id}",sample_id=samples.keys()),
        expand("data/tpm/{sample_id}.tsv", sample_id=list(samples.keys())),
        expand("data/cpm/{sample_id}.tsv", sample_id=list(samples.keys())),
        #expand("pysamstats_output_dir/{sample_id}/", sample_id=list(samples.keys())),
        #expand("comparison/{sample_id}.csv", sample_id= samples.keys()),
        #expand("data/single_counts/{sample_id}.txt", sample_id=samples.keys()),
        #expand("data/vcf_files/{sample_id}.tsv",sample_id=samples.keys()),
        #expand("allcatch_output/{sample_id}/predictions.tsv", sample_id= samples.keys()),
        expand("aggregated_output/{sample}.csv", sample=list(samples.keys()))
        #expand("interactive_output/output_report_{sample}.html",  sample=list(samples.keys()))



rule check_samples:
    input:
        samples_csv= "samples.csv"
    output:
        "check_samples.txt"

    run:
        errors = validate_input(input.samples_csv)
        with open(output[0], "w") as output_file:
            if errors:
                for error in errors:
                    output_file.write(f"Error: {error}\n")
            else:
                output_file.write("Sample format and file existence checks passed")



# Function to get input FASTQ files based on sample_id
def get_input_fastqs(wildcards):
    sample_name = wildcards.sample
    fastq_file = fastq_dataframe[fastq_dataframe['sample_id'] == sample_name]['FASTQ'].values[0]
    return fastq_file


rule fastqc:
    input:
        get_input_fastqs
    output:
        html="fastqc/{sample}.html",
        zip="fastqc/{sample}.zip"
    log:
        "logs/fastqc/{sample}/fastqc.log"

    benchmark:
        "benchmarks/{sample}.fastqc.benchmark.txt"

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
        "multiqc {input} -o {params.output_dir} --force"


rule multiqc_file:
    input:
        "fastqc/{sample}/"
    output:
        multiqc_report="multiqc/{sample}/multiqc_report.html",
        multiqc_fqc= "multiqc/{sample}/multiqc_data/multiqc_fastqc.txt"

    params:
        output_dir= "multiqc/{sample}/"

    benchmark:
        "benchmarks/{sample}.multiqc.benchmark.txt"

    shell:
        "multiqc {input} -o {params.output_dir} --force"



rule install_all:
    input: []
    shell:
        """
            snakemake --cores 2 download_star_ref &&
            snakemake --cores all index &&
            snakemake --cores 1 install_arriba_draw_fusions &&
            snakemake --cores 2 install_allcatchr &&
            snakemake --cores 1 pull_ctat_mutations_singularity_image &&
            snakemake --cores 2 install_ctat_mutations &&
            snakemake --use-singularity --cores 2 run_ctat_genome_lib_builder &&
            snakemake --use-conda --cores 2 install_rnaseq_cnv &&
            snakemake --cores 2 install_fusioncatcher
        """


rule install_fusioncatcher:
    params:
        data_directory=config["rna_fusion_data_directory"]
    shell:
        """
            conda config --add channels defaults &&
            conda config --add channels bioconda &&
            conda config --add channels conda-forge &&
            conda create -n fusioncatcher fusioncatcher &&
            source activate fusioncatcher &&
            download-human-db.sh &&
            mv human_v102 {params.data_directory}
        """



rule download_star_ref:
    input:
        star_directory= config["star_files"]

    shell:
        "cd {input.star_directory} && "
        "wget 'https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz' && "
        "gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz && "
        "wget 'https://ftp.ensembl.org/pub/release-94/gtf/homo_sapiens/Homo_sapiens.GRCh38.94.gtf.gz' && "
        "gunzip Homo_sapiens.GRCh38.94.gtf.gz "

#conda install star=2.7.1a

rule index:
        input:
            fa = config['star_ref'], # provide your reference FASTA file
            gtf = config['star_gtf'] # provide your GTF file
        output:
            directory(config["genome_index"])

        threads: config['threads']

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

    benchmark:
        "benchmarks/{sample_id}.star_aligner.benchmark.txt"

    conda:
        "envs/star.yaml"

    threads: config['threads']

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
        '--readFilesCommand zcat && '
        'rm -r {output.directory}/_STARgenome {output.directory}/_STARpass1 '


rule samtools_index:
    input:
        "STAR_output/{sample}/Aligned.sortedByCoord.out.bam"
    output:
        "STAR_output/{sample}/Aligned.sortedByCoord.out.bam.bai"
    log:
        "logs/samtools_index/{sample}.log"
    params:
        extra=""  # optional params string
    benchmark:
        "benchmarks/{sample}.samtools_index.benchmark.txt"
    threads: config['threads']
    wrapper:
        "v2.6.0/bio/samtools/index"


rule pysamstat:
    input:
        bam="STAR_output/{sample_id}/Aligned.sortedByCoord.out.bam",
        fa= config['star_ref']

    benchmark:
        "benchmarks/{sample_id}.pysamstat.benchmark.txt"

    output:
        directory("pysamstats_output_dir/{sample_id}/"),
        ikzf1="pysamstats_output_dir/{sample_id}/{sample_id}_IKZF1.tsv",
        pax5="pysamstats_output_dir/{sample_id}/{sample_id}_PAX5.tsv",
        coverage="pysamstats_output_dir/{sample_id}/example.coverage.txt"


    shell:
        """
        pysamstats --type variation --chromosome 7 -u --start 50382594 --end 50382595 -f {input.fa} {input.bam} > {output.ikzf1} &&
        pysamstats --type variation --chromosome 9 -u --start 37015167 --end 37015168 -f {input.fa} {input.bam} > {output.pax5} &&
        pysamstats --type coverage --chromosome 7 -u --start 50382594 --end 50382595 {input.bam} > {output.coverage}
        """


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
    benchmark:
        "benchmarks/{sample_id}.arriba.benchmark.txt"
    params:
        genome_build="GRCh38",           # Required when blacklist or known_fusions is set
        default_blacklist=False,         # Optional
        default_known_fusions=True,      # Optional
        sv_file="",                      # File containing information from structural variant analysis
        extra="-i 1,2"                  # Optional parameters
    threads: config['threads']
    wrapper:
        "v2.6.0/bio/arriba"


rule install_arriba_draw_fusions:
    shell:
        """
        wget 'https://github.com/suhrig/arriba/releases/download/v2.4.0/arriba_v2.4.0.tar.gz' &&
        tar -xzf arriba_v2.4.0.tar.gz &&
        cd arriba_v2.4.0 &&
        ./download_references.sh hs37d5viral+GENCODE19
        """


rule run_draw_arriba_fusion:
    input:
        fusions = "fusions/{sample_id}.tsv",
        bam="STAR_output/{sample_id}/Aligned.sortedByCoord.out.bam",
        bai="STAR_output/{sample_id}/Aligned.sortedByCoord.out.bam.bai",
        annotation=config["star_gtf"] ,  # Path to annotation GTF
        r_script="scripts/draw_fusions.R"

    output:
        pdf="fusions/{sample_id}.pdf"

    conda:
        "envs/arriba_draw_fusions.yaml"

    shell:
        '''
        Rscript {input.r_script} \
        --fusions={input.fusions} \
        --alignments={input.bam} \
        --output={output.pdf} \
        --annotation={input.annotation}\
        --cytobands=arriba_v2.4.0/database/cytobands_hg19_hs37d5_GRCh37_v2.4.0.tsv \
        --proteinDomains=arriba_v2.4.0/database/protein_domains_hg19_hs37d5_GRCh37_v2.4.0.gff3
        '''



rule run_fusioncatcher:
    input:
        fastq_directory = config["ctat_input_directory"],
        data_directory = config["rna_fusion_data_directory"],
        left= lambda wildcards: samples[wildcards.sample_id][0],
        right= lambda wildcards: samples[wildcards.sample_id][1]

    output:
        dir=directory("fusioncatcher_output/{sample_id}"),
        file="fusioncatcher_output/{sample_id}/final-list_candidate-fusion-genes.txt"

    conda:
        "envs/fusioncatcher.yaml"

    params:
        sample_id=lambda wildcards: wildcards.sample_id


    benchmark:
        "benchmarks/{sample_id}.fusioncatcher.benchmark.txt"

    shell:
        '''
        fusioncatcher \
        -d {input.data_directory} \
        -i {input.left},{input.right} \
        -o {output.dir}'''


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
    benchmark:
        "benchmarks/{sample_id}.reads_per_gene_to_counts.benchmark.txt"
    shell:
        """
        echo -e "Gene\t{wildcards.sample_id}" > {output.counts}
        awk 'NR > {params.skip_rows} {{print $1 "\t" $2}}' {input.reads_per_gene} >> {output.counts}
        """



rule run_allcatchr_on_single_count_files:
    input:
        r_script = "scripts/run_ALLCatchR.R",
        input_file = 'data/counts/{sample}.tsv'

    benchmark:
        "benchmarks/{sample}.allcatchr.benchmark.txt"

    output:
        "allcatch_output/{sample}/predictions.tsv"

    conda:
        "envs/catchall.yaml"

    shell:
        "Rscript {input.r_script} {input.input_file} {output};"
        "mv predictions.tsv {output}"



rule pull_ctat_mutations_singularity_image:
    shell:
        '''wget "https://data.broadinstitute.org/Trinity/CTAT_SINGULARITY/CTAT_MUTATIONS/__archived_releases/ctat_mutations.v3.2.0.simg"'''


rule install_ctat_mutations:
    params:
        software_url="https://github.com/NCIP/ctat-mutations/archive/refs/tags/CTAT-Mutations-v3.2.0.tar.gz",
        genome_lib_url="https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/__genome_libs_StarFv1.10/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play.tar.gz",
        genome_lib_build_dir= config["ctat_genome_lib_build_dir"],
        mutation_lib_supplement_url="https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/MUTATION_LIB_SUPPLEMENT/GRCh38.mutation_lib_supplement.Jul272020.tar.gz",
        rna_editing_url="https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/MUTATION_LIB_SUPPLEMENT/rna_editing/GRCh38.RNAediting.vcf.gz",
        cravat_url = "https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/MUTATION_LIB_SUPPLEMENT/cravat_lib/cravat.GRCh38.tar.bz2"

    shell:
        "mkdir -p {params.genome_lib_build_dir} && cd {params.genome_lib_build_dir} && "

        "wget {params.genome_lib_url} && "
        "tar xvf GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play.tar.gz && "

        "wget {params.software_url} && "
        "tar xvf CTAT-Mutations-v3.2.0.tar.gz && "
        
        "cd {params.genome_lib_build_dir}/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir && "
        "wget {params.mutation_lib_supplement_url} && "
        "tar xvf GRCh38.mutation_lib_supplement.Jul272020.tar.gz && "
        
        "wget {params.rna_editing_url} && "
        "gunzip GRCh38.RNAediting.vcf.gz && "
        
        "wget {params.cravat_url} && "
        "tar xvf cravat.GRCh38.tar.bz2"



rule run_ctat_genome_lib_builder:
    input:
        genome_lib = config["ctat_genome_lib_build_dir"]
    params:
        genome_lib_build_dir= config["ctat_genome_lib_build_dir"]

    shell:
        '''
        singularity exec -e -B {params.genome_lib_build_dir}/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir \
        -B {params.genome_lib_build_dir}/ctat-mutations-CTAT-Mutations-v3.2.0/mutation_lib_prep \
        ctat_mutations.v3.2.0.simg \
        {params.genome_lib_build_dir}/ctat-mutations-CTAT-Mutations-v3.2.0/mutation_lib_prep/ctat-mutation-lib-integration.py \
        --genome_lib_dir {params.genome_lib_build_dir}/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir \
        '''


# Rule to run ctat-mutations
rule run_ctat_mutations:
    input:
        input_directory= config["ctat_input_directory"]

    params:
        sample_id= lambda wildcards: wildcards.sample_id,
        genome_lib_build_dir= config["ctat_genome_lib_build_dir"],
        left= lambda wildcards: samples_test[sample_id]['left'],
        right= lambda wildcards: samples_test[sample_id]['right'],
        input_directory= config["ctat_input_directory"],
        threads= config['threads']

    benchmark:
        "benchmarks/{sample_id}.ctat_mutations.benchmark.txt"

    output:
        vcf="ctat_output_directory/{sample_id}/{sample_id}.filtered.vcf.gz",
        directory= directory("ctat_output_directory/{sample_id}/")

    threads: config['threads']

    shell:
        '''
        singularity exec -e -B {input.input_directory}:/ctat_input \
        -B {params.genome_lib_build_dir}/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir:/ctat_genome_lib_dir:ro \
        -B {params.genome_lib_build_dir}/ctat-mutations-CTAT-Mutations-v3.2.0\
        ctat_mutations.v3.2.0.simg \
        {params.genome_lib_build_dir}/ctat-mutations-CTAT-Mutations-v3.2.0/ctat_mutations \
        --left /ctat_input/{params.left} \
        --right /ctat_input/{params.right} \
        --sample_id {params.sample_id} \
        --output {output.directory} \
        --cpu {config[threads]} \
        --genome_lib_dir /ctat_genome_lib_dir \
        --boosting_method none \
        --no_cravat &&
        sed -i ';' {output.vcf} 
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



rule calculate_tpm_and_cpm:
    input:
        reads_per_gene="STAR_output/{sample_id}/ReadsPerGene.out.tab",
        r_script="scripts/calculate_tpm_and_cpm.R"
        
    output:
        tpm="data/tpm/{sample_id}.tsv",
        cpm="data/cpm/{sample_id}.tsv"

    shell:
        "Rscript {input.r_script} {input.reads_per_gene} {output.tpm} {output.cpm}"
   
        


rule install_rnaseq_cnv:
    conda:
        "envs/rnaseqcnv.yaml"

    shell:
        "Rscript -e 'devtools::install_github(\"honzee/RNAseqCNV\")'"


# TODO: Check if all the cancer.vcf files need to skip 263 lines. Needs maybe a modified skript
rule prepare_vcf_files:
    input:
        vcf_dir="ctat_output_directory/{sample_id}/{sample_id}.filtered.vcf.gz",
        r_script= "scripts/prepare_vcf_files.R"
    output:
        counts="data/vcf_files/{sample_id}.tsv"

    shell:
        # Execute the R script using Rscript command
        "Rscript {input.r_script} {input.vcf_dir} {output.counts};"

#TODO: Need to rewrite meta.txt for each sample! Execute in R script : run_rnaseq_cnv.R",
rule run_rnaseq_cnv:
    input:
        r_script = "scripts/run_rnaseq_cnv.R",
        config_file = "data/config.txt",
        metadata_file = "data/meta.txt",
        input_counts= "data/single_counts/{sample_id}.txt",
        input_vcf= "data/vcf_files/{sample_id}.tsv"

    output:
        directory=directory("RNAseqCNV_output/{sample_id}/"),
        rna_seq_cnv_alteration_file="RNAseqCNV_output/{sample_id}/estimation_table.tsv",
        rna_seq_cnv_log2foldchange_file="RNAseqCNV_output/{sample_id}/log2_fold_change_per_arm.tsv",
        rna_seq_cnv_manual_an_table_file= "RNAseqCNV_output/{sample_id}/manual_an_table.tsv"

    conda:
        "envs/rnaseqcnv.yaml"

    shell:
        "Rscript {input.r_script} {input.config_file} {input.metadata_file} {wildcards.sample_id};"
        "sleep 0.20;"
        "mv RNAseqCNV_output/estimation_table.tsv {output.directory};"
        "mv RNAseqCNV_output/manual_an_table.tsv {output.directory};"
        "mv RNAseqCNV_output/log2_fold_change_per_arm.tsv {output.directory};"
        "mv RNAseqCNV_output/alteration_matrix.tsv {output.directory};"




rule check_subtype_and_karyotype:
    input:
        prediction_file = "allcatch_output/{sample}/predictions.tsv",
        rna_seq_cnv_estimation_file="RNAseqCNV_output/{sample}/estimation_table.tsv",
        fusioncatcher_file = "fusioncatcher_output/{sample}/final-list_candidate-fusion-genes.txt",
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
            rna_seq_cnv_log2foldchange_file = "RNAseqCNV_output/{sample}/log2_fold_change_per_arm.tsv",
            rna_seq_cnv_manual_an_table_file = "RNAseqCNV_output/{sample}/manual_an_table.tsv",
            star_log_final_out_file = "STAR_output/{sample}/Log.final.out",
            multiqc_fqc_right = "multiqc/{sample}_right/multiqc_data/multiqc_fastqc.txt",
            multiqc_fqc_left = "multiqc/{sample}_left/multiqc_data/multiqc_fastqc.txt",
            comparison_file = "comparison/{sample}.csv"

    output:
        csv="aggregated_output/{sample}.csv"

    shell:
        """
        bash scripts/process_data.sh {output.csv} {input.star_log_final_out_file} {wildcards.sample} {input.multiqc_fqc_left} {input.multiqc_fqc_right} {input.prediction_file} {input.fusioncatcher_file} {input.arriba_file} {input.rna_seq_cnv_manual_an_table_file} {input.rna_seq_cnv_log2foldchange_file} {input.comparison_file}
        """




rule interactive_report:
    input:
        prediction_file="allcatch_output/{sample}/predictions.tsv",
        fusioncatcher_file="fusioncatcher_output/{sample}/final-list_candidate-fusion-genes.txt",
        arriba_file="fusions/{sample}.tsv",
        arriba_file_fusion="fusions/{sample}.pdf",
        rna_seq_cnv_log2foldchange_file="RNAseqCNV_output/{sample}/log2_fold_change_per_arm.tsv",
        rna_seq_cnv_plot="RNAseqCNV_output/{sample}/{sample}_CNV_main_fig.png",
        rna_seq_cnv_manual_an_table_file="RNAseqCNV_output/{sample}/manual_an_table.tsv",
        star_log_final_out_file="STAR_output/{sample}/Log.final.out",
        multiqc_fqc_right="multiqc/{sample}_right/multiqc_report.html",
        multiqc_fqc_left="multiqc/{sample}_left/multiqc_report.html",
        comparison_file="comparison/{sample}.csv",
        pysamstats_files_IKZF1="pysamstats_output_dir/{sample}/{sample}_IKZF1.csv",
        pysamstats_files_PAX5="pysamstats_output_dir/{sample}/{sample}_PAX5_P80R.tsv",
        pysamstats_files_coverage="pysamstats_output_dir/{sample}/example.coverage.txt",
        hotspots="pysamstats_output_dir/{sample}/Hotspots/"

    output:
        html="interactive_output/{sample}/output_report_{sample}.html"

    shell:
        """
        mkdir interactive_output/{wildcards.sample}/fusions &&
        cp {input.arriba_file_fusion} interactive_output/{wildcards.sample}/fusions &&

        mkdir interactive_output/{wildcards.sample}/multiqc_right &&
        cp {input.multiqc_fqc_right} interactive_output/{wildcards.sample}/multiqc_right &&

        mkdir interactive_output/{wildcards.sample}/multiqc_left &&
        cp {input.multiqc_fqc_left} interactive_output/{wildcards.sample}/multiqc_left &&

        mkdir interactive_output/{wildcards.sample}/RNAseqCNV &&
        cp {input.rna_seq_cnv_plot} interactive_output/{wildcards.sample}/RNAseqCNV &&

        python scripts/generate_report.py  {input.prediction_file} {input.fusioncatcher_file} {input.arriba_file} {input.rna_seq_cnv_log2foldchange_file} {input.rna_seq_cnv_manual_an_table_file} {input.star_log_final_out_file}  {input.comparison_file} {input.pysamstats_files_IKZF1} {input.pysamstats_files_PAX5} {input.pysamstats_files_coverage} {input.hotspots} {wildcards.sample} {output.html}
        """

