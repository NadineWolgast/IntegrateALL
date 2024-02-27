configfile: "config.yaml"

# Import custom Python scripts
from scripts.create_sample_dataframe import create_sample_dataframe
from scripts.generate_files import generate_files
from scripts.validate_input import validate_input
from scripts.get_ctat_input_files import get_ctat_input_files


def creatingfolders(specificfolder: str) -> str:
    """
    As the name suggest it will create a folder if the folder does not exist. it will also check if the end is '/' and add
     it if not there. it will return the folder it created

    :param specificfolder: The folder needs to be created
    :return: will not return anything. Either it will create if the folder does not exist or not return anything
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
creatingfolders('Variants_RNA_Seq_Reads')
creatingfolders('data/combined_counts')
creatingfolders('refs/fusioncatcher')
#creatingfolders('refs/STAR')
creatingfolders('refs/GATK')

sample_file = config["sample_file"]
samples = {}
with open(sample_file,"r") as f:
    next(f)
    for line in f:
        sample_id, left, right = line.strip().split(",")
        samples[sample_id] = (left, right)

# Get FASTQs for QC
fastq_dataframe = create_sample_dataframe(sample_file)

#Create config.txt and meta.txt for RNASeqCNV
generate_files(sample_file, 'data/')

absolute_path = config["absolute_path"]

rule all:
    input:
        "check_samples.txt",
        expand("STAR_output/{sample_id}/Aligned.sortedByCoord.out.bam",sample_id=list(samples.keys())),
        expand("STAR_output/{sample_id}/Aligned.sortedByCoord.out.bam.bai", sample_id=list(samples.keys())),
        #expand("Variants_RNA_Seq_Reads/{sample_id}/fixed-rg/{sample_id}.bam", sample_id=list(samples.keys())),
        #expand("Variants_RNA_Seq_Reads/{sample_id}/deduped_bam/{sample_id}.bam.bai", sample_id=list(samples.keys())),
        #expand("Variants_RNA_Seq_Reads/{sample_id}/split/{sample_id}.bam", sample_id=list(samples.keys())),
        #expand("Variants_RNA_Seq_Reads/{sample_id}/recal/{sample_id}_recal.table", sample_id=list(samples.keys())),
        #expand("Variants_RNA_Seq_Reads/{sample_id}/recal/{sample_id}.bam", sample_id=list(samples.keys())),
        expand("Variants_RNA_Seq_Reads/{sample_id}/filter/{sample_id}.snvs.filtered.vcf", sample_id=list(samples.keys())),
        expand("fusions/{sample_id}.pdf",sample_id=samples.keys()),
        expand("fusions/{sample_id}.tsv",sample_id=samples.keys()),

        expand("multiqc/{sample}/multiqc_data/multiqc_fastqc.txt", sample=fastq_dataframe['sample_id']),
        expand("fusioncatcher_output/{sample_id}/final-list_candidate-fusion-genes.txt",sample_id=list(samples.keys())),
        expand("data/vcf_files/GATK/{sample_id}_Gatk.tsv", sample_id=samples.keys()),     
        expand("data/tpm/{sample_id}.tsv", sample_id=list(samples.keys())),
        #expand("data/cpm/{sample_id}.tsv", sample_id=list(samples.keys())),        
        expand("pysamstats_output_dir/{sample_id}/", sample_id=list(samples.keys())),
        expand("Hotspots/{sample_id}",sample_id=list(samples.keys())),
        expand("comparison/{sample_id}.csv", sample_id= samples.keys()),
        expand("data/single_counts/{sample_id}.txt", sample_id=samples.keys()),
        expand("allcatch_output/{sample_id}/predictions.tsv", sample_id= samples.keys()),        
        expand("aggregated_output/{sample}.csv", sample=list(samples.keys())),        
        expand("RNAseqCNV_output/gatk/{sample_id}_gatk", sample_id=samples.keys()),        
        expand("interactive_output/{sample}/output_report_{sample}.html",  sample=list(samples.keys()))




rule check_samples:
    input:
        samples_csv= config["sample_file"]
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
            snakemake --cores 10 download_ref &&
            snakemake --cores 10 index_ref &&
            snakemake --cores all index_star &&
            snakemake --cores 1 install_arriba_draw_fusions &&
            snakemake --use-conda --cores 2 install_allcatchr &&
            snakemake --use-conda --cores 2 install_rnaseq_cnv &&
            snakemake --cores 2 install_fusioncatcher
        """


rule download_ref:
    input:
        star_directory = absolute_path + "/refs/GATK"
    output:
        vcf= "refs/GATK/GRCH38/dbSNP.vcf",
        ref= "refs/GATK/GRCH38/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    shell:
        "cd {input.star_directory} && "
        "wget 'http://141.2.194.197/rnaeditor_annotations/GRCH38.tar.gz' && "
        "tar -xzf GRCH38.tar.gz "


rule index_ref:
    input: "refs/STAR/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    output: "refs/STAR/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai"
    shell:
        "samtools faidx {input}"


rule index_star:
        input:
            fa = absolute_path + "/refs/GATK/GRCH38/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
            gtf = absolute_path + "/refs/GATK/GRCH38/Homo_sapiens.GRCh38.83.gtf"
        output:
            directory(absolute_path + "/refs/GATK/STAR/ensembl_94_100")

        threads: config['threads']

        shell:
            'STAR --runThreadN {threads} '
            '--runMode genomeGenerate '
            '--genomeDir {output} '
            '--genomeFastaFiles {input.fa} '
            '--sjdbGTFfile {input.gtf} '
            '--sjdbOverhang 100'



rule install_arriba_draw_fusions:
    shell:
        '''
        wget 'https://github.com/suhrig/arriba/releases/download/v2.4.0/arriba_v2.4.0.tar.gz' &&
        tar -xzf arriba_v2.4.0.tar.gz &&
        cd arriba_v2.4.0 &&
        ./download_references.sh hs37d5viral+GENCODE19
        '''


rule install_allcatchr:
    conda:
    	"envs/install_allcatchr.yaml"

    shell:
        "Rscript -e 'devtools::install_github(\"ThomasBeder/ALLCatchR\")'"


rule install_rnaseq_cnv:
    conda:
        "envs/rnaseqenv.yaml"

    shell:
        "Rscript -e 'devtools::install_github(\"honzee/RNAseqCNV\")'"
        

rule install_fusioncatcher:
    input:
        data_directory=absolute_path + "/refs/fusioncatcher"
    output:
        directory(absolute_path + "/refs/fusioncatcher/fusioncatcher-master/data/human_v102")

    shell:
        """
        cd  {input.data_directory} &&
        wget 'https://github.com/ndaniel/fusioncatcher/archive/refs/heads/master.zip' && 
        unzip master.zip && 
        cd fusioncatcher-master/data && 
        ./download-human-db.sh
        """


rule run_star_aligner:
    input:
        fastq1 = lambda wildcards: samples[wildcards.sample_id][0],  # Path to left FASTQ from sample sheet
        fastq2 = lambda wildcards: samples[wildcards.sample_id][1],  # Path to right FASTQ from sample sheet
        genome_index = absolute_path + "/refs/GATK/STAR/ensembl_94_100"
    params:
    	

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
        '--genomeDir {input.genome_index} '
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
        bai="STAR_output/{sample_id}/Aligned.sortedByCoord.out.bam.bai",
        fa= absolute_path + "/refs/GATK/GRCH38/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    params:
        out_dir="pysamstats_output_dir/{sample_id}/"

    output:
        pysamstats_output_dir = temporary(directory("pysamstats_output_dir/{sample_id}/")),
        ikzf1="pysamstats_output_dir/{sample_id}/{sample_id}_IKZF1.csv",

    shell:
        """
        mkdir -p {output.pysamstats_output_dir} &&
        pysamstats --type variation --chromosome 7 -u --start 50382593 --end 50382596 -f {input.fa} {input.bam} > {output.ikzf1} &&
        python scripts/run_pysamstats.py  {input.bam} {input.fa} {wildcards.sample_id} {params.out_dir} 
        """



rule get_Hotspots:
    input:
        pysamstats_output_dir="pysamstats_output_dir/{sample_id}/",
        r_script="scripts/Get_Amino_for_Hotspot.R",
        gatk_file = "Variants_RNA_Seq_Reads/{sample_id}/filter/{sample_id}.snvs.filtered.vcf"
    output:
        hotspot_output_dir=directory("Hotspots/{sample_id}") 
    shell:
        """
        mkdir -p {output.hotspot_output_dir} &&
        Rscript {input.r_script} {input.pysamstats_output_dir} {input.gatk_file} {output.hotspot_output_dir}
        """




rule run_arriba:
    input:
        bam="STAR_output/{sample_id}/Aligned.sortedByCoord.out.bam",
        genome=absolute_path + "/refs/GATK/GRCH38/Homo_sapiens.GRCh38.dna.primary_assembly.fa",       # Path to reference genome
        annotation=absolute_path + "/refs/GATK/GRCH38/Homo_sapiens.GRCh38.83.gtf",   # Path to annotation GTF
        custom_blacklist=[]
    output:
        # Approved gene fusions
        fusions="fusions/{sample_id}.tsv",
        # Discarded gene fusions (optional)
        discarded=temporary("fusions/{sample_id}.discarded.tsv")
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


rule run_draw_arriba_fusion:
    input:
        fusions = "fusions/{sample_id}.tsv",
        bam="STAR_output/{sample_id}/Aligned.sortedByCoord.out.bam",
        bai="STAR_output/{sample_id}/Aligned.sortedByCoord.out.bam.bai",
        annotation=absolute_path + "/refs/GATK/GRCH38/Homo_sapiens.GRCh38.83.gtf",
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
        data_directory= absolute_path + "/refs/fusioncatcher/fusioncatcher-master/data/human_v102",
        left= lambda wildcards: samples[wildcards.sample_id][0],
        right= lambda wildcards: samples[wildcards.sample_id][1]

    output:
        dir=directory("fusioncatcher_output/{sample_id}"),
        keep_file="fusioncatcher_output/{sample_id}/final-list_candidate-fusion-genes.txt"

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
        -o {output.dir} \
        -p 20'''


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



rule run_allcatchr:
    input:
        r_script = "scripts/run_ALLCatchR.R",
        input_file = 'data/counts/{sample}.tsv'

    benchmark:
        "benchmarks/{sample}.allcatchr.benchmark.txt"

    output:
        "allcatch_output/{sample}/predictions.tsv"

    conda:
        "envs/install_allcatchr.yaml"
        
    benchmark:
        "benchmarks/{sample_id}.allcatchr.benchmark.txt"

    shell:
        "Rscript {input.r_script} {input.input_file} {output};"
        "mv predictions.tsv {output}"





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
   
        



rule replace_rg:
    input:
        "STAR_output/{sample}/Aligned.sortedByCoord.out.bam"

    output:
        temporary("Variants_RNA_Seq_Reads/{sample}/fixed-rg/{sample}.bam")

    benchmark:
        "benchmarks/{sample_id}.replace_rg.benchmark.txt"

    log:
        "logs/picard/replace_rg/{sample}.log"
    params:
        extra="--RGLB lib1 --RGPL illumina --RGPU {sample} --RGSM {sample}",
        java_opts=""

    resources:
        mem_mb=2048,
    wrapper:
        "v3.3.3/bio/picard/addorreplacereadgroups"


rule markduplicates_bam:
    input:
        bams="Variants_RNA_Seq_Reads/{sample}/fixed-rg/{sample}.bam"

    output:
        bam=temporary("Variants_RNA_Seq_Reads/{sample}/deduped_bam/{sample}.bam"),
        metrics=temporary("Variants_RNA_Seq_Reads/{sample}/deduped_bam/{sample}.metrics.txt")

    benchmark:
        "benchmarks/{sample_id}.markduplicates_bam.benchmark.txt"

    log:
        "logs/picard/dedup_bam/{sample}.log",

    resources:
        mem_mb=5000,
    wrapper:
        "v3.3.3/bio/picard/markduplicates"


rule index_bam:
    input:
        "Variants_RNA_Seq_Reads/{sample}/deduped_bam/{sample}.bam"
    output:
        temporary("Variants_RNA_Seq_Reads/{sample}/deduped_bam/{sample}.bam.bai")
    log:
        "logs/samtools_index/{sample}_deduped.log"
    params:
        extra=""  # optional params string
    benchmark:
        "benchmarks/{sample}.samtools_index.benchmark.txt"

    threads: config['threads']
    wrapper:
        "v2.6.0/bio/samtools/index"


rule splitncigarreads:
    input:
        bam="Variants_RNA_Seq_Reads/{sample}/deduped_bam/{sample}.bam",
        bai="Variants_RNA_Seq_Reads/{sample}/deduped_bam/{sample}.bam.bai",
        ref= "refs/GATK/GRCH38/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    output:
        temporary("Variants_RNA_Seq_Reads/{sample}/split/{sample}.bam")

    benchmark:
        "benchmarks/{sample}.splitncigarreads.benchmark.txt"

    log:
        "logs/gatk/splitNCIGARreads/{sample}.log",
    params:
        extra="",  # optional
        java_opts="",  # optional
    resources:
        mem_mb=5000,
    threads: config['threads']
    wrapper:
        "v3.3.3/bio/gatk/splitncigarreads"




rule gatk_baserecalibrator:
    input:
    	bam="Variants_RNA_Seq_Reads/{sample}/split/{sample}.bam",
        ref= "refs/GATK/GRCH38/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
        dict= "refs/GATK/GRCH38/Homo_sapiens.GRCh38.dna.primary_assembly.dict",
        known= "refs/GATK/GRCH38/dbSNP.vcf",
    output:
        recal_table=temporary("Variants_RNA_Seq_Reads/{sample}/recal/{sample}_recal.table")

    benchmark:
        "benchmarks/{sample}.gatk_baserecalibrator.benchmark.txt"

    log:
        "logs/gatk/baserecalibrator/{sample}.log",
    params:
        extra="",  # optional
        java_opts="",  # optional
    resources:
        mem_mb=5000,
    threads: config['threads']
    wrapper:
        "v3.3.3/bio/gatk/baserecalibrator"


rule gatk_applybqsr:
    input:
        bam="Variants_RNA_Seq_Reads/{sample}/split/{sample}.bam",
        ref="refs/GATK/GRCH38/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
        dict="refs/GATK/GRCH38/Homo_sapiens.GRCh38.dna.primary_assembly.dict",
        recal_table="Variants_RNA_Seq_Reads/{sample}/recal/{sample}_recal.table"
    output:
        bam=temporary("Variants_RNA_Seq_Reads/{sample}/recal/{sample}.bam")

    benchmark:
        "benchmarks/{sample}.gatk_applybqsr.benchmark.txt"

    log:
        "logs/gatk/gatk_applybqsr/{sample}.log",
    params:
        extra="",  # optional
        java_opts="",  # optional
        embed_ref=True,  # embed the reference in cram output
    resources:
        mem_mb=5000,
    threads: config['threads']
    wrapper:
        "v3.3.3/bio/gatk/applybqsr"



rule haplotype_caller:
    input:
        bam= "Variants_RNA_Seq_Reads/{sample}/recal/{sample}.bam",
        ref= "refs/GATK/GRCH38/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
        #known= "refs/STAR/dbSNP.vcf" #optional

    output:
        vcf=temporary("Variants_RNA_Seq_Reads/{sample}/calls/{sample}.vcf")

    benchmark:
        "benchmarks/{sample}.haplotype_caller.benchmark.txt"

    log:
        "logs/gatk/haplotypecaller/{sample}.log",
    params:
        extra="-ERC GVCF --output-mode EMIT_ALL_CONFIDENT_SITES --dont-use-soft-clipped-bases -stand-call-conf 20.0",  # optional
        java_opts="",  # optional
    threads: config['threads']
    resources:
        mem_mb=5000,
    wrapper:
        "v3.3.3/bio/gatk/haplotypecaller"


rule gatk_filter:
    input:
        vcf="Variants_RNA_Seq_Reads/{sample}/calls/{sample}.vcf",
        ref= "refs/GATK/GRCH38/Homo_sapiens.GRCh38.dna.primary_assembly.fa",

    output:
        vcf="Variants_RNA_Seq_Reads/{sample}/filter/{sample}.snvs.filtered.vcf"

    benchmark:
        "benchmarks/{sample}.gatk_filter.benchmark.txt"

    log:
        "logs/gatk/filter/{sample}.snvs.log",
    params:
        filters={"myfilter": "AB < 0.2 || MQ0 > 50"},
        extra="",  # optional arguments, see GATK docs
        java_opts="",  # optional
    resources:
        mem_mb=5000,
    threads: config['threads']
    wrapper:
        "v3.3.3/bio/gatk/variantfiltration"


rule prepare_vcf_files_from_GATK:
    input:
        vcf_dir="Variants_RNA_Seq_Reads/{sample_id}/filter/{sample_id}.snvs.filtered.vcf",
        r_script="scripts/prepare_vcf-files_gatk.R"
    output:
        tsv="data/vcf_files/GATK/{sample_id}_Gatk.tsv"
    conda:
        "envs/gatk.yaml"
    shell:
        "grep ':AD:' {input.vcf_dir} > {input.vcf_dir}_sel && "
        "Rscript {input.r_script} {input.vcf_dir} {output.tsv};"
            


rule run_rnaseq_cnv_gatk:
    input:
        r_script="scripts/run_rnaseq_cnv_gatk.R",
        config_file="data/config.txt",
        metadata_file="data/meta.txt",
        input_counts="data/single_counts/{sample_id}.txt",
        input_tsv="data/vcf_files/GATK/{sample_id}_Gatk.tsv"

    output:
        directory=directory("RNAseqCNV_output/gatk/{sample_id}_gatk/"),
        rna_seq_cnv_alteration_file="RNAseqCNV_output/gatk/{sample_id}_gatk/estimation_table.tsv",
        rna_seq_cnv_log2foldchange_file="RNAseqCNV_output/gatk/{sample_id}_gatk/log2_fold_change_per_arm.tsv",
        rna_seq_cnv_manual_an_table_file="RNAseqCNV_output/gatk/{sample_id}_gatk/manual_an_table.tsv",
        rna_seq_cnv_plot="RNAseqCNV_output/gatk/{sample_id}_gatk/{sample_id}/{sample_id}_CNV_main_fig.png"


    benchmark:
        "benchmarks/{sample}.rnaseq_cnv.benchmark.txt"

    conda:
        "envs/rnaseqenv.yaml"

    shell:
        "Rscript {input.r_script} {input.config_file} {input.metadata_file} {wildcards.sample_id} {output.directory};"



rule check_subtype_and_karyotype:
    input:
        prediction_file = "allcatch_output/{sample}/predictions.tsv",
        rna_seq_cnv_estimation_file="RNAseqCNV_output/gatk/{sample}_gatk/estimation_table.tsv",
        fusioncatcher_file = "fusioncatcher_output/{sample}/final-list_candidate-fusion-genes.txt",
        arriba_file= "fusions/{sample}.tsv",
        chromosome_counts_karyotype_file= "data/annotation/chromosome_counts_vs_subtype.txt",
        anno_gene_fusions_file= "data/annotation/anno_Gene_fusions_vs_Subtypes.txt",
        r_script= "scripts/define_subtype_and_caryotype.R"

    output:
        csv= "comparison/{sample}.csv"
        
    conda:
        "envs/subtype.yaml"

    shell:
        "Rscript {input.r_script} {input.prediction_file} {input.rna_seq_cnv_estimation_file} {input.fusioncatcher_file} {input.arriba_file} {input.chromosome_counts_karyotype_file} {input.anno_gene_fusions_file} {output.csv} "



rule aggregate_output:
    input:
            prediction_file = "allcatch_output/{sample}/predictions.tsv",
            fusioncatcher_file = "fusioncatcher_output/{sample}/final-list_candidate-fusion-genes.txt",
            arriba_file = "fusions/{sample}.tsv",
            rna_seq_cnv_log2foldchange_file = "RNAseqCNV_output/gatk/{sample}_gatk/log2_fold_change_per_arm.tsv",
            rna_seq_cnv_manual_an_table_file = "RNAseqCNV_output/gatk/{sample}_gatk/manual_an_table.tsv",
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
        rna_seq_cnv_log2foldchange_file="RNAseqCNV_output/gatk/{sample}_gatk/log2_fold_change_per_arm.tsv",
        rna_seq_cnv_plot="RNAseqCNV_output/gatk/{sample}_gatk/{sample}/{sample}_CNV_main_fig.png",
        rna_seq_cnv_manual_an_table_file="RNAseqCNV_output/gatk/{sample}_gatk/manual_an_table.tsv",
        star_log_final_out_file="STAR_output/{sample}/Log.final.out",
        multiqc_fqc_right="multiqc/{sample}_right/multiqc_report.html",
        multiqc_fqc_left="multiqc/{sample}_left/multiqc_report.html",
        comparison_file="comparison/{sample}.csv",
        hotspots="Hotspots/{sample}"

    output:
        html="interactive_output/{sample}/output_report_{sample}.html"

    shell:
        """
        mkdir -p interactive_output/{wildcards.sample}/fusions &&
        cp {input.arriba_file_fusion} interactive_output/{wildcards.sample}/fusions &&

        mkdir -p interactive_output/{wildcards.sample}/multiqc_right &&
        cp {input.multiqc_fqc_right} interactive_output/{wildcards.sample}/multiqc_right &&

        mkdir -p interactive_output/{wildcards.sample}/multiqc_left &&
        cp {input.multiqc_fqc_left} interactive_output/{wildcards.sample}/multiqc_left &&

        mkdir -p interactive_output/{wildcards.sample}/RNAseqCNV &&
        cp {input.rna_seq_cnv_plot} interactive_output/{wildcards.sample}/RNAseqCNV &&

        python scripts/generate_report.py  {input.prediction_file} {input.fusioncatcher_file} {input.arriba_file} {input.rna_seq_cnv_log2foldchange_file} {input.rna_seq_cnv_manual_an_table_file} {input.star_log_final_out_file}  {input.comparison_file} {input.hotspots} {wildcards.sample} {output.html}
        """
        
