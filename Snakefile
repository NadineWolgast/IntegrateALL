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
creatingfolders('refs/fusioncatcher')
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
        expand("Variants_RNA_Seq_Reads/{sample_id}/fixed-rg/{sample_id}.bam", sample_id=list(samples.keys())),
        expand("Variants_RNA_Seq_Reads/{sample_id}/deduped_bam/{sample_id}.bam.bai", sample_id=list(samples.keys())),
        expand("Variants_RNA_Seq_Reads/{sample_id}/split/{sample_id}.bam", sample_id=list(samples.keys())),
        expand("Variants_RNA_Seq_Reads/{sample_id}/recal/{sample_id}_recal.table", sample_id=list(samples.keys())),
        expand("Variants_RNA_Seq_Reads/{sample_id}/recal/{sample_id}.bam", sample_id=list(samples.keys())),
        expand("Variants_RNA_Seq_Reads/{sample_id}/filter/{sample_id}.snvs.filtered.vcf", sample_id=list(samples.keys())),
        expand("fusions/{sample_id}.pdf",sample_id=samples.keys()),
        expand("fusions/{sample_id}.tsv",sample_id=samples.keys()),

        expand("multiqc/{sample}/multiqc_data/multiqc_fastqc.txt", sample=fastq_dataframe['sample_id']),
        expand("fusioncatcher_output/{sample_id}/final-list_candidate-fusion-genes.txt",sample_id=list(samples.keys())),
        expand("data/vcf_files/GATK/{sample_id}_Gatk.tsv", sample_id=samples.keys()),
        expand("data/tpm/{sample_id}.tsv", sample_id=list(samples.keys())),
        #expand("data/cpm/{sample_id}.tsv", sample_id=list(samples.keys())),        
        #expand("pysamstats_output_dir/{sample_id}/", sample_id=list(samples.keys())),
        expand("Hotspots/{sample_id}",sample_id=list(samples.keys())),
        #expand("comparison/{sample_id}.csv", sample_id= samples.keys()),
        expand("data/single_counts/{sample_id}.txt", sample_id=samples.keys()),
        expand("allcatch_output/{sample_id}/predictions.tsv", sample_id= samples.keys()),
        #expand("aggregated_output/{sample}.csv", sample=list(samples.keys())),
        expand("RNAseqCNV_output/gatk/{sample_id}_gatk", sample_id=samples.keys()),
        expand("Final_classification/{sample_id}_output_report.csv",sample_id=list(samples.keys())),
        expand("interactive_output/{sample}/output_report_{sample}.html",  sample=list(samples.keys()))



rule check_samples:
    input:
        samples_csv=config["sample_file"]
    output:
        "check_samples.txt"
    script:
        "scripts/check_samples.py"




# Function to get input FASTQ files based on sample_id
def get_input_fastqs(wildcards):
    sample_name = wildcards.sample
    fastq_file = fastq_dataframe[fastq_dataframe['sample_id'] == sample_name]['FASTQ'].values[0]
    return fastq_file



rule install_all:
    input: 
        # Reference files (~17GB total)
        "refs/GATK/GRCH38/dbSNP.vcf",                                            # ~13GB
        "refs/GATK/GRCH38/Homo_sapiens.GRCh38.dna.primary_assembly.fa",         # ~3GB
        "refs/GATK/GRCH38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai",     # ~3KB
        absolute_path + "/refs/GATK/STAR/ensembl_94_100/SA",                     # ~1GB
        "scripts/dbSNP_hg38.rda",                                                # ~50MB
        "scripts/pseudoautosomal_regions_hg38.rda",                              # ~1KB
        # Tool installations (via marker files)
        "logs/install_arriba_draw_fusions.done",                                 # ~10MB
        "logs/install_allcatchr.done",                                           # R packages
        "logs/install_rnaseq_cnv.done",                                          # R packages
        # Large downloads (optional - comment out for testing)
        absolute_path + "/refs/fusioncatcher/fusioncatcher-master/data/human_v102/version.txt"  # ~4.4GB
    message: "All reference files and tools installed successfully! Total size: ~21GB"
    resources:
        # Allow parallel downloads/installations
        mem_mb=2000,
        # Conservative retry logic for network downloads
        tmpdir="/tmp"
    retries: 3

rule download_ref:
    input:
        star_directory = absolute_path + "/refs/GATK"
    output:
        vcf= "refs/GATK/GRCH38/dbSNP.vcf",
        ref= "refs/GATK/GRCH38/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    message: "Downloading reference files (~16GB): FASTA + dbSNP VCF"
    resources:
        mem_mb=1000,
        tmpdir="/tmp"
    retries: 3
    shell:
        """
        mkdir -p {input.star_directory} &&
        cd {input.star_directory} && 
        # Robust download with retry logic and progress bar
        wget --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 3 \
             --progress=bar --show-progress \
             'http://141.2.194.197/rnaeditor_annotations/GRCH38.tar.gz' && 
        # Verify download success before extraction
        if [ -f "GRCH38.tar.gz" ]; then
            echo "✅ Download successful, extracting..."
            tar -xzf GRCH38.tar.gz &&
            echo "✅ Reference files extracted successfully"
        else
            echo "❌ Download failed" && exit 1
        fi
        """


rule download_rda:
    output:
        dbsnp="scripts/dbSNP_hg38.rda",
        par="scripts/pseudoautosomal_regions_hg38.rda"
    message: "Downloading RNAseqCNV reference data (~50MB)"
    resources:
        mem_mb=500
    retries: 3
    shell:
        """
        # Download dbSNP with retry logic
        wget --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 3 \
             --progress=bar --show-progress \
             -O {output.dbsnp} https://github.com/honzee/RNAseqCNV/raw/master/data/dbSNP_hg38.rda &&
        echo "✅ dbSNP_hg38.rda downloaded successfully" &&
        
        # Download pseudoautosomal regions with retry logic
        wget --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 3 \
             --progress=bar --show-progress \
             -O {output.par} https://github.com/honzee/RNAseqCNV/raw/master/data/pseudoautosomal_regions_hg38.rda &&
        echo "✅ pseudoautosomal_regions_hg38.rda downloaded successfully"
        """


rule index_ref:
    input: "refs/GATK/GRCH38/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    output: "refs/GATK/GRCH38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai"
    shell:
        "samtools faidx {input}"


rule index_star:
        input:
            fa = absolute_path + "/refs/GATK/GRCH38/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
            gtf = absolute_path + "/refs/GATK/GRCH38/Homo_sapiens.GRCh38.83.gtf"
        output:
            directory = directory(absolute_path + "/refs/GATK/STAR/ensembl_94_100")

        threads: config['threads']

        conda:
            "envs/star.yaml"

        shell:
            'mkdir -p {output.directory} && '
            'STAR --runThreadN {config[threads]} '
            '--runMode genomeGenerate '
            '--genomeDir {output} '
            '--genomeFastaFiles {input.fa} '
            '--sjdbGTFfile {input.gtf} '
            '--sjdbOverhang 100'



rule install_arriba_draw_fusions:
    output:
        touch("logs/install_arriba_draw_fusions.done")
    message: "Installing Arriba draw_fusions tool (~10MB)"
    benchmark:
        "benchmarks/install_arriba_draw_fusions.benchmark.txt"
    resources:
        mem_mb=1000
    retries: 3
    shell:
        '''
        # Robust download with retry logic and verification
        wget --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 3 \
             --progress=bar --show-progress \
             'https://github.com/suhrig/arriba/releases/download/v2.4.0/arriba_v2.4.0.tar.gz' &&
        
        # Verify download success before extraction
        if [ -f "arriba_v2.4.0.tar.gz" ]; then
            echo "✅ Arriba download successful, extracting..."
            tar -xzf arriba_v2.4.0.tar.gz &&
            echo "✅ Arriba installation completed successfully"
        else
            echo "❌ Arriba download failed" && exit 1
        fi &&
        
        touch {output}
        '''


rule install_allcatchr:
    output:
        touch("logs/install_allcatchr.done")
    message: "Installing ALLCatchR R package from GitHub"
    conda:
        "envs/install_allcatchr.yaml"
    benchmark:
        "benchmarks/install_allcatchr.benchmark.txt"
    resources:
        mem_mb=4000  # R package compilation needs memory
    retries: 2  # GitHub API can be flaky
    shell:
        """
        echo "📦 Installing ALLCatchR package..."
        Rscript -e '
            tryCatch({{
                devtools::install_github("ThomasBeder/ALLCatchR_bcrabl1", 
                                       Ncpus = {config[threads]},
                                       force = TRUE,
                                       upgrade = "never")
                cat("✅ ALLCatchR installation completed successfully\\n")
            }}, error = function(e) {{
                cat("❌ ALLCatchR installation failed:", conditionMessage(e), "\\n")
                quit(status = 1)
            }})'
        touch {output}
        """


rule install_rnaseq_cnv:
    output:
        touch("logs/install_rnaseq_cnv.done")
    message: "Installing RNAseqCNV R package from GitHub"
    conda:
        "envs/rnaseqenv.yaml"
    benchmark:
        "benchmarks/install_rnaseq_cnv.benchmark.txt"
    resources:
        mem_mb=4000  # R package compilation needs memory
    retries: 2  # GitHub API can be flaky
    shell:
        """
        echo "📦 Installing RNAseqCNV package with dependencies..."
        Rscript -e '
            tryCatch({{
                devtools::install_github("honzee/RNAseqCNV", 
                                       dependencies = TRUE, 
                                       Ncpus = {config[threads]},
                                       force = TRUE,
                                       upgrade = "never")
                cat("✅ RNAseqCNV installation completed successfully\\n")
            }}, error = function(e) {{
                cat("❌ RNAseqCNV installation failed:", conditionMessage(e), "\\n")
                quit(status = 1)
            }})'
        touch {output}
        """
        

rule install_fusioncatcher:
    input:
        data_directory=absolute_path + "/refs/fusioncatcher"
    output:
        directory(absolute_path + "/refs/fusioncatcher/fusioncatcher-master/data/human_v102")
    message: "Installing FusionCatcher database - checking for existing references first"
    benchmark:
        "benchmarks/install_fusioncatcher.benchmark.txt"
    resources:
        mem_mb=2000,
        tmpdir="/tmp"
    retries: 2
    shell:
        """
        cd {input.data_directory} &&
        
        # Check if we can copy from existing installation first
        EXISTING_DB="/media/nadine/InternalMaybe/Blast-o-Matic-Fusioninator_cluster/refs/fusioncatcher/fusioncatcher-master/data/human_v102"
        
        if [ -d "$EXISTING_DB" ] && [ -f "$EXISTING_DB/version.txt" ]; then
            echo "🔍 Found existing FusionCatcher database, copying instead of downloading..."
            echo "📁 Source: $EXISTING_DB"
            
            # Download FusionCatcher source code first
            if [ ! -f "master.zip" ]; then
                echo "📦 Downloading FusionCatcher source code..."
                wget --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 3 \
                     --progress=bar --show-progress \
                     'https://github.com/ndaniel/fusioncatcher/archive/refs/heads/master.zip' && 
                unzip -q master.zip
            fi
            
            # Copy existing database (much faster than download)
            echo "📋 Copying FusionCatcher database (~4.4GB)..."
            mkdir -p fusioncatcher-master/data
            cp -r "$EXISTING_DB" fusioncatcher-master/data/ &&
            echo "✅ FusionCatcher database copied successfully! (Saved 30+ minutes download time)"
            
        else
            echo "❌ No existing database found at $EXISTING_DB"
            echo "📦 Falling back to download method..."
            
            # Download FusionCatcher source code with robust retry logic
            echo "📦 Downloading FusionCatcher source code..."
            wget --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 3 \
                 --progress=bar --show-progress \
                 'https://github.com/ndaniel/fusioncatcher/archive/refs/heads/master.zip' && 
            
            # Verify download before extraction
            if [ -f "master.zip" ]; then
                echo "✅ FusionCatcher source download successful, extracting..."
                unzip -q master.zip &&
                echo "✅ FusionCatcher source extracted successfully"
            else
                echo "❌ FusionCatcher source download failed" && exit 1
            fi &&
            
            # Download the large human database (~4.4GB) 
            cd fusioncatcher-master/data && 
            echo "📦 Starting FusionCatcher human database download (~4.4GB)..."
            echo "⏰ This will take 30+ minutes depending on your internet connection"
            echo "💡 Speed issues? Consider copying from existing installation next time!"
            ./download-human-db.sh &&
            echo "✅ FusionCatcher database installation completed successfully!"
        fi
        """


rule fastqc:
    input:
        get_input_fastqs
    output:
        html="fastqc/{sample}.html",
        zip="fastqc/{sample}.zip"
    log:
        "logs/fastqc/{sample}/fastqc.log"

    resources:
        mem_mb=1000

    benchmark:
        "benchmarks/{sample}.fastqc.benchmark.txt"

    wrapper:
        "v3.10.2/bio/fastqc"

rule unzip:
    input:
        sample="fastqc/{sample}.zip"
    output:
        dir=directory("fastqc/{sample}"),
	    sentinel="fastqc/{sample}/unzip_done.sentinel"
    resources:
        mem_mb=1000
    shell:
        "unzip -q {input.sample} -d {output.dir} && touch {output.sentinel}"


rule multiqc_file:
    input:
        "fastqc/{sample}/"
    output:
        multiqc_report="multiqc/{sample}/multiqc_report.html",
        multiqc_fqc= "multiqc/{sample}/multiqc_data/multiqc_fastqc.txt"

    params:
        output_dir= "multiqc/{sample}/"

    resources:
        mem_mb=1000

    benchmark:
        "benchmarks/{sample}.multiqc.benchmark.txt"

    shell:
        "multiqc {input} -o {params.output_dir} --force"


rule run_star_aligner:
    input:
        fastq1 = lambda wildcards: samples[wildcards.sample_id][0],  # Path to left FASTQ from sample sheet
        fastq2 = lambda wildcards: samples[wildcards.sample_id][1],  # Path to right FASTQ from sample sheet
        genome_index = absolute_path + "/refs/GATK/STAR/ensembl_94_100"  

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
        mem_mb=config['star_mem']

    shell:
        'mkdir -p {output.directory} && '
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
        "STAR_output/{sample_id}/Aligned.sortedByCoord.out.bam"
    output:
        "STAR_output/{sample_id}/Aligned.sortedByCoord.out.bam.bai"
    log:
        "logs/samtools_index/{sample_id}.log"
    params:
        extra=""  # optional params string
    benchmark:
        "benchmarks/{sample_id}.samtools_index.benchmark.txt"
    threads: config['threads']
    wrapper:
        "v2.6.0/bio/samtools/index"


rule pysamstat:
    input:
        bam="STAR_output/{sample_id}/Aligned.sortedByCoord.out.bam",
        bai="STAR_output/{sample_id}/Aligned.sortedByCoord.out.bam.bai",
        fa= absolute_path + "/refs/GATK/GRCH38/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

    conda:
        "envs/pysamstat.yaml"
        
    params:
        out_dir="pysamstats_output_dir/{sample_id}/"

    output:
        pysamstats_output_dir = directory("pysamstats_output_dir/{sample_id}/"),
        ikzf1="pysamstats_output_dir/{sample_id}/{sample_id}_IKZF1.csv"

    threads: config['threads']

    resources:
        threads=config['threads'],
        mem_mb=20000

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

    resources:
        mem_mb=10000    
       
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
        extra=""                         # Optional parameters
    threads: config['threads']
    resources:
        mem_mb=20000
    wrapper:
        "v7.2.0/bio/arriba"


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

    resources:
        mem_mb=20000

    shell:
        '''
        Rscript {input.r_script} \
        --fusions={input.fusions} \
        --alignments={input.bam} \
        --output={output.pdf} \
        --annotation={input.annotation}\
        --cytobands=arriba_v2.4.0/database/cytobands_hg38_GRCh38_v2.4.0.tsv \
        --proteinDomains=arriba_v2.4.0/database/protein_domains_hg38_GRCh38_v2.4.0.gff3
        '''



rule run_fusioncatcher:
    input:
        data_directory= absolute_path + "/refs/fusioncatcher/fusioncatcher-master/data/human_v102",
        left= lambda wildcards: samples[wildcards.sample_id][0],
        right= lambda wildcards: samples[wildcards.sample_id][1]

    output:
        dir=directory("fusioncatcher_output/{sample_id}"),
        keep_file="fusioncatcher_output/{sample_id}/final-list_candidate-fusion-genes.txt",
        sentinel="fusioncatcher_output/{sample_id}/fusioncatcher_done.sentinel"

    log:
        "logs/fusioncatcher/{sample_id}.log"

    conda:
        "envs/fusioncatcher.yaml"

    params:
        sample_id=lambda wildcards: wildcards.sample_id

    benchmark:
        "benchmarks/{sample_id}.fusioncatcher.benchmark.txt"
        
    resources:
        threads=config['threads'],
        mem_mb=config['star_mem']

    shell:
        '''
        mkdir -p logs/fusioncatcher
        fusioncatcher \
        -d {input.data_directory} \
        -i {input.left},{input.right} \
        -o {output.dir} \
        --skip-blat \
        --threads={config[threads]} \
        > logs/fusioncatcher/{wildcards.sample_id}.log 2>&1 \
        && touch {output.sentinel}
        '''


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
        input_file = 'data/counts/{sample_id}.tsv'

    output:
        out_file = "allcatch_output/{sample_id}/predictions.tsv"

    conda:
        "envs/install_allcatchr.yaml"

    params:
        out_dir = absolute_path + "/allcatch_output/{sample_id}",
        abs = absolute_path

    threads: config['threads']

    resources:
        threads=config['threads'],
        mem_mb=8000    
        
    benchmark:
        "benchmarks/{sample_id}.allcatchr.benchmark.txt"

    shell:
        "mkdir -p {params.out_dir} && "
        "cd {params.out_dir} && "
        "Rscript {params.abs}/{input.r_script} {params.abs}/{input.input_file} {params.abs}/{output.out_file} "
        

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

    log:
        "logs/picard/replace_rg/{sample}.log",
    params:
        extra="--RGLB lib1 --RGPL illumina --RGPU {sample} --RGSM {sample}",
        java_opts=""

    resources:
        mem_mb=2048,
    wrapper:
        "v3.10.2/bio/picard/addorreplacereadgroups"


rule markduplicates_bam:
    input:
        bams="Variants_RNA_Seq_Reads/{sample}/fixed-rg/{sample}.bam"

    output:
        bam=temporary("Variants_RNA_Seq_Reads/{sample}/deduped_bam/{sample}.bam"),
        metrics=temporary("Variants_RNA_Seq_Reads/{sample}/deduped_bam/{sample}.metrics.txt")
    log:
        "logs/picard/dedup_bam/{sample}.log"

    resources:
        mem_mb=5000
    wrapper:
        "v3.10.2/bio/picard/markduplicates"


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
        "v3.10.2/bio/samtools/index"


rule splitncigarreads:
    input:
        bam="Variants_RNA_Seq_Reads/{sample}/deduped_bam/{sample}.bam",
        bai="Variants_RNA_Seq_Reads/{sample}/deduped_bam/{sample}.bam.bai",
        ref= "refs/GATK/GRCH38/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    output:
        temporary("Variants_RNA_Seq_Reads/{sample}/split/{sample}.bam"),
    log:
        "logs/gatk/splitNCIGARreads/{sample}.log",
    params:
        extra="",  # optional
        java_opts="",  # optional
    resources:
        mem_mb=10000,
    threads: config['threads']
    wrapper:
        "v4.3.0/bio/gatk/splitncigarreads"




rule gatk_baserecalibrator:
    input:
        bam="Variants_RNA_Seq_Reads/{sample}/split/{sample}.bam",
        ref= "refs/GATK/GRCH38/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
        dict= "refs/GATK/GRCH38/Homo_sapiens.GRCh38.dna.primary_assembly.dict",
        known= "refs/GATK/GRCH38/dbSNP.vcf",
    output:
        recal_table=temporary("Variants_RNA_Seq_Reads/{sample}/recal/{sample}_recal.table"),
    log:
        "logs/gatk/baserecalibrator/{sample}.log",
    params:
        extra="",  # optional
        java_opts="",  # optional
    resources:
        mem_mb=10000,
    threads: config['threads']
    wrapper:
        "v4.3.0/bio/gatk/baserecalibrator"


rule gatk_applybqsr:
    input:
        bam="Variants_RNA_Seq_Reads/{sample}/split/{sample}.bam",
        ref="refs/GATK/GRCH38/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
        dict="refs/GATK/GRCH38/Homo_sapiens.GRCh38.dna.primary_assembly.dict",
        recal_table="Variants_RNA_Seq_Reads/{sample}/recal/{sample}_recal.table"
    output:
        bam=temporary("Variants_RNA_Seq_Reads/{sample}/recal/{sample}.bam"),
    log:
        "logs/gatk/gatk_applybqsr/{sample}.log",
    params:
        extra="",  # optional
        java_opts="",  # optional
        embed_ref=True,  # embed the reference in cram output
    resources:
        mem_mb=10000,
    threads: config['threads']
    wrapper:
        "v4.3.0/bio/gatk/applybqsr"

rule haplotype_caller:
    input:
        bam= "Variants_RNA_Seq_Reads/{sample}/recal/{sample}.bam",
        ref= "refs/GATK/GRCH38/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
        #known= "refs/STAR/dbSNP.vcf" #optional

    output:
        vcf=temporary("Variants_RNA_Seq_Reads/{sample}/calls/{sample}.vcf")

    log:
        "logs/gatk/haplotypecaller/{sample}.log"

    params:
        extra="",  # optional
        java_opts="",  # optional

    threads: config['threads']

    resources:
        mem_mb=10000

    wrapper:
        "v4.3.0/bio/gatk/haplotypecaller"

rule gatk_filter:
    input:
        vcf="Variants_RNA_Seq_Reads/{sample}/calls/{sample}.vcf",
        ref= "refs/GATK/GRCH38/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

    output:
        vcf="Variants_RNA_Seq_Reads/{sample}/filter/{sample}.snvs.filtered.vcf"

    log:
        "logs/gatk/filter/{sample}.snvs.log"

    params:
        filters={"myfilter": "FS > 30.0 || QD < 2.0"},
        extra="",  # optional arguments, see GATK docs
        java_opts="",  # optional

    resources:
        mem_mb=10000

    threads: config['threads']

    wrapper:
        "v4.3.0/bio/gatk/variantfiltration"


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
        r_script="scripts/modified_RNASeqCNV_wrapper.R",
        config_file="data/config.txt",
        metadata_file="data/meta.txt",
        input_counts="data/single_counts/{sample_id}.txt",
        input_tsv="data/vcf_files/GATK/{sample_id}_Gatk.tsv"

    output:
        directory=directory("RNAseqCNV_output/gatk/{sample_id}_gatk/"),
        rna_seq_cnv_alteration_file="RNAseqCNV_output/gatk/{sample_id}_gatk/estimation_table.tsv",
        rna_seq_cnv_log2foldchange_file="RNAseqCNV_output/gatk/{sample_id}_gatk/log2_fold_change_per_arm.tsv",
        rna_seq_cnv_manual_an_table_file="RNAseqCNV_output/gatk/{sample_id}_gatk/manual_an_table.tsv",
        rna_seq_cnv_plot="RNAseqCNV_output/gatk/{sample_id}_gatk/{sample_id}/{sample_id}_CNV_main_fig.png",
        ml_input="RNAseqCNV_output/gatk/{sample_id}_gatk/{sample_id}_ml_input.csv",
    conda:
        "envs/rnaseqenv.yaml"
        
    resources:
        threads=config['threads'],
        mem_mb=5000
        
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


rule predict_karyotype:
    input:
        ml_input="RNAseqCNV_output/gatk/{sample}_gatk/{sample}_ml_input.csv",

    output:
        csv="karyotype_prediction/{sample}.csv"

    conda:
        "envs/ml.yaml"

    shell:
        "python scripts/predict_karyotype.py {input.ml_input} {output.csv} {wildcards.sample} "
        

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

rule final_classification:
    input:
        allcatchr_file="allcatch_output/{sample}/predictions.tsv",
        karyotype="karyotype_prediction/{sample}.csv",
        fusioncatcher_file="fusioncatcher_output/{sample}/final-list_candidate-fusion-genes.txt",
        arriba_file="fusions/{sample}.tsv",
        hotspots="Hotspots/{sample}",
        classification_file="data/annotation/Class_test.csv"

    output:
        csv="Final_classification/{sample}_output_report.csv",
        text="Final_classification/{sample}_output_txt.csv",
        driver="Final_classification/{sample}_driver.csv",
        curation="Final_classification/{sample}_curation.csv"

    shell:
        """
        python scripts/make_final_classification.py  {wildcards.sample} {input.allcatchr_file} {input.karyotype} {input.fusioncatcher_file} {input.arriba_file} {input.hotspots} {input.classification_file} {output.csv} {output.text} {output.curation} {output.driver}
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
        hotspots="Hotspots/{sample}",
        karyotype="karyotype_prediction/{sample}.csv",
        text="Final_classification/{sample}_output_txt.csv",
        driver="Final_classification/{sample}_driver.csv"

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
        cp scripts/logo.png interactive_output/{wildcards.sample}/ &&

        python scripts/generate_report.py  {input.prediction_file} {input.fusioncatcher_file} {input.arriba_file} {input.rna_seq_cnv_log2foldchange_file} {input.rna_seq_cnv_manual_an_table_file} {input.star_log_final_out_file}  {input.comparison_file} {input.hotspots} {wildcards.sample} {input.karyotype} {input.text} {input.driver} {output.html}
        """
