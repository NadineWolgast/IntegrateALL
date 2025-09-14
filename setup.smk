"""
IntegrateALL Setup Workflow
===========================

One-time setup for reference files, databases, and tool installations.
Run this ONCE before using the main analysis pipeline.

Usage:
    snakemake --snakefile setup.smk --cores <N>

This will install:
- Reference genome and annotations (~16GB)
- STAR genome index (~1GB) 
- RNAseqCNV reference data (~50MB)
- FusionCatcher database (~4.4GB)
- R packages (ALLCatchR, RNAseqCNV)
- Arriba draw_fusions tool (~10MB)

Total: ~21GB of reference data
"""

configfile: "config.yaml"

# Import absolute path from config
absolute_path = config["absolute_path"]

# ------------------------
# SETUP TARGET
# ------------------------
rule setup_all:
    input:
        # Reference files (~17GB total)
        "refs/GATK/GRCH38/dbSNP.vcf",
        "refs/GATK/GRCH38/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
        "refs/GATK/GRCH38/Homo_sapiens.GRCh38.83.gtf",
        "refs/GATK/GRCH38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai",
        # STAR index (match rule 'index_star' output -> directory)
        absolute_path + "/refs/GATK/STAR/ensembl_94_100",
        # RNAseqCNV RDA files
        "scripts/dbSNP_hg38.rda",
        "scripts/pseudoautosomal_regions_hg38.rda",
        # Tool installs (marker files)
        "logs/install_arriba_draw_fusions.done",
        "logs/install_allcatchr.done",
        "logs/install_rnaseq_cnv.done",
        # FusionCatcher database (match rule 'install_fusioncatcher' output -> directory)
        absolute_path + "/refs/fusioncatcher/fusioncatcher-master/data/human_v102",
        #"logs/install_fusioncatcher.done"
    message: "🎉 All reference files and tools installed successfully! Total size: ~21GB"
    shell:
        """
        echo "✅ IntegrateALL setup completed successfully!"
        echo "📊 Reference data installed: ~21GB"
        echo "🚀 Ready to run main analysis pipeline"
        echo ""
        echo "Next steps:"
        echo "  1. Check your samples.csv file"
        echo "  2. Run: snakemake --cores <N>"
        """

# ------------------------
# DOWNLOADS
# ------------------------
rule download_ref:
    output:
        vcf = protected("refs/GATK/GRCH38/dbSNP.vcf"),
        ref = protected("refs/GATK/GRCH38/Homo_sapiens.GRCh38.dna.primary_assembly.fa"),
        gtf = protected("refs/GATK/GRCH38/Homo_sapiens.GRCh38.83.gtf")
    message: "Downloading reference files (~16GB): FASTA + dbSNP VCF - PROTECTED from deletion"
    resources:
        mem_mb = 1000,
        tmpdir = "/tmp"
    retries: 3
    shell:
        """
        mkdir -p refs/GATK &&
        cd refs/GATK &&
        wget --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 3 \
             --progress=bar --show-progress \
             'http://141.2.194.197/rnaeditor_annotations/GRCH38.tar.gz' &&
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
        dbsnp = protected("scripts/dbSNP_hg38.rda"),
        par = protected("scripts/pseudoautosomal_regions_hg38.rda")
    message: "Downloading RNAseqCNV reference data (~50MB) - PROTECTED from deletion"
    resources:
        mem_mb = 500
    retries: 3
    shell:
        """
        wget --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 3 \
             --progress=bar --show-progress \
             -O {output.dbsnp} https://github.com/honzee/RNAseqCNV/raw/master/data/dbSNP_hg38.rda &&
        echo "✅ dbSNP_hg38.rda downloaded successfully" &&
        wget --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 3 \
             --progress=bar --show-progress \
             -O {output.par} https://github.com/honzee/RNAseqCNV/raw/master/data/pseudoautosomal_regions_hg38.rda &&
        echo "✅ pseudoautosomal_regions_hg38.rda downloaded successfully"
        """

# ------------------------
# INDEXING
# ------------------------
rule index_ref:
    input: "refs/GATK/GRCH38/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    output: protected("refs/GATK/GRCH38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai")
    message: "Indexing reference genome - PROTECTED from deletion"
    shell:
        "samtools faidx {input}"

rule index_star:
    input:
        fa  = "refs/GATK/GRCH38/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
        gtf = "refs/GATK/GRCH38/Homo_sapiens.GRCh38.83.gtf"
    output:
        directory = protected(directory(absolute_path + "/refs/GATK/STAR/ensembl_94_100"))
    message: "Building STAR genome index (~1GB) - PROTECTED from deletion"
    threads: config['threads']
    conda:
        "envs/star.yaml"
    resources:
        mem_mb = config['star_mem'],
        tmpdir = "/tmp"
    shell:
        """
        mkdir -p {output.directory} &&
        STAR --runThreadN {config[threads]} \
             --runMode genomeGenerate \
             --genomeDir {output.directory} \
             --genomeFastaFiles {input.fa} \
             --sjdbGTFfile {input.gtf} \
             --sjdbOverhang 100
        """

# ------------------------
# TOOL INSTALLS
# ------------------------
rule install_arriba_draw_fusions:
    output:
        touch("logs/install_arriba_draw_fusions.done")
    message: "Installing Arriba draw_fusions tool (~10MB)"
    benchmark:
        "benchmarks/install_arriba_draw_fusions.benchmark.txt"
    resources:
        mem_mb = 1000
    retries: 3
    shell:
        '''
        wget --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 3 \
             --progress=bar --show-progress \
             'https://github.com/suhrig/arriba/releases/download/v2.4.0/arriba_v2.4.0.tar.gz' &&
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
        mem_mb = 4000
    retries: 2
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
        mem_mb = 4000
    retries: 2
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
        data_directory = absolute_path + "/refs/fusioncatcher"
    output:
        directory_output = protected(directory(absolute_path + "/refs/fusioncatcher/fusioncatcher-master/data/human_v102")),
        done_marker = "logs/install_fusioncatcher.done"
    message: "Installing FusionCatcher database - PROTECTED from deletion"
    benchmark:
        "benchmarks/install_fusioncatcher.benchmark.txt"
    resources:
        mem_mb = 2000,
        tmpdir = "/tmp"
    retries: 2
    shell:
        """
        # Create logs directory first
        cd {absolute_path}
        mkdir -p logs
        
        # Check if FusionCatcher database already exists
        if [ -d "{output.directory_output}" ] && [ -f "{output.directory_output}/version.txt" ]; then
            echo "✅ FusionCatcher database already exists at {output.directory_output}"
            echo "✅ Skipping download - creating done marker"
            touch {output.done_marker}
            exit 0
        fi

        # If not present, proceed with installation
        cd {input.data_directory} || exit 1

        echo "📦 FusionCatcher database not found locally"
        echo "📦 Downloading FusionCatcher source code..."
        wget --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 3 \
             --progress=bar --show-progress \
             'https://github.com/ndaniel/fusioncatcher/archive/refs/heads/master.zip' &&
        if [ -f "master.zip" ]; then
            echo "✅ FusionCatcher source download successful, extracting..."
            unzip -q master.zip &&
            echo "✅ FusionCatcher source extracted successfully"
        else
            echo "❌ FusionCatcher source download failed" && exit 1
        fi

        cd fusioncatcher-master/data &&
        echo "📦 Starting FusionCatcher human database download (~4.4GB)..."
        ./download-human-db.sh &&
        echo "✅ FusionCatcher database installation completed successfully!" &&
        cd {absolute_path} &&
        touch {output.done_marker}
        """
