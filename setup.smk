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

# Setup target - installs everything needed for the pipeline
rule setup_all:
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
    message: "Indexing reference genome"
    shell:
        "samtools faidx {input}"

rule index_star:
        input:
            fa = absolute_path + "/refs/GATK/GRCH38/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
            gtf = absolute_path + "/refs/GATK/GRCH38/Homo_sapiens.GRCh38.83.gtf"
        output:
            directory = directory(absolute_path + "/refs/GATK/STAR/ensembl_94_100")
        message: "Building STAR genome index (~1GB)"
        threads: config['threads']
        conda:
            "envs/star.yaml"
        resources:
            mem_mb=config['star_mem'],
            tmpdir="/tmp"
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
    message: "Installing FusionCatcher database - checking if already present"
    benchmark:
        "benchmarks/install_fusioncatcher.benchmark.txt"
    resources:
        mem_mb=2000,
        tmpdir="/tmp"
    retries: 2
    shell:
        """
        cd {input.data_directory} &&
        
        # Check if database already exists locally
        if [ -d "fusioncatcher-master/data/human_v102" ] && [ -f "fusioncatcher-master/data/human_v102/version.txt" ]; then
            echo "✅ FusionCatcher database already exists locally - skipping installation"
            exit 0
        fi
        
        echo "📦 FusionCatcher database not found locally"
        echo "💡 Tip: You can manually copy existing databases to speed up installation"
        echo "📦 Proceeding with download and installation..."
        
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
        ./download-human-db.sh &&
        echo "✅ FusionCatcher database installation completed successfully!"
        """