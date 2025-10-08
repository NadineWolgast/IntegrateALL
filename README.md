# IntegrateALL

IntegrateALL is a machine learning based pipeline for multilevel-data extraction and subtype classification in B-cell precursor ALL. It provides the user with an interactive report, containing subtype prediction, virtual karyotype, fusions, quality control, single nucleotide variants using a snakemake workflow management system.

## Features

- **Machine Learning Classification**: Automated subtype prediction for B-cell precursor ALL using ALLCatchR and automated karyotype prediction using KaryALL
- **Comprehensive Analysis**: Fusion detection, CNV analysis, variant calling, and quality control
- **Interactive Reporting**: HTML reports with visualizations and predictions
- **Two-Workflow Architecture**: Separate setup and analysis workflows for improved stability
- **Automatic Classification**: Classification according to WHO-HAEM5 and ICC classification

## Overview

![Pipeline Overview](Pipeline_Overview.png?raw=true)

## Pipeline Architecture

![Rule Flowchart](pipeline_rule_flowchart.png?raw=true)

IntegrateALL uses a **two-workflow architecture** with multiple Snakefiles:

1. **Setup Workflow** (`setup.smk`): One-time installation of reference data and tools (~21GB)
2. **Analysis Workflows**:
   - `Snakefile`: Standard analysis workflow for local/single-node execution
   - `Snakefile.cluster`: Optimized workflow for cluster environments with job-specific resource allocation

---

## Installation

### Prerequisites

- **Linux operating system** (tested on Ubuntu/CentOS)
- **Snakemake** >= 7.3
- **Conda/Mamba** for dependency management
- **50GB free disk space** (21GB for references + analysis space)
- **Minimum 50GB RAM** for STAR alignment

### Step 1: Install Miniconda and Mamba

```bash
# Download and install Miniconda
cd ~/software
wget https://repo.anaconda.com/miniconda/Miniconda3-py311_23.11.0-2-Linux-x86_64.sh
sh Miniconda3-py311_23.11.0-2-Linux-x86_64.sh

# Install mamba for faster dependency resolution
conda install mamba=1.5.8 -n base -c conda-forge
```

### Step 2: Download IntegrateALL

```bash
# Clone the repository
git clone https://github.com/NadineWolgast/IntegrateALL.git
cd IntegrateALL

# Or download and extract the zip file
wget https://github.com/NadineWolgast/IntegrateALL/archive/refs/heads/main.zip
unzip main.zip && cd IntegrateALL-main
```

### Step 3: Create Environment

```bash
# Create and activate the conda environment
conda activate base
mamba env create --name integrateall --file environment.yaml
conda activate integrateall
```

---

## Configuration

Edit the `config.yaml` file to match your system:

```yaml
absolute_path: /absolute/path/to/IntegrateALL   # No trailing slash!
star_mem: 50000    # Minimum 50GB RAM for STAR
threads: 4         # Adjust to your CPU cores
```

**⚠️ Important**: Use absolute paths without trailing slashes to avoid errors.

---

## Setup (Run Once)

Install all required reference data and tools:

```bash
snakemake --snakefile setup.smk --cores 4
```

**This setup downloads (~21GB total):**
- Reference genome and annotations (~16GB)
- STAR genome index (~1GB)
- RNAseqCNV reference data (~50MB)  
- FusionCatcher database (~4.4GB)
- R packages (ALLCatchR, RNAseqCNV)
- Arriba draw_fusions tool (~10MB)

**⏱️ Time**: 1-3 hours (depending on internet speed)

Setup only runs once - subsequent executions skip existing components.

---

## Usage

### Sample Preparation

1. **Prepare sample sheet**: Edit `samples.csv` with your FASTQ file paths:

```csv
sample_id,left,right
Sample_01,/path/to/sample_01_R1.fastq.gz,/path/to/sample_01_R2.fastq.gz
Sample_02,/path/to/sample_02_R1.fastq.gz,/path/to/sample_02_R2.fastq.gz
```

2. **Test with provided data** (optional):

```csv
sample_id,left,right
Test,/path/to/IntegrateALL/data/samples/sub1_new.fq.gz,/path/to/IntegrateALL/data/samples/sub2_new.fq.gz
```

### Running the Analysis

#### Workflow Selection

Choose the appropriate Snakefile for your environment:

- **`Snakefile`**: Standard workflow for local execution or single-node processing
- **`Snakefile.cluster`**: Standard cluster workflow that always runs both Arriba and FusionCatcher for comprehensive fusion detection
- **`Snakefile.conditional`**: Optimized cluster workflow with Arriba-first approach that only runs FusionCatcher if initial classification fails, significantly reducing runtime for samples with clear driver fusions

#### Local Execution

1. **Dry run** (preview jobs):
```bash
snakemake -n --use-conda --conda-frontend conda
```

2. **Full analysis**:
```bash
snakemake --cores 20 --use-conda --conda-frontend conda
```

**Note**: The `--conda-frontend conda` flag prevents libmamba-related errors during conda environment creation.

3. **Single sample**:
```bash
snakemake --cores 4 --use-conda --conda-frontend conda STAR_output/YOUR_SAMPLE_ID/Aligned.sortedByCoord.out.bam
```

### Cluster Execution

**Option 1: Standard cluster workflow (comprehensive fusion detection)**
```bash
snakemake -s Snakefile.cluster --cores 20 --use-conda --conda-frontend conda
```

**Option 2: Optimized cluster workflow (Arriba-first approach)**
```bash
snakemake -s Snakefile.conditional --cores 20 --use-conda --conda-frontend conda
```

> **Recommendation**: Use `Snakefile.conditional` for faster processing when you expect clear driver fusions. Use `Snakefile.cluster` for comprehensive analysis when you need maximum sensitivity or are analyzing challenging samples.

**Option 3: SBATCH submission script**

First, configure the submission script for your cluster:
```bash
# Edit submit_snakemake_job.sh to match your cluster configuration
cp submit_snakemake_job.sh my_submit_script.sh
# Adjust partition name, memory, CPU count, and conda paths
```

Then submit the job:
```bash
sbatch my_submit_script.sh
```

**Option 4: Simple SLURM:**
```bash
srun -c 20 --mem 100G snakemake --cores 20 --use-conda --conda-frontend conda
```

**Option 5: SLURM Executor:**
```bash
snakemake --slurm --default-resources mem_mb=5000 threads=4 slurm_partition=YOUR_PARTITION --jobs 200 --use-conda --conda-frontend conda --keep-going
```

---

## Outputs

### Interactive Reports
- `interactive_output/{SAMPLE_ID}/output_report_{SAMPLE_ID}.html`

### Analysis Results
- **Alignments**: `STAR_output/{SAMPLE_ID}/`
- **Fusions**: `fusions/{SAMPLE_ID}.pdf`, `fusioncatcher_output/{SAMPLE_ID}/`
- **Fusion Intersects**: `data/fusion_intersect/{SAMPLE_ID}.csv`
- **Variants**: `Variants_RNA_Seq_Reads/{SAMPLE_ID}/`
- **CNV Analysis**: `RNAseqCNV_output/{SAMPLE_ID}/`
- **Classification**: `allcatch_output/{SAMPLE_ID}/predictions.tsv`
- **Karyotype Prediction**: `karyotype_prediction/{SAMPLE_ID}.csv`
- **Quality Control**: `qc/fastqc/{SAMPLE_ID}/`, `qc/multiqc/{SAMPLE_ID}/`
- **Individual Classification**: `Final_classification/{SAMPLE_ID}_output_report.csv`
- **Aggregated Summary**: `Final_classification/Aggregated_output_curation.csv`

### Expression Data
- **TPM**: `data/tpm/{SAMPLE_ID}.tsv`
- **CPM**: `data/cpm/{SAMPLE_ID}.tsv`
- **Counts**: `data/counts/{SAMPLE_ID}.tsv`
- **Raw counts**: `STAR_output/{SAMPLE_ID}/ReadsPerGene.out.tab`
- **Combined counts**: `data/combined_counts/ensemble_counts.tsv`, `data/combined_counts/gene_counts.tsv`

---

## Advanced Usage

### Partial Analysis

Run only specific analysis modules by editing the `rule all` section in `Snakefile`:

```python
rule all:
    input:
        "check_samples.txt",
        expand("fusioncatcher_output/{sample_id}/final-list_candidate-fusion-genes.txt", 
               sample_id=list(samples.keys())),
        # Uncomment desired outputs:
        # expand("STAR_output/{sample_id}/Aligned.sortedByCoord.out.bam", sample_id=list(samples.keys())),
        # expand("allcatch_output/{sample_id}/predictions.tsv", sample_id=samples.keys()),
        # expand("interactive_output/{sample}/output_report_{sample}.html", sample=list(samples.keys()))
```

### Environment Management

```bash
# Deactivate environment
conda deactivate

# Reactivate for analysis
conda activate integrateall
```

---

## Troubleshooting

### Common Issues

1. **"Snakefile not found"**: Ensure you're in the IntegrateALL directory
2. **Memory errors**: Increase `star_mem` in `config.yaml` (minimum 50000)
3. **Permission denied**: Check file paths and permissions in `samples.csv`
4. **Environment conflicts**: Remove and recreate the conda environment

### Support

For issues and bug reports, please use the GitHub issue tracker.

---

## Citation

If you use IntegrateALL in your research, please cite:
https://doi.org/10.1101/2025.09.25.673987
