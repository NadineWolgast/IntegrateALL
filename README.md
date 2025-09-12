# IntegrateALL

The IntegrateALL is a diagnostic RNA-Seq Fusion Pipeline for acute lymphoblastic leukemia samples. It provides the user with an interactiv report, containing subtype prediction, virtual karyotype, fusions, quality control, single nucleotide variants using a snakemake workflow management system. 

## Overview
![Blast-o-Matic-Fusionator](Pipeline_Overview.png?raw=true)

## Rule Flowchart
![Blast-o-Matic-Fusionator-Flowchart](pipeline_rule_flowchart.png?raw=true)
##  Prerequisites

The pipeline currently requires at least snakemake-minimal >= 7.3 and mamba, to be able to use mamba for dependency management.
Install miniconda locally and install mamba using the following commands:
```bash
cd ~/path/to/my/software
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
sh Miniconda3-latest-Linux-x86_64.sh
conda install mamba -n base -c conda-forge
```

Then, download the gihub repository of IntegrateALL with the command: 
```bash
git clone https://github.com/NadineKraft/IntegrateALL.git 
```

or download and unpack the zip from https://github.com/NadineKraft/IntegrateALL/archive/refs/heads/main.zip
into the directory from where you want it to run - this can be a different directory as the one where your data is stored.

Snakemake will be installed with all its dependencies in an isolated software environment via

```bash
cd /path/to/IntegrateALL
conda activate base
mamba env create --name integrateall --file environment.yaml
conda activate integrateall
```

You can deactivate the environment when you don't need it anymore with 

```bash
conda deactivate 
```
but keep it activated if you want to execute the next steps.

## Pipeline Setup (Two-Workflow Architecture)

IntegrateALL uses a two-workflow architecture for improved stability and efficiency:

1. **Setup Workflow** (`setup.smk`): One-time installation of reference data and tools
2. **Analysis Workflow** (`Snakefile`): Sample analysis using the installed references

### Step 1: Configuration
Change the path in config.yaml file to point to the **absolute path** where you've installed the pipeline:

```yaml
absolute_path: /absolute/path/to/IntegrateALL   # For example: /home/IntegrateALL
star_mem: 50000 # You can increase the amount of memory, but 50 GB is the minimum
threads: 4 # You can increase the amount of threads, but 4 is the minimum
```
**Don't** put an extra slash after the directory or it will throw an error.


```bash
# NOTE: This old installation method has been replaced by setup.smk
# Please see Step 2 below for the new two-workflow architecture
``` 
**This will need ~60 GB of space and takes ~6 hours**

### Step 2: One-Time Setup (Required Before First Use)

Install all required pipeline tools and reference data with the setup workflow:

```bash
snakemake --snakefile setup.smk --cores 4
``` 

**This setup will download and install (~21GB total):**
- Reference genome and annotations (~16GB)
- STAR genome index (~1GB) 
- RNAseqCNV reference data (~50MB)
- FusionCatcher database (~4.4GB)
- R packages (ALLCatchR, RNAseqCNV)
- Arriba draw_fusions tool (~10MB)

**The setup takes approximately 1-3 hours depending on your internet connection.**

All reference files are write-protected to prevent accidental deletion. The setup only needs to be run once - subsequent runs will skip already installed components.

## Step 3: Test and Run the Analysis Pipeline
Copy or move your FASTQ files into **ONE** directory and change the samples.csv file to point to your actual samples and adjust the sample_id names. 
You can also test the pipeline with the provided samples (sub1_new.fq.gz and sub2_new.fq.gz) in the directory data/samples:

| sample_id   |      left     |  right |
|----------|:-------------:|------:|
| Test |  /path/to/IntegrateALL/data/samples/sub1_new.fq.gz	 | /path/to/IntegrateALL/data/samples/sub2_new.fq.gz |


**After completing the setup in Step 2**, you can run the analysis workflow.

To test and see the analysis pipeline's execution jobs before running you can use:
```bash
  snakemake -n
```
This will list the resulting jobs and reasons for them. If you agree with them, run the analysis with:
```bash
  snakemake --use-conda --cores 20
```
You can adjust the amount of cores to your available amount with **--cores all**. This will allow parallelization and faster execution for multiple jobs. 
This command will invoke the whole analysis for all samples in your samples.csv.

If you want to run only a selection of the pipeline analysis methods you can change the command for example to:

```bash
  snakemake --use-conda --cores 1 allowed_rules run_fusioncatcher
```
But you will need to adjust the **rule all** inside the Snakemake file like this: 

```python
  rule all:
    input:
        "check_samples.txt",
        #expand("STAR_output/{sample_id}/Aligned.sortedByCoord.out.bam.bai", sample_id=list(samples.keys())),
        #expand("fusions/{sample_id}.pdf",sample_id=samples.keys()),
        #expand("STAR_output/{sample_id}/Aligned.sortedByCoord.out.bam",sample_id=list(samples.keys())),
        #expand("multiqc/{sample}/multiqc_data/multiqc_fastqc.txt", sample=fastq_dataframe['sample_id']),
        expand("fusioncatcher_output/{sample_id}/final-list_candidate-fusion-genes.txt",sample_id=list(samples.keys())),
        #expand("RNAseqCNV_output/{sample_id}",sample_id=samples.keys()),
        #expand("data/tpm/{sample_id}.tsv", sample_id=list(samples.keys())),
        #expand("data/cpm/{sample_id}.tsv", sample_id=list(samples.keys())),
        #expand("pysamstats_output_dir/{sample_id}/", sample_id=list(samples.keys())),
        #expand("comparison/{sample_id}.csv", sample_id= samples.keys()),
        #expand("data/single_counts/{sample_id}.txt", sample_id=samples.keys()),
        #expand("data/vcf_files/{sample_id}.tsv",sample_id=samples.keys()),
        #expand("allcatch_output/{sample_id}/predictions.tsv", sample_id= samples.keys()),
        #expand("aggregated_output/{sample}.csv", sample=list(samples.keys())),
        #expand("interactive_output/{sample}/output_report_{sample}.html",  sample=list(samples.keys()))
```

This will run only the analysis Fusioncatcher for all samples in your samples.csv file.

You can also run a single analysis for only one of your samples.
For example, if you want the STAR Alignment output for only one sample you can change **YOUR_SAMPLE_ID** to one of your 
actual sample_ids from the samples.csv file and run the following command:
```bash
  snakemake --use-conda --cores 4 STAR_output/YOUR_SAMPLE_ID/Aligned.sortedByCoord.out.bam
```
You don't need to adjust the Snakemake file for this.

If you want to run the pipeline on a cluster with slurm you can change the command to match your available resources and run it with:
```bash
  srun -c 20 --mem 100G snakemake --use-conda --cores 20
```

Or you can use the slurm executor plugin to let SLURM manage the jobs with:
```bash
  snakemake --slurm --default-resources mem_mb=5000 threads=4 slurm_partition=<YOUR PARTITION> --jobs 200 --use-conda --keep-going
```

The pipeline will output an interactive report for each of your samples in the folder `/path/to/the/pipeline/interactive_output/**YOUR_SAMPLE_ID**/output_report_YOUR_SAMPLE_ID.html` with the necessary result files. 
If you want to process the output further you will find all produced data in the individual output folders:
* Variants_RNA_Seq_Reads/YOUR_SAMPLE_ID (This contains the output generated by GATK)
* allcatch_output/YOUR_SAMPLE_ID/predictions.tsv
* fastqc/YOUR_SAMPLE_ID_left and fastqc/YOUR_SAMPLE_ID_right
* fusions/YOUR_SAMPLE_ID.tsv and YOUR_SAMPLE_ID.pdf (This contains the output generated by ARRIBA)
* fusioncatcher_output/YOUR_SAMPLE_ID
* multiqc/YOUR_SAMPLE_ID_left and multiqc/YOUR_SAMPLE_ID_right
* pysamstats_output_dir/YOUR_SAMPLE_ID
* RNAseqCNV_output/YOUR_SAMPLE_ID
* STAR_output/YOUR_SAMPLE_ID
* aggregated_output/YOUR_SAMPLE_ID.csv (contains a summary of the above results)
* Hotspots/YOUR_SAMPLE_ID

Furthermore, the pipline produces the following files for your downstream analysis:   
* data/tpm/YOUR_SAMPLE_ID.tsv
* data/cpm/YOUR_SAMPLE_ID.tsv
* data/vcf_files/YOUR_SAMPLE_ID.tsv
* data/counts/YOUR_SAMPLE_ID.tsv
* data/reads_per_gene/YOUR_SAMPLE_IDReadsPerGene.out.tab

## NEW: Two-Workflow Architecture

IntegrateALL now uses **two separate workflows** for better stability:

### 1. Setup Workflow (Run Once)
```bash
snakemake --snakefile setup.smk --use-conda --cores all
```
- Downloads and installs all references (~21GB)
- Only needed once or when updating references
- Much faster than the old `install_all` method

### 2. Analysis Workflow (For Each Dataset)  
```bash
snakemake --use-conda --cores all
```
- Runs the main RNA-seq analysis
- No more unexpected reference downloads during analysis
- More stable and faster pipeline starts

**Migration**: If you used the old system, just run the setup workflow once, then use the main workflow as before.

