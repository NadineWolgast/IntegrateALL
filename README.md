# Blast-o-Matic-Fusioninator

Snakemake diagnostic RNA-Seq Fusion Pipline for ALL

![Blast-o-Matic-Fusionator](Pipeline.png?raw=true)


Includes the following rules:

![Rulegraph](rulegraph.svg?raw=true)

##  Prerequisites

For the pipeline to run, we need

    snakemake
    conda / miniconda / anaconda / mamba (highly recommended to use mamba, see below)

to be available.

The pipeline currently requires at least snakemake-minimal >= 7.3 and mamba, to be able to use mamba for dependency management.
In cluster environments, versions are often outdated, or tools not even available. Hence the recommendation to simply install your own miniconda locally for your user. Use wget to obtain the binary (probably you want the generic linux one), and execute it to install locally, log out, log in again to activate it, and also get mamba with:
```bash
cd ~/path/to/my/software
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
sh Miniconda3-latest-Linux-x86_64.sh
conda install mamba -n base -c conda-forge
```

Then, download the gihub repository of Blast-o-Matic-Fusioninator with the command: 
```bash
git clone https://github.com/NadineKraft/Blast-o-Matic-Fusioninator.git 
```
or download and unpack the zip from https://github.com/NadineKraft/Blast-o-Matic-Fusioninator/archive/refs/heads/main.zip
into the directory from where you want it to run - this can be a different directory as the one where your data is stored.
Snakemake will be installed with all dependencies needed to run it in an isolated software environment via

```bash
cd /path/to/Blast-o-Matic-Fusioninator
conda activate base
mamba env create --name bomf --file environment.yaml
conda activate bomf
```

You can deactivate the environment when you don't need it anymore with 

```bash
conda deactivate 
```
but keep it activated if you want to execute the next steps.

## Before you can run the pipeline:
You need to adjust the config.yaml file and install the missing genome librarys



## TL:DR
Change the paths in config.yaml file to absolute paths where you want the references installed:
```yaml
rna_fusion_data_directory: /refs/fusioncatcher/data/human_v102
star_files: /refs/STAR/
star_ref: /refs/STAR/GRCh38.primary_assembly.genome.fa
star_gtf: /refs/STAR/gencode.v32.annotation.gtf
genome_index: /refs/STAR/ensembl_94_100
star_tmp_directory: /refs/STAR/STAR_tmp
genome_lib:  /refs/ctat/ctat_genome_lib_build_dir/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir
ctat_genome_lib_build_dir: /refs/ctat/
```
And install all required pipeline tools with:

```bash
snakemake --use-conda --use-singularity --cores all install_all
```

## Install Fusioncatcher
You need to install and download the fusioncatcher db with:
```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda create -n fusioncatcher fusioncatcher
source activate fusioncatcher
download-human-db.sh
```
Now adjust in config.yaml the rna_fusion_data_directory with the installed path to the downloaded human_v102 directory.
```yaml
rna_fusion_data_directory: /path/to/fusioncatcher/data/human_v102
```
## Install only selected tools
If you don't want to install all tools and references for the pipeline because you already have some of them you can select the missing ones and install them individually:
### CTAT mutations

To install CTAT mutations genome library edit in the config.yaml the path

```yaml
ctat_genome_lib_build_dir: /path/where/the/ctat_genome_library_shall_be_installed # You need the absolute path here
```

and run the following commands:

```bash
snakemake --cores 1 pull_ctat_mutations_singularity_image
snakemake --cores 1 install_ctat_mutations
snakemake --cores 1 run_ctat_genome_lib_builder
```
This will install the Plug-n-Play genome library needed for running CTAT mutations. Keep in mind, that this will need at least 78 GB space. 
If you already have a CTAT mutations Genome library installed, you can adjust the path ctat_genome_lib_build_dir in the config.yaml with the actual path and skip this installation. 

### ALLCatchR
Install the ALLCatchR with the command:
```bash
snakemake --cores 1 install_allcatchr
```

## RNASeqCNV 
Install RNASeqCnv with the command:
```bash
snakemake --use-conda --cores 1 install_rnaseq_cnv
```


### Fusioncatcher
See: https://github.com/ndaniel/fusioncatcher for more information or install and download the fusioncatcher db with:
```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda create -n fusioncatcher fusioncatcher
source activate fusioncatcher
download-human-db.sh
```
Now adjust in config.yaml the rna_fusion_data_directory with the installed path to the downloaded human_v102 directory.
```yaml
rna_fusion_data_directory: /path/to/fusioncatcher/data/human_v102
```

### ARRIBA draw fusions
In order to produce arribas publication-quality visualizations of the transcripts involved in predicted fusions it needs to be installed with the command
```bash
snakemake --cores 1 --use-conda install_arriba_draw_fusions
```
This will download and install arrbia version 2.4.0 and its' database in the same directory as the pipeline.  

### STAR reference files
You can download your STAR reference and gtf file or use the pipeline to get them. You only need to adjust the paths inside the config.yaml to point where they are or should be stored.

If you want the pipeline to download them, you need to edit the path in config.yaml for star_files to the path where the reference files should be downloadet and stored.
Then run 

```bash
  snakemake --cores 1 download_star_ref
```

Adjust now the paths of star_ref and star_gtf in the config.yaml to point to the actual files.
```yaml
star_ref: /path/to/STAR_indexfiles/GRCh38.primary_assembly.genome.fa
star_gtf: /path/to/STAR_indexfiles/gencode.v32.annotation.gtf
```

Generate the STAR genome_index with

```bash
  snakemake --cores all index
```
Now you have all needed reference files and tools to run the pipeline. 


## Run examples
Copy or move your FASTQ files into one directory and change the samples.csv file to point to your actual samples and change the sample_id names. For example:

| sample_id   |      left     |  right |
|----------|:-------------:|------:|
| sample_XYZ |  absolute/path/to/your/sample_XYZ_R1_fastq.gz | absolute/path/to/your/sample_XYZ_R2_fastq.gz |
| another_sample |  absolute/path/to/your/another_sample_R1_fastq.gz | absolute/path/to/your/another_sample_R2_fastq.gz |
| example_sample |  absolute/path/to/your/example_sample_R1_fastq.gz | absolute/path/to/your/example_sample_R2_fastq.gz |

Inside the config.yaml file change this line to point to your actual FASTQ samples directory:
```yaml
ctat_input_directory: /absolute/path/to/your/samples/directory # You need the absolute path here!
```
**Don't** put an extra slash after the directory or CTAT will throw an error.


To test and see the pipelines execution jobs before running the pipeline you can run the command:
```bash
  snakemake -n
```
This will list the resulting jobs and reasons for them. If you agree with them, you can run them with:
```bash
  snakemake --use-singularity --use-conda --cores 1
```
You can adjust the amount of cores to your available amount with **--cores all**. This will allow parallelization and faster execution for multiple jobs. 
This command will invoke the whole analysis for all samples in your samples.csv.

If you want to run only a selection of the pipeline analysis methods you can change the command for example to:

```bash
  snakemake --use-conda --cores 1 allowed_rules run_fusioncatcher
```
But you will need to adjust the **rule all** in the Snakemake file like this: 

```python
  rule all:
    input:
        "check_samples.txt",
        #expand("STAR_output/{sample_id}/Aligned.sortedByCoord.out.bam.bai", sample_id=list(samples.keys())),
        #expand("fusions/{sample_id}.pdf",sample_id=samples.keys()),
        #expand("STAR_output/{sample_id}/Aligned.sortedByCoord.out.bam",sample_id=list(samples.keys())),
        #expand("multiqc/{sample}/multiqc_data/multiqc_fastqc.txt", sample=fastq_dataframe['sample_id']),
        expand("fusioncatcher_output/{sample_id}/final-list_candidate-fusion-genes.txt",sample_id=list(samples.keys())),
        #expand("ctat_output_directory/{sample_id}/{sample_id}.filtered.vcf.gz",sample_id=samples_test.keys()),
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
For example, if you want the CTAT mutations output for only one sample you can change **YOUR_SAMPLE_ID** to one of your 
actual sample_ids from the samples.csv file and run the following command:
```bash
  snakemake --use-singularity --use-conda --cores 1 ctat_output_directory/YOUR_SAMPLE_ID/
```
You don't need to adjust the Snakemake file for this.

If you want to run the pipeline on a cluster with slurm you can change the command to match your available resources and run it with:
```bash
  srun -c 10 --mem 100G snakemake --use-singularity --use-conda --cores 10 --resources threads=100 -j 10
```

The pipeline will output an interactive report for each of your samples in the folder `/path/to/the/pipeline/interactive_output/**YOUR_SAMPLE_ID**/output_report_YOUR_SAMPLE_ID.html` with the necessary result files. 
If you want to process the output further you will find all produced data in the individual output folders:
* ctat_output_directory/YOUR_SAMPLE_ID
* allcatch_output/YOUR_SAMPLE_ID
* fastqc/YOUR_SAMPLE_ID_left and fastqc/YOUR_SAMPLE_ID_right
* fusions/YOUR_SAMPLE_ID.tsv and YOUR_SAMPLE_ID.pdf (This contains the output generated by ARRIBA)
* fusioncatcher_output/YOUR_SAMPLE_ID
* multiqc/YOUR_SAMPLE_ID_left and multiqc/YOUR_SAMPLE_ID_right
* pysamstats_output_dir/YOUR_SAMPLE_ID
* RNAseqCNV_output/YOUR_SAMPLE_ID
* STAR_output/YOUR_SAMPLE_ID
*  aggregated_output/YOUR_SAMPLE_ID.csv (contains a summary of the above results)

Furthermore, the pipline produces the following files for your downstream analysis:   
* data/tpm/YOUR_SAMPLE_ID.tsv
* data/cpm/YOUR_SAMPLE_ID.tsv
* data/vcf_files/YOUR_SAMPLE_ID.tsv
* data/counts/YOUR_SAMPLE_ID.tsv

