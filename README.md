# Blast-o-Matic-Fusioninator

Snakemake diagnostic RNA-Seq Fusion Pipline for ALL

![Blast-o-Matic-Fusionator](Blast-o-Matic-Fusioninator.png?raw=true)


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
mamba env create --name Blast-o-Matic-Fusioninator --file environment.yaml
conda activate Blast-o-Matic-Fusioninator
```

You can deactivate the environment when you don't need it anymore with 

```bash
conda deactivate 
```
but keep it activated if you want to execute the next steps.

## Before you can run the pipeline:
You need to adjust the config.yaml file and install the missing genome librarys

### CTAT mutations

To install CTAT mutations genome library edit in the config.yaml the path ctat_genome_lib_build_dir: /path/where/the/ctat_genome_library_shall_be_installed  and run the following commands:
```bash
snakemake --cores 10 pull_ctat_mutations_singularity_image
snakemake --cores 10 install_ctat_mutations
snakemake --cores 10 run_ctat_genome_lib_builder
```
This will install the Plug-n-Play genome library needed for running CTAT mutations. Keep in mind, that this will need at least 78 GB space. 
if you already have a CTAT mutations Genome library installed, you can adjust the path ctat_genome_lib_build_dir in the config.yaml with the actual path and skip this installation. 

### ALLCatchR
Install the ALLCatchR with the command:
```bash
snakemake --cores 10 install_allcatchr
```

## RNASeqCNV 
Install RNASeqCnv with the command:
```bash
snakemake --use-conda --cores 10 install_rnaseq_cnv
```


### Fusioncatcher
See: https://github.com/ndaniel/fusioncatcher for installing and downloading fusioncatcher db:
```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda create -n fusioncatcher fusioncatcher
source activate fusioncatcher
download-human-db.sh
```
Now adjust in config.yaml the rna_fusion_data_directory with the installed path to the downloaded human_v102 directory.
```bash
rna_fusion_data_directory: /path/to/fusioncatcher/data/human_v102
```

### STAR reference files
You can download your STAR reference and gtf file or use the pipeline to get them. You only need to adjust the paths inside the config.yaml to point where they are or should be stored.

If you want the pipeline to download them, you need to edit the path in config.yaml for star_files to the path where the reference files should be downloadet and stored.
Then run 

```bash
  snakemake --cores 10 download_star_ref
```

Adjust now the paths of star_ref and star_gtf in the config.yaml to point to the actual files.
```bash
star_ref: /path/to/STAR_indexfiles/GRCh38.primary_assembly.genome.fa
star_gtf: /path/to/STAR_indexfiles/gencode.v32.annotation.gtf
```

Generate the STAR genome_index with

```bash
  snakemake --cores 10 index
```
Now you have all needed reference files to run the pipeline. 


## Run examples
Change the samples.csv file to point to your actual samples and change the sample_id names.

To test and see the pipelines execution jobs before running the pipeline you can run the command:
```bash
  snakemake -n
```
To see the resulting jobs and reasons for them. If you agree with them, you can run
them with:
```bash
  snakemake --use-singularity --use-conda --cores 10 
```
You can adjust the amount of cores to your available amount.
This command will invoke the whole analysis for all samples in your samples.csv.
If you want to run only a selection of the pipeline analysis methods you can change the command for example to:

```bash
  snakemake --use-singularity --use-conda --cores 10 allowed_rules run_allcatchr_on_single_count_files
```
This will run only the analysis ALLCatchR for all samples in your samples.csv file.

You can also run an analysis for only one sample.
For example, if you want the CTAT mutations output for only one sample you can change YOUR_SAMPLE_ID to one of your 
actual sample_ids from the samples.csv file and run the following command:
```bash
  snakemake --use-singularity --use-conda --cores 10 ctat_output_directory/YOUR_SAMPLE_ID
```

If you want to run the pipeline on a cluster with slurm for example you can run it with:
```bash
  srun -c 20 --mem 100G snakemake --use-singularity --use-conda --cores 20 --resources threads=200 -j 20
```
