# Blast-o-Matic-Fusioninator

Snakemake diagnostic RNA-Seq Fusion Pipline for ALL




## Setup
See: https://snakemake.readthedocs.io/en/stable/getting_started/installation.html for installing snakemake

# TODO: write about installing mamba
```bash
conda activate base
mamba env create --name Blast-o-Matic-Fusioninator --file environment.yaml
conda activate Blast-o-Matic-Fusioninator
```

### CTAT mutations

See: https://github.com/NCIP/ctat-mutations/wiki/CTAT-mutations-installation#ctat-mutations-genome-lib-installation for installing the CTAT mutations genome library

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

### STAR reference files
You can download your STAR reference and gtf file or use the pipeline to get them. You only need to adjust the paths inside the config.yaml to point where they are or should be stored.

### Adjust config.yaml
```bash
sample_file: samples.csv  # Change left and right paths to your FASTQ samples
                            and adjust the sample_id

rna_fusion_data_directory: /path/to/fusioncatcher/data/human_v102

star_ref: /path/to/STAR_indexfiles/GRCh38.primary_assembly.genome.fa
star_gtf: /path/to/STAR_indexfiles/gencode.v32.annotation.gtf
star_mem: 80000
threads: 20
genome_index: /path/to/genome_index/STAR/ensembl_94_100
star_tmp_directory: /path/to/STAR_tmp


ctat_input_directory: /path/to/data/samples
genome_lib: /path/to/ctat-mutationst/ctat_genome_lib_build_dir

# singularity
use-singularity: True
latency-wait: 60

```
## Run example

To deploy this project run

```bash
  snakemake --use-conda --cores 10 
```
