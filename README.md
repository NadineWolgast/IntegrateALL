# Blast-o-Matic-Fusioninator
Snakemake diagnostic RNA-Seq Fusion Pipline for ALL

Setup:
See: https://snakemake.readthedocs.io/en/stable/getting_started/installation.html for installing snakemake
conda activate base
mamba env create --name Blast-o-Matic-Fusioninator --file environment.yaml
conda activate Blast-o-Matic-Fusioninator

Edit samples.csv with the paths to your FASTQ files.
Edit config.yaml with the correct paths.


Run example:
snakemake --use-conda --cores 10 
