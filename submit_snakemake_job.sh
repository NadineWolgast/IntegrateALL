#!/bin/bash
#SBATCH --job-name=IntegrateALL
#SBATCH --partition=normal         # Adjust partition name for your cluster
#SBATCH --time=48:00:00            # Adjust time limit as needed
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8          # Adjust CPU count as needed
#SBATCH --mem=32G                  # Adjust memory as needed
#SBATCH --output=snakemake_%j.out
#SBATCH --error=snakemake_%j.err

# Print job info
echo "Starting IntegrateALL pipeline on $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Working directory: $(pwd)"

# Load required modules (adjust as needed for your cluster)
# module load conda
# module load snakemake

# Activate conda environment (adjust path as needed for your system)
# Option 1: If using miniconda/anaconda
# source ~/miniconda3/etc/profile.d/conda.sh
# Option 2: If using module system
# module load conda
conda activate integrateall  # Change to your environment name

# Check if environment is activated
echo "Active conda environment: $CONDA_DEFAULT_ENV"

# Set Snakemake parameters
CORES=8
SNAKEFILE="Snakefile.cluster"
CONFIG="config.yaml"

# Run Snakemake with cluster-optimized settings
echo "Starting Snakemake pipeline..."
snakemake \
    -s $SNAKEFILE \
    --configfile $CONFIG \
    --cores $CORES \
    --latency-wait 10 \
    --rerun-incomplete \
    --keep-going \
    --printshellcmds \
    --reason \
    --stats snakemake_stats.json \
    --use-conda \
    --conda-frontend mamba \
    2>&1 | tee snakemake_$(date +%Y%m%d_%H%M%S).log

# Print completion info
echo "Pipeline finished on $(date)"
echo "Exit code: $?"