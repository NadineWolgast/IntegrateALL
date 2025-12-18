# Singularity/Apptainer Guide for IntegrateALL

This guide explains how to build and run IntegrateALL using Singularity containers on HPC clusters.

**Note:** Singularity is now called Apptainer in newer versions, but the commands are the same.

## Table of Contents
- [Why Singularity?](#why-singularity)
- [Building the Container](#building-the-container)
- [Running on HPC Clusters](#running-on-hpc-clusters)
- [SLURM Integration](#slurm-integration)
- [Troubleshooting](#troubleshooting)

---

## Why Singularity?

- ✅ **Designed for HPC**: No root privileges needed to run containers
- ✅ **Security**: Runs as your user (no privilege escalation)
- ✅ **Performance**: Direct access to host filesystems
- ✅ **Compatibility**: Can convert Docker images
- ✅ **Cluster-friendly**: Works with SLURM, PBS, SGE

---

## Building the Container

### Option 1: Build from Definition File (Recommended)

```bash
# On your local machine or build node (requires sudo/root)
cd IntegrateALL

# Build the container
sudo singularity build integrateall.sif integrateall.def

# This creates integrateall.sif (~2-3 GB)
```

### Option 2: Convert from Docker Image

If you already have the Docker image:

```bash
# Pull from Docker Hub (if you've pushed it there)
singularity pull docker://yourusername/integrateall:latest

# Or convert local Docker image
sudo singularity build integrateall.sif docker-daemon://integrateall:latest
```

### Option 3: Build on Cluster with Fakeroot

Some clusters allow fakeroot builds without sudo:

```bash
# Check if fakeroot is available
singularity --version

# Build with fakeroot
singularity build --fakeroot integrateall.sif integrateall.def
```

---

## Running on HPC Clusters

### Step 1: Transfer Container to Cluster

```bash
# Copy the .sif file to your cluster
scp integrateall.sif username@cluster.domain:/path/to/project/

# Or build on cluster if fakeroot is available
```

### Step 2: Setup (One-Time)

Create the reference data directory on the cluster:

```bash
# On cluster
cd /path/to/project
mkdir -p refs

# Run setup
singularity exec \
  --bind $(pwd)/refs:/IntegrateALL/refs \
  --bind $(pwd)/config.yaml:/IntegrateALL/config.yaml \
  --bind $(pwd):/IntegrateALL \
  --pwd /IntegrateALL \
  integrateall.sif \
  bash -c "source /opt/conda/etc/profile.d/conda.sh && \
           conda activate integrateall && \
           snakemake --snakefile setup.smk --cores 4 --use-conda"
```

**Important:** Setup downloads ~21GB of data. Consider:
- Using a shared filesystem for refs (don't duplicate per user)
- Running setup on a data transfer node
- Using high-bandwidth node for downloads

### Step 3: Prepare Your Data

```bash
# Create directory structure
mkdir -p data/samples

# Copy FASTQ files
cp /path/to/fastq/*.fastq.gz data/samples/

# Update samples.csv with cluster paths
cat > samples.csv << EOF
sample_id,R1,R2
sample1,/IntegrateALL/data/samples/sample1_R1.fastq.gz,/IntegrateALL/data/samples/sample1_R2.fastq.gz
EOF
```

### Step 4: Run Analysis

```bash
# Interactive test
singularity exec \
  --bind $(pwd):/IntegrateALL \
  --pwd /IntegrateALL \
  integrateall.sif \
  bash -c "source /opt/conda/etc/profile.d/conda.sh && \
           conda activate integrateall && \
           snakemake -n --cores 8"  # Dry run

# Full analysis
singularity exec \
  --bind $(pwd):/IntegrateALL \
  --pwd /IntegrateALL \
  integrateall.sif \
  bash -c "source /opt/conda/etc/profile.d/conda.sh && \
           conda activate integrateall && \
           snakemake --cores 8 --use-conda"
```

---

## SLURM Integration

### Option 1: Single SLURM Job

Submit the entire pipeline as one job:

```bash
#!/bin/bash
#SBATCH --job-name=integrateall
#SBATCH --output=integrateall_%j.log
#SBATCH --error=integrateall_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --partition=long

# Load singularity module (if needed)
module load singularity

# Set working directory
cd /path/to/project

# Run pipeline
singularity exec \
  --bind $(pwd):/IntegrateALL \
  --pwd /IntegrateALL \
  integrateall.sif \
  bash -c "source /opt/conda/etc/profile.d/conda.sh && \
           conda activate integrateall && \
           snakemake --cores $SLURM_CPUS_PER_TASK --use-conda"
```

Submit:
```bash
sbatch run_pipeline.slurm
```

### Option 2: Snakemake with SLURM Executor

Use Snakemake's SLURM executor to submit individual rules as jobs:

```bash
#!/bin/bash
#SBATCH --job-name=integrateall_controller
#SBATCH --output=controller_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=72:00:00

module load singularity

cd /path/to/project

# Snakemake submits jobs via SLURM
singularity exec \
  --bind $(pwd):/IntegrateALL \
  --pwd /IntegrateALL \
  integrateall.sif \
  bash -c "source /opt/conda/etc/profile.d/conda.sh && \
           conda activate integrateall && \
           snakemake --snakefile Snakefile.cluster \
                     --executor slurm \
                     --jobs 10 \
                     --default-resources slurm_account=your_account \
                     --use-conda"
```

### Option 3: Wrapper Script for Cluster

Create `run_singularity_cluster.sh`:

```bash
#!/bin/bash
# IntegrateALL Singularity Wrapper for HPC Clusters

set -e

# Configuration
CONTAINER="integrateall.sif"
CORES=16
MEMORY="64G"
TIME="48:00:00"
PARTITION="long"

# Check if container exists
if [ ! -f "$CONTAINER" ]; then
    echo "Error: Container $CONTAINER not found"
    exit 1
fi

# Function to submit SLURM job
submit_job() {
    local job_name=$1
    local command=$2

    cat << EOF | sbatch
#!/bin/bash
#SBATCH --job-name=$job_name
#SBATCH --output=${job_name}_%j.log
#SBATCH --cpus-per-task=$CORES
#SBATCH --mem=$MEMORY
#SBATCH --time=$TIME
#SBATCH --partition=$PARTITION

cd \$(pwd)

singularity exec \\
  --bind \$(pwd):/IntegrateALL \\
  --pwd /IntegrateALL \\
  $CONTAINER \\
  bash -c "source /opt/conda/etc/profile.d/conda.sh && \\
           conda activate integrateall && \\
           $command"
EOF
}

# Parse command
case "$1" in
    setup)
        echo "Submitting setup job..."
        submit_job "integrateall_setup" "snakemake --snakefile setup.smk --cores $CORES --use-conda"
        ;;
    run)
        echo "Submitting analysis job..."
        submit_job "integrateall_run" "snakemake --cores $CORES --use-conda"
        ;;
    *)
        echo "Usage: $0 {setup|run}"
        exit 1
        ;;
esac
```

Usage:
```bash
chmod +x run_singularity_cluster.sh
./run_singularity_cluster.sh setup
./run_singularity_cluster.sh run
```

---

## Advanced Configuration

### Binding Multiple Directories

```bash
singularity exec \
  --bind /scratch/user/data:/IntegrateALL/data \
  --bind /work/refs:/IntegrateALL/refs \
  --bind /home/user/project:/IntegrateALL/output \
  integrateall.sif \
  snakemake --cores 8
```

### Using Overlay Filesystems

For writable conda environments on read-only filesystems:

```bash
# Create overlay image
singularity overlay create --size 4096 overlay.img

# Use overlay
singularity exec \
  --overlay overlay.img \
  --bind $(pwd):/IntegrateALL \
  integrateall.sif \
  snakemake --cores 8
```

### GPU Support (if needed)

```bash
singularity exec \
  --nv \
  --bind $(pwd):/IntegrateALL \
  integrateall.sif \
  snakemake --cores 8
```

---

## Shared Reference Setup

For multi-user clusters, set up shared references:

```bash
# As admin or in shared directory
export SHARED_REFS=/shared/data/integrateall/refs

mkdir -p $SHARED_REFS

# Run setup once
singularity exec \
  --bind $SHARED_REFS:/IntegrateALL/refs \
  --bind config.yaml:/IntegrateALL/config.yaml \
  integrateall.sif \
  bash -c "source /opt/conda/etc/profile.d/conda.sh && \
           conda activate integrateall && \
           snakemake --snakefile setup.smk --cores 8 --use-conda"

# Make read-only
chmod -R a-w $SHARED_REFS

# Users bind the shared refs
singularity exec \
  --bind $SHARED_REFS:/IntegrateALL/refs:ro \
  --bind $(pwd):/IntegrateALL \
  integrateall.sif \
  snakemake --cores 8
```

---

## Troubleshooting

### Issue 1: "No space left on device"

**Problem:** Singularity cache fills up `/tmp`.

**Solution:**
```bash
# Set cache to larger partition
export SINGULARITY_CACHEDIR=/scratch/$USER/singularity_cache
mkdir -p $SINGULARITY_CACHEDIR
```

### Issue 2: "Permission denied" for conda

**Problem:** Conda tries to write to read-only container.

**Solution:** Use overlay or bind a writable directory:
```bash
singularity exec \
  --bind $(pwd)/.conda:/opt/conda \
  integrateall.sif \
  snakemake --cores 8
```

### Issue 3: Container is too large

**Problem:** 2-3GB .sif file is hard to transfer.

**Solution:**
```bash
# Build on cluster instead of transferring
# Or use singularity pull from repository:
singularity pull library://user/integrateall:latest
```

### Issue 4: "Module not found" errors

**Problem:** Conda environment not activated.

**Solution:** Always wrap commands:
```bash
bash -c "source /opt/conda/etc/profile.d/conda.sh && \
         conda activate integrateall && \
         your_command"
```

### Issue 5: Snakemake can't create conda envs

**Problem:** No write access to conda directory.

**Solution:** Set conda prefix to writable location:
```bash
export CONDA_PKGS_DIRS=$HOME/.conda/pkgs
mkdir -p $CONDA_PKGS_DIRS

singularity exec \
  --bind $CONDA_PKGS_DIRS:/opt/conda/pkgs \
  integrateall.sif \
  snakemake --use-conda --conda-prefix $(pwd)/.snakemake/conda
```

---

## Best Practices for Clusters

1. **Use scratch space** for intermediate files:
   ```bash
   export TMPDIR=/scratch/$USER/tmp
   mkdir -p $TMPDIR
   ```

2. **Set resource limits** in Snakefile.cluster

3. **Monitor jobs**:
   ```bash
   squeue -u $USER
   tail -f logs/star/*.log
   ```

4. **Clean up** after completion:
   ```bash
   # Keep only final results
   snakemake --delete-temp-output
   ```

5. **Test with dry-run** first:
   ```bash
   snakemake -n --cores 8
   ```

---

## Converting Between Docker and Singularity

### Docker → Singularity

```bash
# Method 1: Direct conversion
singularity build integrateall.sif docker://integrateall:latest

# Method 2: From local Docker
docker save integrateall:latest -o integrateall.tar
singularity build integrateall.sif docker-archive://integrateall.tar
```

### Singularity → Docker

```bash
# Not directly possible, but can rebuild from same Dockerfile
```

---

## Performance Tips

1. **Use local scratch for I/O intensive tasks**
2. **Bind specific directories** (not entire filesystem)
3. **Disable Apptainer sandboxing** if allowed: `--no-home`
4. **Use SSD/NVMe** for working directories
5. **Enable parallel execution** with `--jobs` in Snakemake

---

## Support

For Singularity-specific issues:
- Singularity documentation: https://sylabs.io/docs/
- Apptainer documentation: https://apptainer.org/docs/

For IntegrateALL issues:
- GitHub: https://github.com/NadineWolgast/IntegrateALL/issues

---

## Example: Complete Cluster Workflow

```bash
# 1. Build container locally
sudo singularity build integrateall.sif integrateall.def

# 2. Transfer to cluster
scp integrateall.sif user@cluster:/work/user/

# 3. On cluster: Setup
cd /work/user
mkdir -p refs data/samples

singularity exec \
  --bind $(pwd):/IntegrateALL \
  integrateall.sif \
  bash -c "source /opt/conda/etc/profile.d/conda.sh && \
           conda activate integrateall && \
           snakemake --snakefile setup.smk --cores 8 --use-conda"

# 4. Copy your FASTQ files
cp /path/to/data/*.fastq.gz data/samples/

# 5. Run analysis via SLURM
sbatch run_pipeline.slurm

# 6. Monitor
squeue -u $USER
tail -f slurm-*.out

# 7. Collect results
ls -lh Final_classification/
ls -lh interactive_output/
```
