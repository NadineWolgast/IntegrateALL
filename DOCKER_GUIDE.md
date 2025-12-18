# Docker Guide for IntegrateALL

This guide explains how to build and run IntegrateALL using Docker containers.

## Table of Contents
- [Quick Start](#quick-start)
- [Building the Docker Image](#building-the-docker-image)
- [Running the Pipeline](#running-the-pipeline)
- [Data Management](#data-management)
- [Production Deployment](#production-deployment)
- [Troubleshooting](#troubleshooting)

---

## Quick Start

```bash
# 1. Build the Docker image
docker build -t integrateall:latest .

# 2. Run setup (one-time, ~21GB download)
docker run -v $(pwd)/refs:/IntegrateALL/refs \
           -v $(pwd)/config.yaml:/IntegrateALL/config.yaml \
           integrateall:latest \
           snakemake --snakefile setup.smk --cores 4 --use-conda

# 3. Run analysis
docker run -v $(pwd)/data:/IntegrateALL/data \
           -v $(pwd)/refs:/IntegrateALL/refs \
           -v $(pwd)/config.yaml:/IntegrateALL/config.yaml \
           -v $(pwd)/samples.csv:/IntegrateALL/samples.csv \
           -v $(pwd)/output:/IntegrateALL/output \
           integrateall:latest \
           snakemake --cores 4 --use-conda
```

---

## Building the Docker Image

### Option 1: Use the Existing Dockerfile

The repository includes a basic Dockerfile:

```bash
cd IntegrateALL
docker build -t integrateall:latest .
```

**Build time**: ~30-60 minutes depending on your connection and CPU.

### Option 2: Improved Multi-Stage Dockerfile (Recommended)

For production use, create an optimized Dockerfile:

```dockerfile
# Dockerfile.production
FROM mambaorg/micromamba:1.5.8 AS builder

# Copy environment file
COPY environment.yaml /tmp/environment.yaml

# Create conda environment
RUN micromamba create -n integrateall -f /tmp/environment.yaml -y && \
    micromamba clean --all --yes

# Final stage
FROM ubuntu:22.04

# Install minimal dependencies
RUN apt-get update && apt-get install -y \
    wget bzip2 curl git cmake \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Copy conda environment from builder
COPY --from=builder /opt/conda /opt/conda

# Set working directory
WORKDIR /IntegrateALL

# Copy pipeline files
COPY Snakefile Snakefile.cluster setup.smk config.yaml ./
COPY scripts/ ./scripts/
COPY envs/ ./envs/
COPY data/annotation/ ./data/annotation/

# Set environment variables
ENV PATH=/opt/conda/envs/integrateall/bin:$PATH
ENV CONDA_PREFIX=/opt/conda/envs/integrateall

# Default command
CMD ["/bin/bash"]
```

Build the production image:

```bash
docker build -f Dockerfile.production -t integrateall:production .
```

### Option 3: Build with Docker Compose

Create a `docker-compose.yml`:

```yaml
version: '3.8'

services:
  integrateall:
    build:
      context: .
      dockerfile: Dockerfile
    image: integrateall:latest
    container_name: integrateall_pipeline
    volumes:
      - ./data:/IntegrateALL/data
      - ./refs:/IntegrateALL/refs
      - ./output:/IntegrateALL/output
      - ./config.yaml:/IntegrateALL/config.yaml
      - ./samples.csv:/IntegrateALL/samples.csv
    working_dir: /IntegrateALL
    command: bash
    # Resource limits (adjust based on your system)
    deploy:
      resources:
        limits:
          cpus: '8'
          memory: 60G
        reservations:
          cpus: '4'
          memory: 50G
```

Build and run:

```bash
docker-compose build
docker-compose run integrateall
```

---

## Running the Pipeline

### Step 1: Prepare Configuration

**On your host machine**, prepare the configuration files:

```bash
# Create data directory structure
mkdir -p data/samples refs output

# Edit config.yaml
cat > config.yaml << EOF
absolute_path: /IntegrateALL
sample_file: samples.csv
threads: 4
star_mem: 50000
gatk_mem: 8000
# ... other settings
EOF

# Create samples.csv
cat > samples.csv << EOF
sample_id,R1,R2
sample1,/IntegrateALL/data/samples/sample1_R1.fastq.gz,/IntegrateALL/data/samples/sample1_R2.fastq.gz
sample2,/IntegrateALL/data/samples/sample2_R1.fastq.gz,/IntegrateALL/data/samples/sample2_R2.fastq.gz
EOF
```

**⚠️ Important**: Use container paths (`/IntegrateALL/...`) in `samples.csv`, not host paths!

### Step 2: Run Setup (One-Time)

Install reference data inside the container:

```bash
docker run --rm \
  -v $(pwd)/refs:/IntegrateALL/refs \
  -v $(pwd)/config.yaml:/IntegrateALL/config.yaml \
  integrateall:latest \
  bash -c "source /opt/conda/etc/profile.d/conda.sh && \
           conda activate IALL && \
           snakemake --snakefile setup.smk --cores 4 --use-conda --conda-frontend conda"
```

This downloads ~21GB of reference data to `./refs/`.

### Step 3: Run Analysis

Run the complete pipeline:

```bash
docker run --rm \
  -v $(pwd)/data:/IntegrateALL/data \
  -v $(pwd)/refs:/IntegrateALL/refs \
  -v $(pwd)/config.yaml:/IntegrateALL/config.yaml \
  -v $(pwd)/samples.csv:/IntegrateALL/samples.csv \
  -v $(pwd):/IntegrateALL \
  -w /IntegrateALL \
  integrateall:latest \
  bash -c "source /opt/conda/etc/profile.d/conda.sh && \
           conda activate IALL && \
           snakemake --cores 4 --use-conda"
```

### Step 4: Interactive Mode

For debugging or manual work:

```bash
docker run -it --rm \
  -v $(pwd):/IntegrateALL \
  -w /IntegrateALL \
  integrateall:latest \
  bash

# Inside the container:
conda activate IALL
snakemake --cores 4 --use-conda
```

---

## Data Management

### Volume Mounts

Mount these directories for persistent data:

| Host Path | Container Path | Purpose |
|-----------|----------------|---------|
| `./data/samples/` | `/IntegrateALL/data/samples/` | Input FASTQ files |
| `./refs/` | `/IntegrateALL/refs/` | Reference data (~21GB) |
| `./output/` | `/IntegrateALL/` | All pipeline outputs |
| `./config.yaml` | `/IntegrateALL/config.yaml` | Configuration file |
| `./samples.csv` | `/IntegrateALL/samples.csv` | Sample metadata |

### Directory Structure Inside Container

```
/IntegrateALL/
├── data/
│   ├── samples/           # Input FASTQ files
│   ├── vcf_files/         # Generated VCF files
│   └── tpm/               # Gene expression data
├── refs/                  # Reference data (mounted)
│   ├── GATK/
│   └── fusioncatcher/
├── STAR_output/           # Alignment results
├── fusions/               # Fusion detection results
├── Hotspots/              # Hotspot analysis
├── Final_classification/  # Classification results
├── interactive_output/    # HTML reports
└── logs/                  # Log files
```

### Handling Large Data

For large datasets, use bind mounts for better performance:

```bash
# Instead of copying, mount directly
docker run --rm \
  -v /path/to/large/fastq/dir:/IntegrateALL/data/samples:ro \
  -v $(pwd)/output:/IntegrateALL \
  integrateall:latest \
  snakemake --cores 8
```

---

## Production Deployment

### 1. Docker Image Registry

Push to a registry for deployment:

```bash
# Tag the image
docker tag integrateall:latest myregistry.com/integrateall:v1.0

# Push to registry
docker push myregistry.com/integrateall:v1.0

# Pull on production server
docker pull myregistry.com/integrateall:v1.0
```

### 2. Resource Management

Specify resource limits:

```bash
docker run --rm \
  --cpus="8" \
  --memory="60g" \
  --memory-swap="80g" \
  -v $(pwd):/IntegrateALL \
  integrateall:latest \
  snakemake --cores 8 --resources mem_mb=60000
```

### 3. Running Multiple Samples in Parallel

Use Snakemake's job scheduling:

```bash
docker run --rm \
  --cpus="16" \
  --memory="100g" \
  -v $(pwd):/IntegrateALL \
  integrateall:latest \
  bash -c "source /opt/conda/etc/profile.d/conda.sh && \
           conda activate IALL && \
           snakemake --cores 16 --jobs 4 --use-conda"
```

### 4. Batch Processing Script

Create a wrapper script `run_pipeline.sh`:

```bash
#!/bin/bash
set -e

SAMPLES_FILE="$1"
OUTPUT_DIR="$2"
THREADS="${3:-8}"

if [ -z "$SAMPLES_FILE" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Usage: $0 <samples.csv> <output_dir> [threads]"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Run pipeline
docker run --rm \
  --cpus="$THREADS" \
  --memory="60g" \
  -v $(pwd)/data:/IntegrateALL/data \
  -v $(pwd)/refs:/IntegrateALL/refs \
  -v "$SAMPLES_FILE":/IntegrateALL/samples.csv \
  -v "$OUTPUT_DIR":/IntegrateALL/output \
  -v $(pwd)/config.yaml:/IntegrateALL/config.yaml \
  integrateall:latest \
  bash -c "source /opt/conda/etc/profile.d/conda.sh && \
           conda activate IALL && \
           snakemake --cores $THREADS --use-conda 2>&1 | tee /IntegrateALL/output/pipeline.log"

echo "Pipeline completed! Results in: $OUTPUT_DIR"
```

Usage:

```bash
chmod +x run_pipeline.sh
./run_pipeline.sh samples.csv ./results 8
```

### 5. Cluster Integration

For HPC clusters with Docker/Singularity support:

```bash
# Convert Docker image to Singularity
singularity build integrateall.sif docker://integrateall:latest

# Run on cluster
singularity exec \
  --bind /data:/IntegrateALL/data \
  --bind /refs:/IntegrateALL/refs \
  integrateall.sif \
  snakemake --snakefile Snakefile.cluster --cores 16
```

---

## Troubleshooting

### Issue 1: Permission Errors

**Problem**: Files created by container have wrong ownership.

**Solution**: Run container with your user ID:

```bash
docker run --rm \
  --user $(id -u):$(id -g) \
  -v $(pwd):/IntegrateALL \
  integrateall:latest \
  snakemake --cores 4
```

### Issue 2: Out of Memory

**Problem**: STAR alignment fails with "out of memory".

**Solution**: Increase Docker memory limit:

```bash
# Check current limit
docker info | grep Memory

# Increase in Docker Desktop: Settings > Resources > Memory
# Or specify in docker run:
docker run --rm --memory="60g" --memory-swap="80g" ...
```

### Issue 3: Conda Environment Not Activated

**Problem**: Tools not found inside container.

**Solution**: Always activate the environment:

```bash
docker run --rm \
  integrateall:latest \
  bash -c "source /opt/conda/etc/profile.d/conda.sh && \
           conda activate IALL && \
           snakemake --version"
```

### Issue 4: Reference Data Missing

**Problem**: Pipeline fails with "reference not found".

**Solution**: Ensure refs directory is mounted and setup was run:

```bash
# Check if refs exist
docker run --rm \
  -v $(pwd)/refs:/IntegrateALL/refs \
  integrateall:latest \
  ls -lh /IntegrateALL/refs/

# Re-run setup if needed
docker run --rm \
  -v $(pwd)/refs:/IntegrateALL/refs \
  -v $(pwd)/config.yaml:/IntegrateALL/config.yaml \
  integrateall:latest \
  bash -c "source /opt/conda/etc/profile.d/conda.sh && \
           conda activate IALL && \
           snakemake --snakefile setup.smk --cores 4 --use-conda"
```

### Issue 5: Pipeline Hangs

**Problem**: Pipeline stops responding.

**Solution**: Check logs and use dry-run:

```bash
# Dry-run to check what will be executed
docker run --rm \
  -v $(pwd):/IntegrateALL \
  integrateall:latest \
  bash -c "source /opt/conda/etc/profile.d/conda.sh && \
           conda activate IALL && \
           snakemake -n --cores 4"

# Check specific rule logs
docker run --rm \
  -v $(pwd):/IntegrateALL \
  integrateall:latest \
  cat logs/star/sample1.log
```

---

## Advanced Configuration

### Custom Dockerfile with Pre-Downloaded References

To include references in the image (not recommended for development, but useful for distribution):

```dockerfile
FROM integrateall:latest

# Copy pre-downloaded references (if distributing)
COPY refs/ /IntegrateALL/refs/

# Make read-only to prevent accidental modification
RUN chmod -R a-w /IntegrateALL/refs/
```

### Environment Variables

Set environment variables for the pipeline:

```bash
docker run --rm \
  -e THREADS=8 \
  -e STAR_MEM=60000 \
  -v $(pwd):/IntegrateALL \
  integrateall:latest \
  bash -c "source /opt/conda/etc/profile.d/conda.sh && \
           conda activate IALL && \
           snakemake --cores \$THREADS --use-conda"
```

### Logging

Capture all output:

```bash
docker run --rm \
  -v $(pwd):/IntegrateALL \
  integrateall:latest \
  bash -c "source /opt/conda/etc/profile.d/conda.sh && \
           conda activate IALL && \
           snakemake --cores 4 --use-conda" \
  2>&1 | tee pipeline_$(date +%Y%m%d_%H%M%S).log
```

---

## Best Practices

1. **Always mount volumes**: Don't copy large data into containers
2. **Use specific tags**: Tag images with versions, not just `latest`
3. **Limit resources**: Specify CPU and memory limits to prevent system overload
4. **Check logs**: Always review logs in case of failures
5. **Clean up**: Remove stopped containers and unused images regularly
   ```bash
   docker container prune
   docker image prune
   ```
6. **Backup references**: The 21GB reference data should be backed up separately
7. **Version control**: Keep track of which Docker image was used for each analysis

---

## Support

For issues specific to Docker deployment, please check:
- [GitHub Issues](https://github.com/NadineWolgast/IntegrateALL/issues)
- [Main README](README.md) for pipeline-specific questions
- Docker documentation: https://docs.docker.com/

---

## Citation

If you use IntegrateALL in your research, please cite:
[Citation information to be added]
