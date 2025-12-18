# Use the official Ubuntu 20.04 base image
FROM ubuntu:20.04

# Set the working directory
WORKDIR /IntegrateALL

# Copy the Snakemake pipeline and environment file into the container
COPY . /IntegrateALL

# Set environment variable to prevent interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# Install necessary packages
RUN apt-get update && apt-get install -y \
    wget \
    bzip2 \
    curl \
    unzip \
    git \
    cmake \
    python3 \
    python3-pip \
    && apt-get clean

# Download and install Miniforge (conda-forge only, no ToS required)
RUN wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh && \
    bash Miniforge3-Linux-x86_64.sh -b -p /opt/conda && \
    rm Miniforge3-Linux-x86_64.sh

# Set the PATH environment variable to include conda
ENV PATH=/opt/conda/bin:$PATH

# Mamba is already included in Miniforge, configure conda to use only conda-forge
RUN conda config --add channels conda-forge && \
    conda config --add channels bioconda && \
    conda config --set channel_priority strict

# Install Snakemake using mamba
RUN mamba install -c bioconda -c conda-forge snakemake

# Create and activate the environment using mamba
RUN mamba env create --name IALL --file environment.yaml

# Initialize conda for bash
RUN conda init bash

# Create activation script for easy environment setup
RUN echo '#!/bin/bash\n\
source /opt/conda/etc/profile.d/conda.sh\n\
conda activate IALL\n\
exec "$@"' > /usr/local/bin/activate-env.sh && \
    chmod +x /usr/local/bin/activate-env.sh

# Set the entrypoint to activate conda environment
ENTRYPOINT ["/usr/local/bin/activate-env.sh"]

# Default command is bash
CMD ["bash"]

