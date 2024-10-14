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

# Download and install Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda && \
    rm Miniconda3-latest-Linux-x86_64.sh

# Set the PATH environment variable to include conda
ENV PATH=/opt/conda/bin:$PATH

# Install mamba in the base environment
RUN conda install mamba -n base -c conda-forge

# Install Snakemake using mamba
RUN mamba install -c bioconda -c conda-forge snakemake

# Create and activate the environment using mamba
RUN mamba env create --name IALL --file environment.yaml

# Initialize conda for bash
RUN conda init bash


# Set up the entrypoint to activate the conda environment
ENTRYPOINT ["/bin/bash", "-c", "source /opt/conda/etc/profile.d/conda.sh && conda activate IALL && exec bash"]

# Ensure conda environments are activated properly
CMD ["bash"]

