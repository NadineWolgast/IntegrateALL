#!/usr/bin/env python3
"""
Utility functions for IntegrateALL Snakemake pipeline.
"""

import os
from scripts.create_sample_dataframe import create_sample_dataframe
from scripts.generate_files import generate_files
from scripts.validate_input import validate_input
from scripts.get_ctat_input_files import get_ctat_input_files


def creating_folders(specific_folder: str) -> str:
    """
    Create a folder if it does not exist and ensure it ends with '/'.
    
    :param specific_folder: The folder path to be created
    :return: The created folder path with trailing '/'
    """
    if specific_folder != '':
        if specific_folder[-1] != '/':
            specific_folder = specific_folder + '/'

        specific_folder = os.path.expanduser(specific_folder)
        if not os.path.exists(specific_folder):
            os.makedirs(specific_folder)
    return specific_folder


def setup_directories():
    """Setup required directories for the pipeline."""
    creating_folders('STAR_output')
    creating_folders('Variants_RNA_Seq_Reads')
    creating_folders('refs/fusioncatcher')
    creating_folders('refs/GATK')


def load_samples(sample_file: str) -> dict:
    """
    Load samples from CSV file into dictionary.
    
    :param sample_file: Path to samples.csv file
    :return: Dictionary with sample_id -> (R1, R2) pairs
    """
    samples = {}
    with open(sample_file, "r") as f:
        next(f)  # Skip header
        for line in f:
            sample_id, left, right = line.strip().split(",")
            samples[sample_id] = (left, right)
    return samples


def initialize_pipeline(config):
    """
    Initialize the IntegrateALL pipeline.
    
    :param config: Snakemake config dictionary
    :return: Tuple of (samples_dict, absolute_path)
    """
    # Setup directories
    setup_directories()
    
    # Load samples
    sample_file = config["sample_file"]
    samples = load_samples(sample_file)
    
    # Generate RNASeqCNV config files
    generate_files(sample_file, 'data/')
    
    # Get absolute path
    absolute_path = config["absolute_path"]
    
    return samples, absolute_path