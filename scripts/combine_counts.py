#!/usr/bin/env python3
"""
Optimized script to combine count data across all samples.
Creates ensemble counts and gene symbol counts matrices with dated outputs.
"""

import pandas as pd
import datetime
import os
import sys
import re
from pathlib import Path

def combine_count_files(count_files, gtf_file, ensemble_output, gene_output):
    """
    Combine count files into ensemble and gene symbol matrices.
    
    Args:
        count_files: List of count TSV files
        gtf_file: GTF file for gene ID to symbol mapping
        ensemble_output: Output path for ensemble counts matrix
        gene_output: Output path for gene symbol counts matrix
    """
    
    print(f"📊 Combining count data from {len(count_files)} samples")
    
    # Generate current date string
    current_date = datetime.datetime.now().strftime("%Y%m%d")
    
    try:
        # Create output directory
        os.makedirs(os.path.dirname(ensemble_output), exist_ok=True)
        
        # Create dated output filenames
        ensemble_dated = f"data/combined_counts/ensemble_counts_{current_date}.tsv"
        gene_dated = f"data/combined_counts/gene_counts_{current_date}.tsv"
        
        # Load gene symbol mapping from GTF file
        print("Loading gene symbol mapping from GTF file...")
        gene_mapping = load_gene_mapping(gtf_file)
        print(f"Loaded {len(gene_mapping):,} gene mappings")
        
        # Load and combine count data
        print("Loading and combining count data...")
        combined_data, sample_names = load_count_data(count_files)
        print(f"Combined data for {len(combined_data):,} genes across {len(sample_names)} samples")
        
        # Create ensemble counts matrix
        print("Creating ensemble counts matrix...")
        ensemble_df = create_count_matrix(combined_data, sample_names)
        ensemble_df.to_csv(ensemble_dated, sep='\t')
        
        # Create gene symbol counts matrix
        print("Creating gene symbol counts matrix...")
        gene_df = create_gene_symbol_matrix(combined_data, sample_names, gene_mapping)
        gene_df.to_csv(gene_dated, sep='\t')
        
        # Create symbolic links for Snakemake outputs
        create_symlinks(ensemble_dated, gene_dated, ensemble_output, gene_output)
        
        print(f"✅ Successfully combined count data:")
        print(f"   - Ensemble counts: {ensemble_dated} ({ensemble_df.shape[0]:,} genes, {ensemble_df.shape[1]} samples)")
        print(f"   - Gene symbol counts: {gene_dated} ({gene_df.shape[0]:,} genes, {gene_df.shape[1]} samples)")
        print(f"   - Date: {current_date}")
        
    except Exception as e:
        print(f"❌ Error combining count files: {str(e)}", file=sys.stderr)
        sys.exit(1)

def load_gene_mapping(gtf_file):
    """Load gene ID to gene symbol mapping from GTF file."""
    gene_mapping = {}
    
    with open(gtf_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            if line.startswith('#'):
                continue
            if '\tgene\t' in line:
                # Extract gene_id and gene_name from attributes
                gene_id_match = re.search(r'gene_id "([^"]+)"', line)
                gene_name_match = re.search(r'gene_name "([^"]+)"', line)
                
                if gene_id_match and gene_name_match:
                    gene_id = gene_id_match.group(1)
                    gene_name = gene_name_match.group(1)
                    gene_mapping[gene_id] = gene_name
            
            # Progress indicator for large files
            if line_num % 100000 == 0:
                print(f"   Processed {line_num:,} GTF lines...")
    
    return gene_mapping

def load_count_data(count_files):
    """Load and combine count data from all sample files."""
    combined_data = {}
    sample_names = []
    
    for count_file in count_files:
        # Extract sample name from filename
        sample_name = Path(count_file).stem
        sample_names.append(sample_name)
        
        print(f"   Processing {sample_name}...")
        
        # Load count data efficiently
        try:
            df = pd.read_csv(count_file, sep='\t', header=None, names=['Gene', 'Count'], 
                           dtype={'Gene': str, 'Count': int})
            
            for _, row in df.iterrows():
                gene_id = row['Gene']
                count = int(row['Count']) if pd.notna(row['Count']) else 0
                
                if gene_id not in combined_data:
                    combined_data[gene_id] = {}
                
                combined_data[gene_id][sample_name] = count
                
        except Exception as e:
            print(f"   Warning: Error processing {count_file}: {e}")
            continue
    
    return combined_data, sample_names

def create_count_matrix(combined_data, sample_names):
    """Create count matrix DataFrame."""
    df = pd.DataFrame.from_dict(combined_data, orient='index')
    df.index.name = 'Gene'
    df = df.fillna(0).astype(int)  # Fill missing values and ensure integer counts
    df = df.reindex(columns=sample_names)  # Ensure consistent column order
    return df

def create_gene_symbol_matrix(combined_data, sample_names, gene_mapping):
    """Create gene symbol count matrix by aggregating ENSG IDs."""
    gene_symbol_data = {}
    
    for gene_id, counts in combined_data.items():
        # Map ENSG ID to gene symbol
        gene_symbol = gene_mapping.get(gene_id, gene_id)  # Use original ID if no mapping found
        
        if gene_symbol in gene_symbol_data:
            # If gene symbol already exists, sum the counts (multiple ENSG to same symbol)
            for sample, count in counts.items():
                gene_symbol_data[gene_symbol][sample] = gene_symbol_data[gene_symbol].get(sample, 0) + count
        else:
            gene_symbol_data[gene_symbol] = counts.copy()
    
    # Create DataFrame
    df = pd.DataFrame.from_dict(gene_symbol_data, orient='index')
    df.index.name = 'Gene'
    df = df.fillna(0).astype(int)
    df = df.reindex(columns=sample_names)
    return df

def create_symlinks(ensemble_dated, gene_dated, ensemble_output, gene_output):
    """Create symbolic links for Snakemake expected outputs."""
    # Remove existing symlinks if they exist
    for output_file in [ensemble_output, gene_output]:
        if os.path.exists(output_file) or os.path.islink(output_file):
            os.remove(output_file)
    
    # Create new symlinks
    os.symlink(os.path.basename(ensemble_dated), ensemble_output)
    os.symlink(os.path.basename(gene_dated), gene_output)

# Snakemake execution
if __name__ == "__main__":
    # Use snakemake object inputs and outputs
    count_files = list(snakemake.input.count_files)
    gtf_file = str(snakemake.input.gtf_file)
    
    ensemble_output = str(snakemake.output.ensemble_counts)
    gene_output = str(snakemake.output.gene_counts)
    
    combine_count_files(count_files, gtf_file, ensemble_output, gene_output)