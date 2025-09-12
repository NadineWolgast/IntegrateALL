#!/usr/bin/env python3
"""
Optimized script to process STAR ReadsPerGene.out.tab files.
Replaces awk with faster Python/pandas implementation.
"""

import pandas as pd
import sys
import os

def process_reads_per_gene(reads_per_gene_file, output_file, skip_rows=4):
    """
    Process STAR ReadsPerGene output, extracting gene names and counts.
    
    Args:
        reads_per_gene_file: STAR ReadsPerGene.out.tab file
        output_file: Output file for processed counts
        skip_rows: Number of rows to skip (default: 4)
    """
    
    try:
        # Read STAR ReadsPerGene output (skip header lines, use first 2 columns)
        counts = pd.read_csv(reads_per_gene_file, sep='\t', skiprows=skip_rows, 
                           usecols=[0, 1], names=['Gene', 'Count'], 
                           dtype={'Gene': str, 'Count': int})
        
        # Create output directory
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        
        # Write output with tab separation
        counts.to_csv(output_file, sep='\t', index=False, header=False)
        
        print(f"✅ Processed {len(counts):,} genes from {reads_per_gene_file}")
        
    except Exception as e:
        print(f"❌ Error processing {reads_per_gene_file}: {str(e)}", file=sys.stderr)
        sys.exit(1)

# Snakemake execution
if __name__ == "__main__":
    # Use snakemake object inputs and outputs
    reads_per_gene_file = str(snakemake.input.reads_per_gene)
    output_file = str(snakemake.output.counts)
    
    if not os.path.exists(reads_per_gene_file):
        print(f"❌ Input file not found: {reads_per_gene_file}", file=sys.stderr)
        sys.exit(1)
    
    process_reads_per_gene(reads_per_gene_file, output_file)