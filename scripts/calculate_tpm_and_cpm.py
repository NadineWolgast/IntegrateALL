#!/usr/bin/env python3
"""
Optimized script to calculate TPM and CPM from STAR ReadsPerGene output.
Replaces R script with faster Python/pandas implementation.
"""

import pandas as pd
import numpy as np
import sys
import os

def calculate_tpm_and_cpm(reads_per_gene_file, tpm_output, cpm_output, cds_length_file="data/annotation/cds_length.tsv"):
    """
    Calculate TPM (Transcripts Per Million) and CPM (Counts Per Million) values.
    
    Args:
        reads_per_gene_file: STAR ReadsPerGene.out.tab file
        tpm_output: Output file for TPM values
        cds_length_file: File with gene lengths for TPM calculation
        cpm_output: Output file for CPM values
    """
    
    print(f"📊 Processing: {reads_per_gene_file}")
    
    try:
        # Read STAR ReadsPerGene output (skip first 4 lines, use first 2 columns)
        print("Reading gene counts...")
        counts = pd.read_csv(reads_per_gene_file, sep='\t', skiprows=4, 
                           usecols=[0, 1], names=['Gene', 'Count'], 
                           dtype={'Gene': str, 'Count': int})
        
        # Read CDS length data for TPM calculation
        print("Reading CDS length data...")
        cds_lengths = pd.read_csv(cds_length_file, sep='\t')
        
        # Match genes with their lengths
        print("Matching genes with CDS lengths...")
        merged_data = counts.merge(cds_lengths, left_on='Gene', right_on='ESNG', how='left')
        
        # Fill missing CDS lengths with median
        median_length = merged_data['cds_length'].median()
        merged_data = merged_data.copy()  # Avoid pandas chained assignment warning
        merged_data['cds_length'].fillna(median_length, inplace=True)
        
        print(f"Using median CDS length ({median_length:.0f}) for {merged_data['cds_length'].isna().sum()} missing entries")
        
        # Calculate TPM (Transcripts Per Million)
        print("Calculating TPM values...")
        # TPM = (count / gene_length) * 1e6 / sum(count / gene_length for all genes)
        rpk = merged_data['Count'] / merged_data['cds_length']  # Reads per kilobase
        scaling_factor = rpk.sum() / 1e6
        tpm = rpk / scaling_factor
        
        # Calculate CPM (Counts Per Million)  
        print("Calculating CPM values...")
        total_counts = merged_data['Count'].sum()
        cpm = (merged_data['Count'] / total_counts) * 1e6
        
        # Create output directories
        os.makedirs(os.path.dirname(tpm_output), exist_ok=True)
        os.makedirs(os.path.dirname(cpm_output), exist_ok=True)
        
        # Write TPM output
        tpm_df = pd.DataFrame({
            'Gene': merged_data['Gene'], 
            'tpm': tpm
        })
        tpm_df.to_csv(tpm_output, sep='\t', index=False, float_format='%.6f')
        
        # Write CPM output
        cpm_df = pd.DataFrame({
            'Gene': merged_data['Gene'],
            'CPM': cpm
        })
        cpm_df.to_csv(cpm_output, sep='\t', index=False, float_format='%.6f')
        
        print(f"✅ Successfully calculated:")
        print(f"   - TPM for {len(tpm_df):,} genes → {tpm_output}")
        print(f"   - CPM for {len(cpm_df):,} genes → {cpm_output}")
        print(f"   - Total counts processed: {total_counts:,}")
        
    except Exception as e:
        print(f"❌ Error calculating TPM/CPM: {str(e)}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python calculate_tpm_and_cpm.py <reads_per_gene_file> <tpm_output> <cpm_output>")
        sys.exit(1)
    
    reads_per_gene_file = sys.argv[1]
    tpm_output = sys.argv[2] 
    cpm_output = sys.argv[3]
    
    if not os.path.exists(reads_per_gene_file):
        print(f"❌ Input file not found: {reads_per_gene_file}", file=sys.stderr)
        sys.exit(1)
    
    calculate_tpm_and_cpm(reads_per_gene_file, tpm_output, cpm_output)