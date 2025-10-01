#!/usr/bin/env python3
"""
Fusion Intersect Analysis Script

This script compares fusion calls from Arriba and FusionCatcher to identify
intersecting fusions, including those reported in reverse orientation.
The goal is to identify novel/non-established fusions for subtype analysis.

Author: Nadine Wolgast
"""

import pandas as pd
import sys
import os
import logging
from pathlib import Path

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def normalize_gene_pair(gene1, gene2):
    """
    Normalize gene pair to handle reverse orientations.
    Returns both orientations as a tuple for comparison.
    """
    # Remove @ symbol from IGH genes for comparison
    gene1_norm = gene1.replace('@', '')
    gene2_norm = gene2.replace('@', '')
    
    # Return both orientations for matching
    orientation1 = (gene1_norm, gene2_norm)
    orientation2 = (gene2_norm, gene1_norm)
    
    return orientation1, orientation2


def parse_arriba_file(arriba_file):
    """
    Parse Arriba fusion file and extract gene pairs with spanning reads.
    """
    logger.info(f"Parsing Arriba file: {arriba_file}")
    
    if not os.path.exists(arriba_file) or os.path.getsize(arriba_file) == 0:
        logger.warning(f"Arriba file is empty or doesn't exist: {arriba_file}")
        return {}
    
    try:
        # Read Arriba file and find header line (first non-comment line)
        with open(arriba_file, 'r') as f:
            lines = f.readlines()
        
        # Find first non-comment line (header)
        header_idx = 0
        for i, line in enumerate(lines):
            if not line.startswith('#'):
                header_idx = i
                break
        
        # Read from header line
        df = pd.read_csv(arriba_file, sep='\t', skiprows=header_idx)
        if df.empty:
            logger.warning("Arriba file is empty")
            return {}
        
        arriba_fusions = {}
        
        for _, row in df.iterrows():
            # Parse gene names from Arriba format
            fusion_desc = str(row.iloc[0]).strip()  # Full fusion description
            gene1_col = str(row.iloc[1]).strip()  # gene1 column
            
            # Arriba has two formats:
            # 1. "GENE1(id),GENE2(id)" in first column, gene1 in second column
            # 2. "GENE1" in first column, "GENE2" in second column
            
            if ',' in fusion_desc and '(' in fusion_desc:
                # Format 1: Extract both genes from first column
                parts = fusion_desc.split(',')
                if len(parts) >= 2:
                    gene1 = parts[0].split('(')[0].strip()
                    gene2_raw = parts[1].strip()
                    gene2 = gene2_raw.split('(')[0].strip()
                else:
                    continue  # Skip if we can't parse
            else:
                # Format 2: Gene1 in first column, Gene2 in second column  
                gene1 = fusion_desc
                gene2 = gene1_col
            
            # Get confidence level and filter out low confidence fusions
            confidence = str(row.iloc[14]).strip().lower() if len(row) > 14 else ""
            if confidence == "low":
                continue
            
            # Get spanning reads (split_reads1 + split_reads2)
            split_reads1 = int(row.iloc[9]) if len(row) > 9 else 0  # split_reads1 column
            split_reads2 = int(row.iloc[10]) if len(row) > 10 else 0  # split_reads2 column
            spanning_reads = split_reads1 + split_reads2
            
            # Store fusion with original orientation
            fusion_key = f"{gene1}::{gene2}"
            arriba_fusions[fusion_key] = {
                'gene1': gene1,
                'gene2': gene2,
                'spanning_reads': spanning_reads,
                'source': 'Arriba'
            }
        
        logger.info(f"Found {len(arriba_fusions)} fusions in Arriba")
        return arriba_fusions
        
    except Exception as e:
        logger.error(f"Error parsing Arriba file: {e}")
        return {}


def parse_fusioncatcher_file(fc_file):
    """
    Parse FusionCatcher fusion file and extract gene pairs with spanning reads.
    """
    logger.info(f"Parsing FusionCatcher file: {fc_file}")
    
    if not os.path.exists(fc_file) or os.path.getsize(fc_file) == 0:
        logger.warning(f"FusionCatcher file is empty or doesn't exist: {fc_file}")
        return {}
    
    try:
        # Read FusionCatcher file (tab-separated, skip first line which is header)
        df = pd.read_csv(fc_file, sep='\t', skiprows=1, header=None)
        if df.empty:
            logger.warning("FusionCatcher file is empty")
            return {}
        
        # Define column indices based on FusionCatcher format
        fc_fusions = {}
        
        for _, row in df.iterrows():
            if len(row) < 6:  # Ensure we have enough columns
                continue
                
            gene1 = str(row[0]).strip()  # Gene_1_symbol(5end_fusion_partner)
            gene2 = str(row[1]).strip()  # Gene_2_symbol(3end_fusion_partner)
            fusion_description = str(row[2]).strip()  # Fusion_description
            spanning_reads = int(row[5])  # Spanning_unique_reads
            
            # Skip fusions that are marked as "banned" in the description
            if 'banned' in fusion_description.lower():
                continue
            
            # Store fusion with original orientation
            fusion_key = f"{gene1}::{gene2}"
            fc_fusions[fusion_key] = {
                'gene1': gene1,
                'gene2': gene2,
                'spanning_reads': spanning_reads,
                'source': 'FusionCatcher'
            }
        
        logger.info(f"Found {len(fc_fusions)} fusions in FusionCatcher")
        return fc_fusions
        
    except Exception as e:
        logger.error(f"Error parsing FusionCatcher file: {e}")
        return {}


def find_fusion_intersects(arriba_fusions, fc_fusions):
    """
    Find intersecting fusions between Arriba and FusionCatcher.
    Handles both same and reverse orientations.
    """
    logger.info("Finding fusion intersects...")
    
    intersects = []
    
    # Create normalized lookup for FusionCatcher fusions
    fc_normalized = {}
    for fc_key, fc_data in fc_fusions.items():
        gene1, gene2 = fc_data['gene1'], fc_data['gene2']
        norm1, norm2 = normalize_gene_pair(gene1, gene2)
        fc_normalized[norm1] = fc_data
        fc_normalized[norm2] = fc_data
    
    # Check each Arriba fusion for matches
    for arriba_key, arriba_data in arriba_fusions.items():
        gene1, gene2 = arriba_data['gene1'], arriba_data['gene2']
        norm1, norm2 = normalize_gene_pair(gene1, gene2)
        
        # Check if this fusion exists in FusionCatcher (either orientation)
        fc_match = None
        orientation = "same"
        
        if norm1 in fc_normalized:
            fc_match = fc_normalized[norm1]
            if fc_match['gene1'] != gene1 or fc_match['gene2'] != gene2:
                orientation = "reverse"
        elif norm2 in fc_normalized:
            fc_match = fc_normalized[norm2]
            orientation = "reverse"
        
        if fc_match:
            intersect = {
                'arriba_gene1': arriba_data['gene1'],
                'arriba_gene2': arriba_data['gene2'],
                'arriba_spanning_reads': arriba_data['spanning_reads'],
                'fc_gene1': fc_match['gene1'],
                'fc_gene2': fc_match['gene2'],
                'fc_spanning_reads': fc_match['spanning_reads'],
                'orientation': orientation,
                'fusion_pair': f"{gene1}::{gene2}"  # Use Arriba orientation as reference
            }
            intersects.append(intersect)
    
    logger.info(f"Found {len(intersects)} intersecting fusions")
    return intersects


def write_intersect_results(intersects, output_file, sample_id):
    """
    Write intersecting fusions to CSV file.
    """
    logger.info(f"Writing results to: {output_file}")
    
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    # Convert to DataFrame for easy CSV writing
    if intersects:
        df = pd.DataFrame(intersects)
        
        # Add sample ID column
        df['sample_id'] = sample_id
        
        # Reorder columns
        df = df[['sample_id', 'fusion_pair', 'arriba_gene1', 'arriba_gene2', 
                'arriba_spanning_reads', 'fc_gene1', 'fc_gene2', 
                'fc_spanning_reads', 'orientation']]
        
        # Sort by total spanning reads (descending)
        df['total_spanning_reads'] = df['arriba_spanning_reads'] + df['fc_spanning_reads']
        df = df.sort_values('total_spanning_reads', ascending=False)
        df = df.drop('total_spanning_reads', axis=1)
        
        df.to_csv(output_file, index=False)
        logger.info(f"Wrote {len(intersects)} intersecting fusions to {output_file}")
    else:
        # Create empty file with headers
        empty_df = pd.DataFrame(columns=['sample_id', 'fusion_pair', 'arriba_gene1', 'arriba_gene2', 
                                       'arriba_spanning_reads', 'fc_gene1', 'fc_gene2', 
                                       'fc_spanning_reads', 'orientation'])
        empty_df.to_csv(output_file, index=False)
        logger.info(f"No intersecting fusions found. Created empty file: {output_file}")


def main():
    """
    Main function to run fusion intersect analysis.
    """
    if len(sys.argv) != 5:
        print("Usage: python fusion_intersect_analysis.py <sample_id> <arriba_file> <fusioncatcher_file> <output_file>")
        sys.exit(1)
    
    sample_id = sys.argv[1]
    arriba_file = sys.argv[2]
    fusioncatcher_file = sys.argv[3]
    output_file = sys.argv[4]
    
    logger.info(f"Starting fusion intersect analysis for sample: {sample_id}")
    
    # Parse input files
    arriba_fusions = parse_arriba_file(arriba_file)
    fc_fusions = parse_fusioncatcher_file(fusioncatcher_file)
    
    # Find intersects
    intersects = find_fusion_intersects(arriba_fusions, fc_fusions)
    
    # Write results
    write_intersect_results(intersects, output_file, sample_id)
    
    logger.info("Fusion intersect analysis completed successfully")


if __name__ == "__main__":
    main()