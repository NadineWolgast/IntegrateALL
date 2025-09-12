#!/usr/bin/env python3
"""
Optimized script to aggregate final classification results across all samples.
Combines txt, curation, and report files with dated outputs.
"""

import pandas as pd
import datetime
import os
import sys
from pathlib import Path

def aggregate_classification_files(txt_files, curation_files, report_files, 
                                 txt_output, curation_output, report_output):
    """
    Aggregate classification files with sample IDs and create dated outputs.
    
    Args:
        txt_files: List of txt classification files
        curation_files: List of curation files  
        report_files: List of report files
        txt_output: Output path for aggregated txt file
        curation_output: Output path for aggregated curation file
        report_output: Output path for aggregated report file
    """
    
    print(f"📊 Aggregating classification results from {len(txt_files)} samples")
    
    # Generate current date string
    current_date = datetime.datetime.now().strftime("%Y%m%d")
    
    try:
        # Create output directory
        os.makedirs(os.path.dirname(txt_output), exist_ok=True)
        
        # Create dated output filenames
        txt_dated = f"Final_classification/Aggregated_output_txt_{current_date}.csv"
        curation_dated = f"Final_classification/Aggregated_output_curation_{current_date}.csv"
        report_dated = f"Final_classification/Aggregated_output_report_{current_date}.csv"
        
        # Process txt files
        print("Aggregating txt classification files...")
        aggregate_files(txt_files, txt_dated, '_output_txt.csv')
        
        # Process curation files  
        print("Aggregating curation files...")
        aggregate_files(curation_files, curation_dated, '_curation.csv')
        
        # Process report files
        print("Aggregating report files...")
        aggregate_files(report_files, report_dated, '_output_report.csv')
        
        # Create symbolic links for Snakemake outputs
        create_symlinks(txt_dated, curation_dated, report_dated,
                       txt_output, curation_output, report_output)
        
        print(f"✅ Successfully aggregated classification results:")
        print(f"   - Txt file: {txt_dated}")
        print(f"   - Curation file: {curation_dated}")
        print(f"   - Report file: {report_dated}")
        print(f"   - Date: {current_date}")
        
    except Exception as e:
        print(f"❌ Error aggregating classification files: {str(e)}", file=sys.stderr)
        sys.exit(1)

def aggregate_files(file_list, output_path, suffix_to_remove):
    """Aggregate a list of CSV files with sample IDs."""
    if not file_list:
        # Create empty file if no input files
        Path(output_path).touch()
        return
    
    with open(output_path, 'w') as out_file:
        # Write header with Sample_ID column
        with open(file_list[0], 'r') as first_file:
            header = first_file.readline().strip()
            out_file.write(f"Sample_ID,{header}\n")
        
        # Process each file
        for file_path in file_list:
            sample_id = Path(file_path).stem.replace(suffix_to_remove.replace('.csv', ''), '')
            
            with open(file_path, 'r') as in_file:
                next(in_file)  # Skip header
                for line in in_file:
                    line = line.strip()
                    if line:  # Only write non-empty lines
                        out_file.write(f"{sample_id},{line}\n")

def create_symlinks(txt_dated, curation_dated, report_dated,
                   txt_output, curation_output, report_output):
    """Create symbolic links for Snakemake expected outputs."""
    # Remove existing symlinks if they exist
    for output_file in [txt_output, curation_output, report_output]:
        if os.path.exists(output_file) or os.path.islink(output_file):
            os.remove(output_file)
    
    # Create new symlinks
    os.symlink(os.path.basename(txt_dated), txt_output)
    os.symlink(os.path.basename(curation_dated), curation_output)
    os.symlink(os.path.basename(report_dated), report_output)

# Snakemake execution
if __name__ == "__main__":
    # Use snakemake object inputs and outputs
    txt_files = list(snakemake.input.txt_files)
    curation_files = list(snakemake.input.curation_files)  
    report_files = list(snakemake.input.report_files)
    
    txt_output = str(snakemake.output.txt_agg)
    curation_output = str(snakemake.output.curation_agg)
    report_output = str(snakemake.output.report_agg)
    
    aggregate_classification_files(txt_files, curation_files, report_files,
                                 txt_output, curation_output, report_output)