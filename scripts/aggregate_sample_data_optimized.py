#!/usr/bin/env python3
"""
Optimized sample data aggregation script for IntegrateALL pipeline.
Processes multiple input files and generates a comprehensive sample report.
Replaces the inefficient bash-based process_data.sh script.
"""

import pandas as pd
import sys
import os
import logging
from pathlib import Path

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def read_comparison_file(comparison_file):
    """Read the comparison file content."""
    try:
        with open(comparison_file, 'r') as f:
            return f.read().strip()
    except Exception as e:
        logger.warning(f"Could not read comparison file {comparison_file}: {e}")
        return ""


def extract_star_uniquely_mapped_reads(star_log_file):
    """Extract uniquely mapped reads count from STAR log file."""
    try:
        with open(star_log_file, 'r') as f:
            for line in f:
                if 'Uniquely mapped reads number' in line:
                    parts = line.strip().split('\t')
                    if len(parts) >= 2:
                        return parts[1].strip()
        return "N/A"
    except Exception as e:
        logger.warning(f"Error reading STAR log file {star_log_file}: {e}")
        return "N/A"


def read_fastqc_data(fastqc_file):
    """Read FastQC data from MultiQC output."""
    try:
        # This reads the MultiQC FastQC summary data
        # The exact format depends on MultiQC output structure
        with open(fastqc_file, 'r') as f:
            content = f.read().strip()
            if content:
                return content
            else:
                return "FastQC data available in MultiQC report"
    except Exception as e:
        logger.warning(f"Could not read FastQC file {fastqc_file}: {e}")
        return "FastQC data not available"


def read_allcatchr_predictions(prediction_file):
    """Read ALLCatchR prediction results."""
    try:
        df = pd.read_csv(prediction_file, sep='\t')
        if df.empty:
            return []
        
        results = []
        for _, row in df.iterrows():
            results.append({
                'prediction': row.get('Prediction', 'N/A'),
                'score': row.get('Score', 'N/A'), 
                'confidence': row.get('Confidence', 'N/A')
            })
        return results
    except Exception as e:
        logger.warning(f"Error reading ALLCatchR predictions {prediction_file}: {e}")
        return []


def read_fusioncatcher_data(fusioncatcher_file):
    """Read FusionCatcher fusion data."""
    try:
        df = pd.read_csv(fusioncatcher_file, sep='\t', skiprows=1, header=None)
        if df.empty:
            return []
            
        fusions = []
        for _, row in df.iterrows():
            if len(row) >= 10:
                fusions.append({
                    'gene_5p': row[0],
                    'pos_5p': row[8] if len(row) > 8 else 'N/A',
                    'gene_3p': row[1], 
                    'pos_3p': row[9] if len(row) > 9 else 'N/A',
                    'spanning_reads': row[5] if len(row) > 5 else 'N/A'
                })
        return fusions
    except Exception as e:
        logger.warning(f"Error reading FusionCatcher data {fusioncatcher_file}: {e}")
        return []


def read_arriba_data(arriba_file):
    """Read Arriba fusion data."""
    try:
        df = pd.read_csv(arriba_file, sep='\t', skiprows=1, header=None)
        if df.empty:
            return []
            
        fusions = []
        for _, row in df.iterrows():
            if len(row) >= 12:
                fusions.append({
                    'gene_5p': row[0],
                    'pos_5p': row[4] if len(row) > 4 else 'N/A',
                    'gene_3p': row[1],
                    'pos_3p': row[5] if len(row) > 5 else 'N/A', 
                    'discordant_mates': row[11] if len(row) > 11 else 'N/A'
                })
        return fusions
    except Exception as e:
        logger.warning(f"Error reading Arriba data {arriba_file}: {e}")
        return []


def read_rnaseqcnv_karyotype_data(manual_table_file):
    """Read RNAseqCNV karyotype analysis data."""
    try:
        df = pd.read_csv(manual_table_file, sep='\t')
        if df.empty:
            return []
            
        karyotypes = []
        for _, row in df.iterrows():
            karyotypes.append({
                'gender': row.get('gender', 'N/A'),
                'chrom_n': row.get('chrom_n', 'N/A'),
                'alterations': row.get('alterations', 'N/A'),
                'status': row.get('status', 'N/A'),
                'comments': row.get('comments', 'N/A')
            })
        return karyotypes
    except Exception as e:
        logger.warning(f"Error reading RNAseqCNV karyotype data {manual_table_file}: {e}")
        return []


def read_rnaseqcnv_arm_calls(log2fc_file):
    """Read RNAseqCNV chromosome arm calls."""
    try:
        df = pd.read_csv(log2fc_file, sep='\t')
        if df.empty:
            return []
            
        arm_calls = []
        for _, row in df.iterrows():
            arm_calls.append({
                'chr': row.iloc[0] if len(row) > 0 else 'N/A',
                'arm': row.iloc[1] if len(row) > 1 else 'N/A', 
                'log2fc': row.iloc[2] if len(row) > 2 else 'N/A'
            })
        return arm_calls
    except Exception as e:
        logger.warning(f"Error reading RNAseqCNV arm calls {log2fc_file}: {e}")
        return []


def generate_aggregated_report(sample_id, comparison_content, uniquely_mapped_reads,
                              fastqc_left, fastqc_right, allcatchr_data, 
                              fusioncatcher_data, arriba_data, karyotype_data, 
                              arm_calls_data):
    """Generate the aggregated report content."""
    
    report_lines = []
    
    # Add comparison file content first
    if comparison_content:
        report_lines.append(comparison_content)
    
    # Add sequencing summary
    report_lines.append(f"The transcriptome sequencing of {sample_id} produced {uniquely_mapped_reads} uniquely aligned sequencing reads, enabling quantification of protein coding genes.")
    
    # Add quality metrics
    report_lines.append("\\nQuality metrics (fastQC / MultiQC) indicated:")
    report_lines.append(fastqc_left)
    report_lines.append(fastqc_right)
    
    # Add ALLCatchR subtype allocation
    report_lines.append(f"\\nALLCatchR allocated for sample {sample_id} the following molecular subtype:")
    report_lines.append("subtype prediction\tscore\tconfidence")
    
    for prediction in allcatchr_data:
        report_lines.append(f"{prediction['prediction']}\t{prediction['score']}\t{prediction['confidence']}")
    
    # Add fusion results
    report_lines.append("\\nfusioncatcher / ARRIBA identified the following driver fusion candidates:")
    
    # FusionCatcher results
    report_lines.append("\\nFusioncatcher:")
    report_lines.append("5' gene name\t5' chr.position\t3' gene name\t3'chr. position\tfusion unique spanning reads")
    
    for fusion in fusioncatcher_data:
        report_lines.append(f"{fusion['gene_5p']}\t{fusion['pos_5p']}\t{fusion['gene_3p']}\t{fusion['pos_3p']}\t{fusion['spanning_reads']}")
    
    # Arriba results
    report_lines.append("\\nARRIBA:")
    report_lines.append("5' gene name\t5' chr.position\t3' gene name\t3'chr. position\tdiscordant_mates")
    
    for fusion in arriba_data:
        report_lines.append(f"{fusion['gene_5p']}\t{fusion['pos_5p']}\t{fusion['gene_3p']}\t{fusion['pos_3p']}\t{fusion['discordant_mates']}")
    
    # RNASeqCNV karyotype results
    report_lines.append("\\nRNASeqCNV identified the following karyotype:")
    report_lines.append("gender\tchrom_n\talterations\tstatus\tcomments")
    
    for karyo in karyotype_data:
        report_lines.append(f"{karyo['gender']}\t{karyo['chrom_n']}\t{karyo['alterations']}\t{karyo['status']}\t{karyo['comments']}")
    
    # Chromosome arm calls
    report_lines.append("Chromosome arm calls")
    report_lines.append("chr\tarm\tlog2FC per arm")
    
    for arm_call in arm_calls_data:
        report_lines.append(f"{arm_call['chr']}\t{arm_call['arm']}\t{arm_call['log2fc']}")
    
    return "\\n".join(report_lines)


def main(sample_id, star_log_file, fastqc_left_file, fastqc_right_file, 
         prediction_file, fusioncatcher_file, arriba_file, manual_table_file,
         log2fc_file, comparison_file, output_file):
    """
    Main function to aggregate sample data.
    """
    
    logger.info(f"📊 Aggregating data for sample: {sample_id}")
    
    # Create output directory
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    # Read all input data efficiently
    logger.info("📖 Reading input data files...")
    
    comparison_content = read_comparison_file(comparison_file)
    uniquely_mapped_reads = extract_star_uniquely_mapped_reads(star_log_file)
    fastqc_left = read_fastqc_data(fastqc_left_file)
    fastqc_right = read_fastqc_data(fastqc_right_file)
    allcatchr_data = read_allcatchr_predictions(prediction_file)
    fusioncatcher_data = read_fusioncatcher_data(fusioncatcher_file)
    arriba_data = read_arriba_data(arriba_file)
    karyotype_data = read_rnaseqcnv_karyotype_data(manual_table_file)
    arm_calls_data = read_rnaseqcnv_arm_calls(log2fc_file)
    
    # Generate aggregated report
    report_content = generate_aggregated_report(
        sample_id, comparison_content, uniquely_mapped_reads,
        fastqc_left, fastqc_right, allcatchr_data,
        fusioncatcher_data, arriba_data, karyotype_data, arm_calls_data
    )
    
    # Write output file
    with open(output_file, 'w') as f:
        f.write(report_content)
    
    logger.info(f"✅ Aggregated report generated: {output_file}")
    logger.info(f"📈 Processed {len(allcatchr_data)} ALLCatchR predictions, "
               f"{len(fusioncatcher_data)} FusionCatcher fusions, "
               f"{len(arriba_data)} Arriba fusions")


if __name__ == "__main__":
    if len(sys.argv) != 12:
        print("Usage: python aggregate_sample_data_optimized.py <sample_id> <star_log_file> "
              "<fastqc_left_file> <fastqc_right_file> <prediction_file> <fusioncatcher_file> "
              "<arriba_file> <manual_table_file> <log2fc_file> <comparison_file> <output_file>")
        sys.exit(1)
    
    args = sys.argv[1:]
    main(*args)