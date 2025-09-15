#!/usr/bin/env python3
"""
Optimized interactive report generator for IntegrateALL pipeline.
Creates HTML reports by processing multiple input files and copying required assets.
"""

import os
import sys
import shutil
import pandas as pd
from pathlib import Path
import logging
from datetime import datetime

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def copy_file_safely(src, dst):
    """Copy file with error handling."""
    try:
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        shutil.copy2(src, dst)
        return True
    except Exception as e:
        logger.warning(f"Could not copy {src} to {dst}: {e}")
        return False


def read_star_log(star_log_file):
    """Extract uniquely mapped reads from STAR log file."""
    try:
        with open(star_log_file, 'r') as f:
            for line in f:
                if 'Uniquely mapped reads number' in line:
                    parts = line.strip().split('\t')
                    if len(parts) >= 2:
                        return parts[1].strip()
        return "N/A"
    except Exception as e:
        logger.warning(f"Error reading STAR log: {e}")
        return "N/A"


def read_allcatchr_predictions(prediction_file):
    """Read ALLCatchR predictions efficiently."""
    try:
        df = pd.read_csv(prediction_file, sep='\t')
        if df.empty:
            return "N/A", "N/A", "N/A"
        row = df.iloc[0]
        return row.get('Prediction', 'N/A'), row.get('Score', 'N/A'), row.get('Confidence', 'N/A')
    except Exception as e:
        logger.warning(f"Error reading ALLCatchR predictions: {e}")
        return "N/A", "N/A", "N/A"


def read_fusioncatcher_data(fusioncatcher_file):
    """Read FusionCatcher data efficiently."""
    fusion_data = []
    try:
        df = pd.read_csv(fusioncatcher_file, sep='\t', skiprows=1, header=None)
        for _, row in df.iterrows():
            if len(row) >= 10:
                fusion_data.append({
                    'gene_5p': row[0],
                    'pos_5p': row[8] if len(row) > 8 else 'N/A',
                    'gene_3p': row[1],
                    'pos_3p': row[9] if len(row) > 9 else 'N/A',
                    'spanning_reads': row[5] if len(row) > 5 else 'N/A'
                })
    except Exception as e:
        logger.warning(f"Error reading FusionCatcher data: {e}")
    
    return fusion_data


def read_arriba_data(arriba_file):
    """Read Arriba data efficiently."""
    fusion_data = []
    try:
        df = pd.read_csv(arriba_file, sep='\t', skiprows=1, header=None)
        for _, row in df.iterrows():
            if len(row) >= 12:
                fusion_data.append({
                    'gene_5p': row[0],
                    'pos_5p': row[4] if len(row) > 4 else 'N/A',
                    'gene_3p': row[1],
                    'pos_3p': row[5] if len(row) > 5 else 'N/A',
                    'discordant_mates': row[11] if len(row) > 11 else 'N/A'
                })
    except Exception as e:
        logger.warning(f"Error reading Arriba data: {e}")
    
    return fusion_data


def read_rnaseqcnv_data(manual_table_file, log2fc_file):
    """Read RNAseqCNV data efficiently."""
    karyotype_data = []
    arm_calls_data = []
    
    # Read karyotype data
    try:
        df = pd.read_csv(manual_table_file, sep='\t')
        for _, row in df.iterrows():
            karyotype_data.append({
                'gender': row.get('gender', 'N/A'),
                'chrom_n': row.get('chrom_n', 'N/A'),
                'alterations': row.get('alterations', 'N/A'),
                'status': row.get('status', 'N/A'),
                'comments': row.get('comments', 'N/A')
            })
    except Exception as e:
        logger.warning(f"Error reading RNAseqCNV karyotype data: {e}")
    
    # Read arm calls data
    try:
        df = pd.read_csv(log2fc_file, sep='\t')
        for _, row in df.iterrows():
            arm_calls_data.append({
                'chr': row.iloc[0] if len(row) > 0 else 'N/A',
                'arm': row.iloc[1] if len(row) > 1 else 'N/A',
                'log2fc': row.iloc[2] if len(row) > 2 else 'N/A'
            })
    except Exception as e:
        logger.warning(f"Error reading RNAseqCNV arm calls data: {e}")
    
    return karyotype_data, arm_calls_data


def generate_report_content(sample_id, star_uniquely_mapped, allcatchr_data, 
                          fusioncatcher_data, arriba_data, karyotype_data, arm_calls_data,
                          final_classification_text):
    """Generate the main report content."""
    
    subtype, score, confidence = allcatchr_data
    
    content = []
    
    # Add final classification content
    try:
        with open(final_classification_text, 'r') as f:
            content.append(f.read().strip())
    except Exception as e:
        logger.warning(f"Could not read final classification file: {e}")
        content.append(f"Final classification analysis for sample {sample_id}")
    
    # Add sequencing summary
    content.append(f"\nThe transcriptome sequencing of {sample_id} produced {star_uniquely_mapped} uniquely aligned sequencing reads, enabling quantification of protein coding genes.")
    
    # Quality metrics placeholder
    content.append("\nQuality metrics (fastQC / MultiQC) indicated:")
    # Note: FastQC content would need to be extracted from MultiQC report
    
    # ALLCatchR results
    content.append(f"\nALLCatchR allocated for sample {sample_id} the following molecular subtype:")
    content.append("subtype prediction\tscore\tconfidence")
    content.append(f"{subtype}\t{score}\t{confidence}")
    
    # Fusion results
    content.append("\nfusioncatcher / ARRIBA identified the following driver fusion candidates:")
    
    # FusionCatcher results
    content.append("\nFusioncatcher:")
    content.append("5' gene name\t5' chr.position\t3' gene name\t3'chr. position\tfusion unique spanning reads")
    for fusion in fusioncatcher_data:
        content.append(f"{fusion['gene_5p']}\t{fusion['pos_5p']}\t{fusion['gene_3p']}\t{fusion['pos_3p']}\t{fusion['spanning_reads']}")
    
    # Arriba results
    content.append("\nARRIBA:")
    content.append("5' gene name\t5' chr.position\t3' gene name\t3'chr. position\tdiscordant_mates")
    for fusion in arriba_data:
        content.append(f"{fusion['gene_5p']}\t{fusion['pos_5p']}\t{fusion['gene_3p']}\t{fusion['pos_3p']}\t{fusion['discordant_mates']}")
    
    # RNASeqCNV results
    content.append("\nRNASeqCNV identified the following karyotype:")
    content.append("gender\tchrom_n\talterations\tstatus\tcomments")
    for karyo in karyotype_data:
        content.append(f"{karyo['gender']}\t{karyo['chrom_n']}\t{karyo['alterations']}\t{karyo['status']}\t{karyo['comments']}")
    
    content.append("Chromosome arm calls")
    content.append("chr\tarm\tlog2FC per arm")
    for arm_call in arm_calls_data:
        content.append(f"{arm_call['chr']}\t{arm_call['arm']}\t{arm_call['log2fc']}")
    
    return "\n".join(content)


def generate_html_report(content, sample_id, output_file):
    """Generate HTML report with embedded content."""
    
    html_template = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>IntegrateALL Report - {sample_id}</title>
        <style>
            body {{ 
                font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
                margin: 40px;
                background-color: #f8f9fa;
                color: #333;
            }}
            .header {{ 
                text-align: center;
                padding: 30px 0;
                background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
                color: white;
                border-radius: 10px;
                margin-bottom: 30px;
                box-shadow: 0 4px 15px rgba(0,0,0,0.1);
            }}
            .header h1 {{ 
                margin: 0;
                font-size: 2.5em;
                font-weight: 300;
            }}
            .header h2 {{ 
                margin: 10px 0 0 0;
                font-size: 1.2em;
                opacity: 0.9;
            }}
            .content {{ 
                background: white;
                padding: 30px;
                border-radius: 10px;
                box-shadow: 0 2px 10px rgba(0,0,0,0.1);
                white-space: pre-wrap;
                font-family: 'Courier New', monospace;
                line-height: 1.6;
            }}
            .timestamp {{
                text-align: center;
                margin-top: 20px;
                color: #666;
                font-size: 0.9em;
            }}
        </style>
    </head>
    <body>
        <div class="header">
            <h1>🧬 IntegrateALL Report</h1>
            <h2>Sample: {sample_id}</h2>
        </div>
        
        <div class="content">
{content}
        </div>
        
        <div class="timestamp">
            Report generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
        </div>
    </body>
    </html>
    """
    
    with open(output_file, 'w') as f:
        f.write(html_template)


def main(sample_id, prediction_file, fusioncatcher_file, arriba_file, arriba_pdf,
         rnaseqcnv_log2fc_file, rnaseqcnv_plot, rnaseqcnv_manual_table, 
         star_log_file, multiqc_report, final_classification_text, output_html):
    """
    Main function to generate interactive report.
    """
    
    logger.info(f"📊 Generating interactive report for sample: {sample_id}")
    
    # Create output directory
    output_dir = os.path.dirname(output_html)
    os.makedirs(output_dir, exist_ok=True)
    
    # Copy required assets
    assets_copied = 0
    
    # Copy fusion PDF
    fusion_dir = os.path.join(output_dir, "fusions")
    if copy_file_safely(arriba_pdf, os.path.join(fusion_dir, os.path.basename(arriba_pdf))):
        assets_copied += 1
    
    # Copy MultiQC report
    qc_dir = os.path.join(output_dir, "qc")
    if copy_file_safely(multiqc_report, os.path.join(qc_dir, os.path.basename(multiqc_report))):
        assets_copied += 1
    
    # Copy RNAseqCNV plot
    rnaseqcnv_dir = os.path.join(output_dir, "RNAseqCNV")
    if copy_file_safely(rnaseqcnv_plot, os.path.join(rnaseqcnv_dir, os.path.basename(rnaseqcnv_plot))):
        assets_copied += 1
    
    # Copy logo if it exists
    logo_src = "scripts/logo.png"
    if os.path.exists(logo_src):
        if copy_file_safely(logo_src, os.path.join(output_dir, "logo.png")):
            assets_copied += 1
    
    logger.info(f"📁 Copied {assets_copied} asset files")
    
    # Read all input data
    logger.info("📖 Reading input data files...")
    
    star_uniquely_mapped = read_star_log(star_log_file)
    allcatchr_data = read_allcatchr_predictions(prediction_file)
    fusioncatcher_data = read_fusioncatcher_data(fusioncatcher_file)
    arriba_data = read_arriba_data(arriba_file)
    karyotype_data, arm_calls_data = read_rnaseqcnv_data(rnaseqcnv_manual_table, rnaseqcnv_log2fc_file)
    
    # Generate report content
    content = generate_report_content(
        sample_id, star_uniquely_mapped, allcatchr_data,
        fusioncatcher_data, arriba_data, karyotype_data, arm_calls_data,
        final_classification_text
    )
    
    # Generate HTML report
    generate_html_report(content, sample_id, output_html)
    
    logger.info(f"✅ Interactive report generated: {output_html}")
    logger.info(f"🔗 Report includes {len(fusioncatcher_data)} FusionCatcher and {len(arriba_data)} Arriba fusions")


if __name__ == "__main__":
    if len(sys.argv) != 12:
        print("Usage: python generate_interactive_report.py <sample_id> <prediction_file> "
              "<fusioncatcher_file> <arriba_file> <arriba_pdf> <rnaseqcnv_log2fc_file> "
              "<rnaseqcnv_plot> <rnaseqcnv_manual_table> <star_log_file> <multiqc_report> "
              "<final_classification_text> <output_html>")
        sys.exit(1)
    
    args = sys.argv[1:]
    main(*args)