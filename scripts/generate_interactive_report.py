#!/usr/bin/env python3
"""
Optimized interactive report generator for IntegrateALL pipeline.
Creates HTML reports with DataTables, navigation, and interactive elements.
"""

import pandas as pd
import sys
import os
import shutil
from collections import defaultdict
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def safe_read_csv(filepath, delimiter='\t', **kwargs):
    """Safely read CSV file with error handling."""
    try:
        return pd.read_csv(filepath, delimiter=delimiter, **kwargs)
    except Exception as e:
        logger.warning(f"Could not read {filepath}: {e}")
        return pd.DataFrame()


def safe_read_file(filepath):
    """Safely read text file with error handling."""
    try:
        with open(filepath, 'r') as f:
            return f.read()
    except Exception as e:
        logger.warning(f"Could not read {filepath}: {e}")
        return ""


def copy_file_safely(src, dst):
    """Copy file with error handling."""
    try:
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        shutil.copy2(src, dst)
        return True
    except Exception as e:
        logger.warning(f"Could not copy {src} to {dst}: {e}")
        return False


def generate_file_paths(sample_id):
    """Generate all file paths for the report (relative to HTML file)."""
    return {
        'arriba_pdf': f'"fusions/{os.path.basename(f"fusions/{sample_id}.pdf")}"',
        'multiqc': f'"qc/multiqc_report.html"',
        'rnaseqcnv_plot': f'"RNAseqCNV/{sample_id}_CNV_main_fig.png"'
    }


def create_fusion_table(file_path, table_type="arriba"):
    """Create HTML table for fusion data (Arriba or FusionCatcher)."""
    if table_type == "arriba":
        headers = ["5' gene name", "5' chr.position", "3' gene name", "3'chr. position", 
                  "discordant_mates", "confidence", "reading frame", "fusion transcript"]
        table_id = "arriba_table"
    else:  # fusioncatcher
        headers = ["5' gene name", "5' chr.position", "3' gene name", "3'chr. position", 
                  "fusion unique spanning reads", "fusion description", "fusion sequence", "predicted effect"]
        table_id = "fusioncatcher_table"

    table_content = f"<table id='{table_id}' class='my-table-class'>\n"
    table_content += "<thead><tr>" + "".join(f"<th>{h}</th>" for h in headers) + "</tr></thead>\n"
    table_content += "<tbody>\n"

    try:
        with open(file_path, 'r') as file:
            next(file)  # Skip header
            for line in file:
                columns = line.strip().split('\t')
                
                if table_type == "arriba":
                    if len(columns) >= 28:  # Full Arriba format
                        row_data = [columns[0], columns[4], columns[1], columns[5], 
                                   columns[11], columns[14], columns[15], columns[27]]
                    elif len(columns) >= 6:  # Test data format
                        row_data = [columns[0], columns[4] if len(columns) > 4 else 'N/A',
                                   columns[1], columns[5] if len(columns) > 5 else 'N/A',
                                   columns[2] if len(columns) > 2 else 'N/A', 'N/A', 'N/A', 'N/A']
                    else:
                        continue
                else:  # fusioncatcher
                    if len(columns) > 15:
                        row_data = [columns[0], columns[8], columns[1], columns[9], 
                                   columns[5], columns[2], columns[14], columns[15]]
                    else:
                        continue
                        
                table_content += "<tr>" + "".join(f"<td>{data}</td>" for data in row_data) + "</tr>\n"
                
    except Exception as e:
        logger.warning(f"Error reading {table_type} file {file_path}: {e}")

    table_content += "</tbody>\n</table>\n"
    return table_content


def generate_hotspot_tables(hotspot_dir):
    """Generate HTML tables from hotspots directory."""
    if not os.path.exists(hotspot_dir) or not os.listdir(hotspot_dir):
        return "<p>No hotspot found</p>"
    
    table_dict = defaultdict(list)
    
    for filename in os.listdir(hotspot_dir):
        file_path = os.path.join(hotspot_dir, filename)
        
        if filename.endswith('_output_file.csv'):
            table_data = safe_read_csv(file_path, delimiter=',')
            if not table_data.empty:
                html_table = table_data.to_html(classes='my-table-class no-sort', index=False)
                display_name = filename.split('_output_file.csv')[0].replace('_', ': ')
                table_dict[display_name].append(f"<h3>{display_name}</h3>{html_table}")

        elif filename.endswith('_with_bases.tsv'):
            columns = ["chrom", "pos", "ref", "reads_pp", "mismatches_pp", "deletions_pp", 
                      "insertions_pp", "A_pp", "C_pp", "T_pp", "G_pp", "N_pp", "new_base"]
            table_data = safe_read_csv(file_path, delimiter=' ', skiprows=1, names=columns)
            if not table_data.empty:
                html_table = table_data.to_html(classes='my-table-class no-sort', index=False)
                display_name = filename.split('_with_bases.tsv')[0].replace('_', ': ')
                table_dict[display_name + " variants"].append(f"<h3>{display_name} variants</h3>{html_table}")

        elif filename.endswith('_gatk_result.tsv'):
            table_data = safe_read_csv(file_path, delimiter='\t')
            if not table_data.empty:
                html_table = table_data.to_html(classes='my-table-class no-sort', index=False)
                display_name = filename.split('_gatk_result.tsv')[0].replace('_', ': ')
                table_dict[display_name + " GATK Result"].append(f"<h3>{display_name} GATK Result</h3>{html_table}")
    
    return "<br>".join([content for key in sorted(table_dict.keys()) for content in table_dict[key]])


def process_star_log(star_log_file):
    """Process STAR log file into structured HTML tables."""
    try:
        star_log_data = safe_read_csv(star_log_file, delimiter='\t')
        if star_log_data.empty:
            return "<p>STAR log data unavailable</p>", ""
            
        star_log_data.columns = ['Parameter', 'Value']
        
        # Find section breaks
        section_breaks = star_log_data[star_log_data.iloc[:, 0].str.endswith(':')].index
        first_section_end = section_breaks[0] if len(section_breaks) > 0 else len(star_log_data)
        
        # Process first part
        first_part = star_log_data.iloc[:first_section_end, :].copy()
        first_part.iloc[:, 0] = first_part.iloc[:, 0].str.replace('|', '')
        first_part_filtered = first_part[first_part.iloc[:, 0] != 'Started job on ']
        first_part_html = first_part_filtered.to_html(classes='my-table-class no-sort', index=False)
        
        # Process second part into sections
        second_part = star_log_data.iloc[first_section_end:, :]
        grouped_sections = {}
        current_section = None
        
        for _, row in second_part.iterrows():
            text = row.iloc[0]
            if text.endswith(':'):
                current_section = text
                grouped_sections[current_section] = []
            elif current_section is not None:
                grouped_sections[current_section].append(row)
        
        table_segments = ''
        for section, rows in grouped_sections.items():
            if rows:  # Only create table if there are rows
                section_data = pd.DataFrame(rows, columns=second_part.columns)
                table_html = section_data.to_html(classes='my-table-class no-sort', index=False)
                table_segments += f'<h2>{section}</h2>{table_html}'
        
        return first_part_html, table_segments
        
    except Exception as e:
        logger.warning(f"Error processing STAR log: {e}")
        return "<p>STAR log data unavailable</p>", ""


def get_custom_css():
    """Return the custom CSS styling."""
    return """
    <style>
        nav {
            position: fixed;
            top: 0;
            left: 0;
            width: 100%;
            background-color: #333;
            padding: 10px 0;
            z-index: 999;
        }
        
        .navbar {
            display: flex;
            justify-content: space-between;
            align-items: center;
            background-color: #333;
            padding: 10px 20px;
            box-sizing: border-box;
            position: fixed;
            top: 0;
            width: 100%;
            z-index: 1000;
        }
        
        .nav-links {
            display: flex;
            align-items: center;
            flex-grow: 1;
        }
        
        .nav-links a {
            color: white;
            padding: 14px 20px;
            text-decoration: none;
            text-align: center;
        }
        
        .nav-links a:hover {
            background-color: #ddd;
            color: black;
        }
        
        .logo {
            flex-shrink: 0;
        }
        
        .logo img {
            height: 50px;
            max-width: 100%;
        }
                
        body {
            margin: 0;
            padding: 0;
            font-family: Calibri, sans-serif; 
        }
        
        h1::before {
            content: "";
            display: block;
            height: 70px;
            margin-top: -70px;
            visibility: hidden;
        }
                
        .my-table-class {
            font-family: Calibri, sans-serif;
            border-collapse: collapse;
            width: 100%;
        }
        .my-table-class th, .my-table-class td {
            border: 1px solid #ddd;
            padding: 8px;
            text-align: left;
        }
        .my-table-class th {
            background-color: #f2f2f2;
            color: #333;
        }
        .my-table-class tr:nth-child(even) {
            background-color: #f9f9f9;
        }
        .my-table-class tr:hover {
            background-color: #f1f1f1;
        }
        
        .center {
          display: block;
          margin-left: auto;
          margin-right: auto;
          width: 50%;
        } 
    </style>
    """


def generate_report(sample_id, prediction_file, fusioncatcher_file, arriba_file, arriba_pdf, 
                   rnaseqcnv_log2fc_file, rnaseqcnv_plot, rnaseqcnv_manual_table,
                   star_log_file, multiqc_report, final_classification_text, output_file):
    """Generate interactive HTML report."""
    
    logger.info(f"📊 Generating interactive report for sample: {sample_id}")
    
    # Create output directory and copy required assets
    output_dir = os.path.dirname(output_file)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        
        # Copy required assets
        assets_copied = 0
        
        # Copy MultiQC report
        multiqc_dir = os.path.join(output_dir, "qc")
        if copy_file_safely(multiqc_report, os.path.join(multiqc_dir, "multiqc_report.html")):
            assets_copied += 1
        
        # Copy RNAseqCNV plot
        rnaseqcnv_dir = os.path.join(output_dir, "RNAseqCNV")
        if copy_file_safely(rnaseqcnv_plot, os.path.join(rnaseqcnv_dir, os.path.basename(rnaseqcnv_plot))):
            assets_copied += 1
        
        # Copy Arriba PDF
        fusion_dir = os.path.join(output_dir, "fusions")
        if copy_file_safely(arriba_pdf, os.path.join(fusion_dir, os.path.basename(arriba_pdf))):
            assets_copied += 1
        
        # Copy logo if it exists
        logo_src = "scripts/logo.png"
        if os.path.exists(logo_src):
            if copy_file_safely(logo_src, os.path.join(output_dir, "logo.png")):
                assets_copied += 1
        
        logger.info(f"📁 Copied {assets_copied} asset files")
    
    # Generate file paths
    paths = generate_file_paths(sample_id)
    
    # Read and process all data files
    prediction_data = safe_read_csv(prediction_file, delimiter='\t')
    prediction_subset = prediction_data.iloc[:, :14] if not prediction_data.empty else pd.DataFrame()
    cell_origin_subset = prediction_data.iloc[:, 15:21] if prediction_data.shape[1] > 15 else pd.DataFrame()
    
    final_text = safe_read_csv(final_classification_text, delimiter='\t')
    final_classification_html = final_text.to_html(classes='my-table-class no-sort', index=False) if not final_text.empty else f"<p>Final classification for sample {sample_id}</p>"
    
    rnaseq_manual_table = safe_read_csv(rnaseqcnv_manual_table, delimiter='\t')
    rnaseq_manual_html = rnaseq_manual_table.to_html(classes='my-table-class no-sort', index=False) if not rnaseq_manual_table.empty else "<p>RNAseqCNV manual table unavailable</p>"
    
    rnaseq_log2fc = safe_read_csv(rnaseqcnv_log2fc_file, delimiter='\t')
    rnaseq_log2fc_html = rnaseq_log2fc.to_html(classes='my-table-class', index=False) if not rnaseq_log2fc.empty else "<p>RNAseqCNV log2FC data unavailable</p>"
    
    # Process STAR log
    star_first_part, star_segments = process_star_log(star_log_file)
    
    # Create fusion tables
    arriba_table = create_fusion_table(arriba_file, "arriba")
    fusioncatcher_table = create_fusion_table(fusioncatcher_file, "fusioncatcher")
    
    # Create data tables
    prediction_html = prediction_subset.to_html(classes='my-table-class no-sort', index=False) if not prediction_subset.empty else "<p>No prediction data available</p>"
    cell_origin_html = cell_origin_subset.to_html(classes='my-table-class no-sort', index=False) if not cell_origin_subset.empty else "<p>No cell of origin data available</p>"
    
    # Get CSS
    custom_css = get_custom_css()

    # Construct HTML output
    html_output = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>IntegrateALL Report - {sample_id}</title>
        <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.11.5/css/jquery.dataTables.css">
        <script type="text/javascript" charset="utf8" src="https://code.jquery.com/jquery-3.6.0.min.js"></script>
        <script type="text/javascript" charset="utf8" src="https://cdn.datatables.net/1.11.5/js/jquery.dataTables.js"></script>
        {custom_css}
    </head>
    <body>
        <nav class="navbar">
            <div class="nav-links">
                <a href="#section1">Overview</a>
                <a href="#section2">Prediction</a>
                <a href="#section3">MultiQC</a>
                <a href="#section4">STAR Results</a>
                <a href="#section5">RNASeqCNV Plot</a>
                <a href="#section6">RNASeqCNV Table</a>
                <a href="#section7">ARRIBA Fusions</a>
                <a href="#section8">ARRIBA Plots</a>
                <a href="#section9">Fusioncatcher Fusions</a> 
            </div>
            <a class="logo" href="https://www.catchall-kfo5010.com/">
                <img src="logo.png" alt="Logo">
            </a>
        </nav>
        <h1>IntegrateALL</h1>
        <h1 id="section1">IntegrateALL classification</h1>
        {final_classification_html}
        
        <h1 id="section2">ALLCatchR Prediction Data</h1>
        {prediction_html}
        <h1 id="section2">ALLCatchR Cell of origin</h1>
        {cell_origin_html}
        
        <h1 id="section3">MultiQC Results</h1>
        <iframe src={paths['multiqc']} name="multiqc_path" width="100%" height="600" frameborder="0"></iframe>
        
        <h1 id="section4">STAR Results</h1>
        {star_first_part}
        {star_segments}
        
        <h1 id="section5">RNASeq-CNV Plot</h1>
        <img src={paths['rnaseqcnv_plot']} alt="RNASeqCNV plot" width="900" height="600" class="center">
        
        <h2 id="section6">RNASeq-CNV Manual An Table</h2>
        {rnaseq_manual_html}
        
        <h2 id="section6">RNASeq-CNV Log2FC Table</h2>        
        {rnaseq_log2fc_html}

        <h1 id="section7">ARRIBA Fusions Table</h1>
        {arriba_table}
        
        <h1 id="section8">ARRIBA Fusions Plots</h1>
        <iframe src={paths['arriba_pdf']} name="arriba_fusion_pdf_path" width="100%" height="600" frameborder="0"></iframe>
        
        <h1 id="section9">Fusioncatcher Fusions</h1>
        {fusioncatcher_table}

        <script>
            $(document).ready(function() {{
                $('.my-table-class').not('.no-sort').DataTable();
            }});
        </script>
    </body>
    </html>
    """

    # Write to HTML file
    output_dir = os.path.dirname(output_file)
    if output_dir:  # Only create directory if there is a directory path
        os.makedirs(output_dir, exist_ok=True)
    with open(output_file, 'w') as file:
        file.write(html_output)
    
    logger.info(f"✅ Interactive report generated: {output_file}")


if __name__ == "__main__":
    if len(sys.argv) != 13:
        print("Usage: python generate_interactive_report.py <sample_id> <prediction_file> "
              "<fusioncatcher_file> <arriba_file> <arriba_pdf> <rnaseqcnv_log2fc_file> "
              "<rnaseqcnv_plot> <rnaseqcnv_manual_table> <star_log_file> <multiqc_report> "
              "<final_classification_text> <output_html>")
        sys.exit(1)
    
    args = sys.argv[1:]
    generate_report(*args)