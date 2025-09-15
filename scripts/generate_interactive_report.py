#!/usr/bin/env python3
"""
Interactive report generator for IntegrateALL pipeline.
Creates HTML reports with DataTables, navigation, and interactive elements.
"""

import pandas as pd
import sys
import os
from collections import defaultdict


def read_file_content(filename):
    """Read content from file."""
    with open(filename, 'r') as file:
        return file.read()


def write_to_file(filename, content):
    """Write content to file."""
    with open(filename, 'w') as file:
        file.write(content)


def extract_uniquely_mapped_reads(filename):
    """Extract uniquely mapped reads from STAR log file."""
    with open(filename, 'r') as file:
        for line in file:
            if 'Uniquely mapped reads number' in line:
                return line.split('\t')[1].strip()
    return ''


def generate_pdf_path(sample_id):
    """Generate path to Arriba fusion PDF."""
    arriba_fusion_pdf_path = f'"fusions/{sample_id}.pdf"'
    return arriba_fusion_pdf_path


def generate_multiqc_path(sample_id):
    """Generate path to MultiQC report."""
    multiqc_path = f'"qc/multiqc/{sample_id}/multiqc_report.html"'
    return multiqc_path


def generate_RNASeq_cnv_png_path(sample_id):
    """Generate path to RNASeqCNV plot."""
    RNASeq_cnv_png_path = f'"RNAseqCNV_output/gatk/{sample_id}_gatk/{sample_id}/{sample_id}_CNV_main_fig.png"'
    return RNASeq_cnv_png_path


def arriba_file_to_html_table(arriba_file):
    """Convert Arriba file to HTML table."""
    table_content = "<table id='arriba_table'>"
    table_content += "<table border='1' class='my-table-class'>\\n"
    table_content += ("<thead><tr><th>5' gene name</th><th>5' chr.position</th><th>3' gene name</th><th>3'chr. "
                      "position</th><th>discordant_mates</th><th>confidence</th><th>reading frame</th><th>fusion "
                      "transcript</th></tr></thead>\\n")

    table_content += "<tbody>\\n"

    try:
        with open(arriba_file, 'r') as file:
            next(file)  # Skip header
            for line in file:
                columns = line.strip().split('\\t')
                if len(columns) >= 28:  # Full Arriba format
                    table_content += "<tr>"
                    table_content += f"<td>{columns[0]}</td><td>{columns[4]}</td><td>{columns[1]}</td><td>{columns[5]}</td><td>{columns[11]}</td><td>{columns[14]}</td><td>{columns[15]}</td><td>{columns[27]}</td>"
                    table_content += "</tr>\\n"
                elif len(columns) >= 6:  # Test data format
                    table_content += "<tr>"
                    table_content += f"<td>{columns[0]}</td><td>{columns[4] if len(columns) > 4 else 'N/A'}</td><td>{columns[1]}</td><td>{columns[5] if len(columns) > 5 else 'N/A'}</td><td>{columns[2] if len(columns) > 2 else 'N/A'}</td><td>N/A</td><td>N/A</td><td>N/A</td>"
                    table_content += "</tr>\\n"
    except Exception as e:
        print(f"Error reading Arriba file: {e}")

    table_content += "</tbody>\\n"
    table_content += "</table>\\n"
    return table_content


def fusionctacher_file_to_html_table(fusioncatcher_file):
    """Convert FusionCatcher file to HTML table."""
    table_content = "<table id='fusioncatcher_table' class='my-table-class'>\\n"
    table_content += ("<thead><tr><th>5' gene name</th><th>5' chr.position</th><th>3' gene name</th><th>3'chr. "
                      "position</th><th>fusion unique spanning reads</th><th>fusion description</th>"
                      "<th>fusion sequence</th><th>predicted effect</th></tr></thead>\\n")
    table_content += "<tbody>\\n"

    try:
        with open(fusioncatcher_file, 'r') as file:
            next(file)  # Skip header
            for line in file:
                columns = line.strip().split('\\t')
                if len(columns) > 15:
                    table_content += "<tr>"
                    table_content += f"<td>{columns[0]}</td><td>{columns[8]}</td><td>{columns[1]}</td><td>{columns[9]}</td><td>{columns[5]}</td><td>{columns[2]}</td><td>{columns[14]}</td><td>{columns[15]}</td>"
                    table_content += "</tr>\\n"
    except Exception as e:
        print(f"Error reading FusionCatcher file: {e}")

    table_content += "</tbody>\\n"
    table_content += "</table>\\n"
    return table_content


def generate_html_table(directory):
    """Generate HTML tables from hotspots directory."""
    table_dict = defaultdict(list)
    html_content = ""
    
    if os.path.exists(directory) and os.listdir(directory):
        sorted_content = ""
        files = os.listdir(directory)
        if len(files) > 1:
            for filename in os.listdir(directory):
                if filename.endswith('_output_file.csv'):
                    file_path = os.path.join(directory, filename)
                    table_data = pd.read_csv(file_path, delimiter=',')
                    html_table = table_data.to_html(classes='my-table-class no-sort', index=False)
                    display_name = filename.split('_output_file.csv')[0].replace('_', ': ')
                    html_content += f"<h3>{display_name}</h3>"
                    html_content += html_table
                    table_dict[display_name].append(f"<h3>{display_name}</h3>" + html_table)

                elif filename.endswith('_with_bases.tsv'):
                    file_path = os.path.join(directory, filename)
                    columns = ["chrom", "pos", "ref", "reads_pp", "mismatches_pp", "deletions_pp", "insertions_pp", "A_pp",
                               "C_pp", "T_pp", "G_pp", "N_pp", "new_base"]
                    table_data = pd.read_csv(file_path, delimiter=' ', skiprows=1, names=columns)
                    html_table = table_data.to_html(classes='my-table-class no-sort', index=False)
                    display_name = filename.split('_with_bases.tsv')[0].replace('_', ': ')
                    html_content += f"<h3>{display_name} variants</h3>"
                    html_content += html_table
                    table_dict[display_name + " variants"].append(f"<h3>{display_name} variants</h3>" + html_table)

                elif filename.endswith('_gatk_result.tsv'):
                    file_path = os.path.join(directory, filename)
                    table_data = pd.read_csv(file_path, delimiter='\\t')
                    html_table = table_data.to_html(classes='my-table-class no-sort', index=False)
                    display_name = filename.split('_gatk_result.tsv')[0].replace('_', ': ')
                    html_content += f"<h3>{display_name} GATK Result</h3>"
                    html_content += html_table
                    table_dict[display_name + " GATK Result"].append(f"<h3>{display_name} GATK Result</h3>" + html_table)
        else:
            sorted_content = "<p>No hotspot found</p>"
    
    for key in sorted(table_dict.keys()):
        sorted_content += "<br>".join(table_dict[key])
    return sorted_content


def generate_report(sample_id, prediction_file, fusioncatcher_file, arriba_file, arriba_pdf, 
                   rnaseqcnv_log2fc_file, rnaseqcnv_plot, rnaseqcnv_manual_table,
                   star_log_file, multiqc_report, final_classification_text, output_file):
    """Generate interactive HTML report."""
    
    # Read all input files
    try:
        prediction_data = pd.read_csv(prediction_file, delimiter='\\t')
        prediction_data_subset = prediction_data.iloc[:, :14]
        cell_of_origin_data_subset = prediction_data.iloc[:, 15:21] if prediction_data.shape[1] > 15 else pd.DataFrame()
    except Exception as e:
        print(f"Error reading prediction file: {e}")
        prediction_data_subset = pd.DataFrame()
        cell_of_origin_data_subset = pd.DataFrame()

    try:
        final_text = pd.read_csv(final_classification_text, delimiter='\\t')
        html_table_text = final_text.to_html(classes='my-table-class no-sort', index=False)
    except Exception as e:
        print(f"Error reading final classification: {e}")
        html_table_text = f"<p>Final classification for sample {sample_id}</p>"

    try:
        star_log_data = pd.read_csv(star_log_file, delimiter='\\t')
        star_log_data.columns = ['Parameter', 'Value']
        
        # Split STAR log into sections
        first_section_end = star_log_data[star_log_data.iloc[:, 0].str.endswith(':')].index[0] if not star_log_data[star_log_data.iloc[:, 0].str.endswith(':')].empty else len(star_log_data)
        first_part = star_log_data.iloc[:first_section_end, :]
        first_part.loc[:, first_part.columns[0]] = first_part[first_part.columns[0]].str.replace('|', '')
        second_part = star_log_data.iloc[first_section_end:, :]
        
        first_part_without_start = first_part[first_part[first_part.columns[0]] != 'Started job on |']
        first_part_html = first_part_without_start.to_html(classes='my-table-class no-sort', index=False)
        
        # Group second part by sections
        grouped_sections = {}
        current_section = None
        
        for index, row in second_part.iterrows():
            text = row[second_part.columns[0]]
            if text.endswith(':'):
                current_section = text
                grouped_sections[current_section] = []
            elif current_section is not None:
                grouped_sections[current_section].append(row)
        
        table_segments = ''
        for section, rows in grouped_sections.items():
            section_data = pd.DataFrame(rows, columns=second_part.columns)
            table_html = section_data.to_html(classes='my-table-class no-sort', index=False)
            table_segments += f'<h2>{section}</h2>'
            table_segments += table_html
            
    except Exception as e:
        print(f"Error reading STAR log: {e}")
        first_part_html = "<p>STAR log data unavailable</p>"
        table_segments = ""

    try:
        rna_seq_cnv_manual_an_table_file = pd.read_csv(rnaseqcnv_manual_table, delimiter='\\t')
        rna_seq_cnv_manual_an_table_file_html_table = rna_seq_cnv_manual_an_table_file.to_html(classes='my-table-class no-sort', index=False)
    except Exception as e:
        print(f"Error reading RNAseqCNV manual table: {e}")
        rna_seq_cnv_manual_an_table_file_html_table = "<p>RNAseqCNV manual table unavailable</p>"

    try:
        rna_seq_cnv_log2foldchange_file = pd.read_csv(rnaseqcnv_log2fc_file, delimiter='\\t')
        rna_seq_cnv_log2foldchange_file_html_table = rna_seq_cnv_log2foldchange_file.to_html(classes='my-table-class', index=False)
    except Exception as e:
        print(f"Error reading RNAseqCNV log2fc file: {e}")
        rna_seq_cnv_log2foldchange_file_html_table = "<p>RNAseqCNV log2FC data unavailable</p>"

    # Generate path references
    arriba_fusion_pdf_path = generate_pdf_path(sample_id)
    rnaseqcnv_png_path = generate_RNASeq_cnv_png_path(sample_id)
    multiqc_path = generate_multiqc_path(sample_id)
    
    # Generate tables
    arriba_html_table = arriba_file_to_html_table(arriba_file)
    fusioncatcher_html_table = fusionctacher_file_to_html_table(fusioncatcher_file)
    
    # Create HTML tables for data
    html_table = prediction_data_subset.to_html(classes='my-table-class no-sort', index=False) if not prediction_data_subset.empty else "<p>No prediction data available</p>"
    cell_of_origin_table = cell_of_origin_data_subset.to_html(classes='my-table-class no-sort', index=False) if not cell_of_origin_data_subset.empty else "<p>No cell of origin data available</p>"

    # Custom CSS for styling
    custom_css = """
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

    html_table_with_style = custom_css + html_table
    cell_of_origin_table_with_style = custom_css + cell_of_origin_table
    arriba_html_table_with_style = custom_css + arriba_html_table

    # Construct HTML output
    html_output = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>IntegrateALL Report - {sample_id}</title>
        <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.11.5/css/jquery.dataTables.css">
        <script type="text/javascript" charset="utf8" src="https://code.jquery.com/jquery-3.6.0.min.js"></script>
        <script type="text/javascript" charset="utf8" src="https://cdn.datatables.net/1.11.5/js/jquery.dataTables.js"></script>
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
        {html_table_text}
        
        <h1 id="section2">ALLCatchR Prediction Data</h1>
        {html_table_with_style}
        <h1 id="section2">ALLCatchR Cell of origin</h1>
        {cell_of_origin_table_with_style}
        
        <h1 id="section3">MultiQC Results</h1>
        <iframe src={multiqc_path} name="multiqc_path" width="100%" height="600" frameborder="0"></iframe>
        
        <h1 id="section4">STAR Results</h1>
        {first_part_html}
        {table_segments}
        
        <h1 id="section5">RNASeq-CNV Plot</h1>
        <img src={rnaseqcnv_png_path} alt="RNASeqCNV plot" width="900" height="600" class="center">
        
        <h2 id="section6">RNASeq-CNV Manual An Table</h2>
        {rna_seq_cnv_manual_an_table_file_html_table}
        
        <h2 id="section6">RNASeq-CNV Log2FC Table</h2>        
        {rna_seq_cnv_log2foldchange_file_html_table}

        <h1 id="section7">ARRIBA Fusions Table</h1>
        {arriba_html_table_with_style}
        
        <h1 id="section8">ARRIBA Fusions Plots</h1>
        <iframe src={arriba_fusion_pdf_path} name="arriba_fusion_pdf_path" width="100%" height="600" frameborder="0"></iframe>
        
        <h1 id="section9">Fusioncatcher Fusions</h1>
        {fusioncatcher_html_table}

        <script>
            $(document).ready(function() {{
                $('.my-table-class').not('.no-sort').DataTable();
            }});
        </script>
    </body>
    </html>
    """

    # Write to HTML file
    with open(output_file, 'w') as file:
        file.write(html_output)


if __name__ == "__main__":
    if len(sys.argv) != 13:
        print("Usage: python generate_interactive_report.py <sample_id> <prediction_file> "
              "<fusioncatcher_file> <arriba_file> <arriba_pdf> <rnaseqcnv_log2fc_file> "
              "<rnaseqcnv_plot> <rnaseqcnv_manual_table> <star_log_file> <multiqc_report> "
              "<final_classification_text> <output_html>")
        sys.exit(1)
    
    args = sys.argv[1:]
    generate_report(*args)