import pandas as pd
import sys
import os


def read_file_content(filename):
    with open(filename, 'r') as file:
        return file.read()


def write_to_file(filename, content):
    with open(filename, 'w') as file:
        file.write(content)


def extract_uniquely_mapped_reads(filename):
    with open(filename, 'r') as file:
        for line in file:
            if 'Uniquely mapped reads number' in line:
                return line.split('\t')[1].strip()
    return ''


def generate_pdf_path(sample_id):
    arriba_fusion_pdf_path = f'"fusions/{sample_id}.pdf"'
    return arriba_fusion_pdf_path


def generate_multiqc_right_path():
    multiqc_right_path = f'"multiqc_right/multiqc_report.html"'
    return multiqc_right_path


def generate_multiqc_left_path():
    multiqc_left_path = f'"multiqc_left/multiqc_report.html"'
    return multiqc_left_path


def generate_RNASeq_cnv_png_path(sample_id):
    RNASeq_cnv_png_path = f'"RNAseqCNV/{sample_id}_CNV_main_fig.png"'
    return RNASeq_cnv_png_path


def generate_Hotspot_path(sample_id):
    Hotspot_path = f'"hotspots/"'
    return Hotspot_path

def arriba_file_to_html_table(arriba_file):
    table_content = "<table id='arriba_table'>"
    table_content += "<table border='1' class='my-table-class'>\n"
    table_content += ("<thead><tr><th>5’ gene name</th><th>5’ chr.position</th><th>3’ gene name</th><th>3’chr. "
                      "position</th><th>discordant_mates</th><th>confidence</th><th>reading frame</th><th>fusion "
                      "transcript</th></tr></thead>\n")  # Header der Tabelle

    table_content += "<tbody>\n"

    with open(arriba_file, 'r') as file:
        next(file)
        for line in file:
            columns = line.strip().split('\t')
            table_content += "<tr>"
            table_content += f"<td>{columns[0]}</td><td>{columns[4]}</td><td>{columns[1]}</td><td>{columns[5]}</td><td>{columns[11]}</td><td>{columns[14]}</td><td>{columns[15]}</td><td>{columns[27]}</td>"
            table_content += "</tr>\n"

    table_content += "</tbody>\n"
    table_content += "</table>\n"
    return table_content


def fusionctacher_file_to_html_table(fusioncatcher_file):
    table_content = "<table id='fusioncatcher_table' class='my-table-class'>\n"
    table_content += ("<thead><tr><th>5’ gene name</th><th>5’ chr.position</th><th>3’ gene name</th><th>3’chr. "
                      "position</th><th>fusion unique spanning reads</th><th>fusion description</th>"
                      "<th>fusion sequence</th><th>predicted effect</th></tr></thead>\n")
    table_content += "<tbody>\n"

    with open(fusioncatcher_file, 'r') as file:
        next(file)
        for line in file:
            columns = line.strip().split('\t')
            table_content += "<tr>"
            table_content += f"<td>{columns[0]}</td><td>{columns[8]}</td><td>{columns[1]}</td><td>{columns[9]}</td><td>{columns[5]}</td><td>{columns[2]}</td><td>{columns[14]}</td><td>{columns[15]}</td>"
            table_content += "</tr>\n"

    table_content += "</tbody>\n"
    table_content += "</table>\n"
    return table_content


def generate_html_table(directory):
    html_content = ""
    print(os.getcwd())
    # Überprüfe, ob das Verzeichnis existiert und ob Dateien darin enthalten sind
    if os.path.exists(directory) and os.listdir(directory):
       # print("directory", directory)
        #html_content += "<h2>Pysamstats Hotspots</h2>"

        # Iteriere über die Dateien im Verzeichnis
        for filename in os.listdir(directory):
            # print("filename", filename)
            # Filtere die Dateien mit der Endung '_output:file.csv' oder '_with_bases.csv'
            if filename.endswith('_output_file.csv'):
                print("in if")
                file_path = os.path.join(directory, filename)

                # Lese die CSV-Datei
                table_data = pd.read_csv(file_path, delimiter=',')
                print("table_data", table_data)
                # Konvertiere die Daten in eine HTML-Tabelle
                html_table = table_data.to_html(classes='my-table-class no-sort', index=False)
                display_name = filename.split('_output_file.csv')[0].replace('_', ': ')
                # Füge die HTML-Tabelle zum HTML-Inhalt hinzu
                html_content += f"<h3>{display_name}</h3>"
                html_content += html_table

            elif filename.endswith('_with_bases.tsv'):
                print("in elif")
                file_path = os.path.join(directory, filename)
                columns = ["chrom", "pos", "ref", "reads_pp", "mismatches_pp", "deletions_pp", "insertions_pp", "A_pp",
                           "C_pp", "T_pp", "G_pp", "N_pp", "new_base"]

                # Lese die CSV-Datei
                table_data = pd.read_csv(file_path, delimiter=' ', skiprows=1, names=columns)
                print("_with_bases table_data", table_data)
                # Konvertiere die Daten in eine HTML-Tabelle
                html_table = table_data.to_html(classes='my-table-class no-sort', index=False)
                display_name = filename.split('_with_bases.tsv')[0].replace('_', ': ')
                # Füge die HTML-Tabelle zum HTML-Inhalt hinzu
                html_content += f"<h3>{display_name} variants</h3>"
                html_content += html_table

    return html_content


def generate_report(prediction_file, fusioncatcher_file, arriba_file,
                    rna_seq_cnv_log2foldchange_file, rna_seq_cnv_manual_an_table_file,
                    star_log_final_out_file,
                    comparison_file, pysamstats_files_IKZF1, pysamstats_files_PAX5,
                    pysamstats_files_coverage, hotspots, sample_id, output_file):
    # Read CSV/TSV files
    prediction_data = pd.read_csv(prediction_file, delimiter='\t')  # Example for TSV
    pysamstats_files_IKZF1 = pd.read_csv(pysamstats_files_IKZF1, delimiter='\t')
    pysamstats_files_PAX5 = pd.read_csv(pysamstats_files_PAX5, delimiter='\t')
    pysamstats_files_coverage = pd.read_csv(pysamstats_files_coverage, delimiter='\t')
    star_log_final_out_file = pd.read_csv(star_log_final_out_file, delimiter='\t')

    #print("hotspots", hotspots)
    hotspot_table = generate_html_table(hotspots)
    #print("hotspot table", hotspot_table)
    prediction_data_subset = prediction_data.iloc[:, :9]
    arriba_fusion_pdf_path = generate_pdf_path(sample_id)
    rnaseqcnv_png_path = generate_RNASeq_cnv_png_path(sample_id)
    multiqc_left_path = generate_multiqc_left_path()
    multiqc_right_path = generate_multiqc_right_path()
    arriba_html_table = arriba_file_to_html_table(arriba_file)
    fusioncatcher_html_table = fusionctacher_file_to_html_table(fusioncatcher_file)
    comparison_file = pd.read_csv(comparison_file, delimiter='\t', skiprows=1)
    comparison_html_table = comparison_file.to_html(classes='my-table-class no-sort', index=False)

    rna_seq_cnv_manual_an_table_file = pd.read_csv(rna_seq_cnv_manual_an_table_file, delimiter='\t')
    rna_seq_cnv_manual_an_table_file_html_table = rna_seq_cnv_manual_an_table_file.to_html(classes='my-table-class no-sort',
                                                                                           index=False)
    rna_seq_cnv_log2foldchange_file = pd.read_csv(rna_seq_cnv_log2foldchange_file, delimiter='\t')
    rna_seq_cnv_log2foldchange_file_html_table = rna_seq_cnv_log2foldchange_file.to_html(classes='my-table-class',
                                                                                         index=False)

    html_table = prediction_data_subset.to_html(classes='my-table-class no-sort', index=False)
    pysamstats_IKZF1_html_table = pysamstats_files_IKZF1.to_html(classes='my-table-class no-sort', index=False)
    pysamstats_PAX5_html_table = pysamstats_files_PAX5.to_html(classes='my-table-class no-sort', index=False)
    pysamstats_coverage_html_table = pysamstats_files_coverage.to_html(classes='my-table-class no-sort', index=False)

    first_section_end = star_log_final_out_file[star_log_final_out_file.iloc[:, 0].str.endswith(':')].index[0]
    first_part = star_log_final_out_file.iloc[:first_section_end, :]
    first_part.loc[:, first_part.columns[0]] = first_part[first_part.columns[0]].str.replace('|', '')
    star_log_final_out_file.columns = ['Parameter', 'Value']
    second_part = star_log_final_out_file.iloc[first_section_end:, :]
    second_part.loc[:, second_part.columns[0]] = second_part[second_part.columns[0]].str.replace('|', '')
    first_part_without_start = first_part[first_part[first_part.columns[0]] != 'Started job on |']
    first_part_html = first_part_without_start.to_html(classes='my-table-class no-sort', index=False)

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

        nav a {
            display: inline-block;
            padding: 10px 20px;
            text-decoration: none;
            color: #fff;
            transition: background-color 0.3s ease;
        }

        nav a:hover {
            background-color: #555;
        }

        body {
            margin-top: 60px; 
            padding-top: 20px; 
            font-family: Calibri, sans-serif; 
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
    arriba_html_table_with_style = custom_css + arriba_html_table
    pysamstats_IKZF1_html_table_with_style = custom_css + pysamstats_IKZF1_html_table
    pysamstats_PAX5_html_table_with_style = custom_css + pysamstats_PAX5_html_table
    pysamstats_coverage_html_table_with_style = custom_css + pysamstats_coverage_html_table

    # Construct HTML output
    html_output = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>Report</title>
        <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.11.5/css/jquery.dataTables.css">
        <script type="text/javascript" charset="utf8" src="https://code.jquery.com/jquery-3.6.0.min.js"></script>
        <script type="text/javascript" charset="utf8" src="https://cdn.datatables.net/1.11.5/js/jquery.dataTables.js"></script>
    </head>
    </head>
    <body>
        <nav>
            <a href="#section1">Overview</a>
            <a href="#section2">Prediction</a>
            <a href="#section3">MultiQC Right</a>
            <a href="#section4">MultiQC Left</a>
            <a href="#section5">STAR Results</a>
            <a href="#section6">Pysamstats</a>
            <a href="#section7">RNASeqCNV Plot</a>
            <a href="#section8">RNASeqCNV Table</a>
            <a href="#section10">ARRIBA Fusions</a>
            <a href="#section11">ARRIBA Plots</a>
            <a href="#section12">Fusioncatcher Fusions</a>            
        </nav>

        
        <h1 id="section1">Blast-o-Matic-Fusioninator Results</h1>
        {comparison_html_table}
        
        <h1 id="section2">ALLCatchR Prediction Data</h1>
        {html_table_with_style}
        
        <h1 id="section3">MultiQC Right FASTQ Results</h1>
        <iframe src={multiqc_right_path} name="multiqc_right_path" width="100%" height="600" frameborder="0"></iframe>
        
        <h1 id="section4">MultiQC Left FASTQ Results</h1>
        <iframe src={multiqc_left_path} name="multiqc_left_path" width="100%" height="600" frameborder="0"></iframe>
        
        <h1 id="section5">STAR Results</h1>
        {first_part_html}
        {table_segments}
        
   
        <h1 id="section6">Pysamstats Hotspots</h1>
        <h2>Es wurden 34 Positionen abgefragt mit folgenden Ergebnissen:</h2>
        {hotspot_table}
        <h2> IKZF1</h2>
        {pysamstats_IKZF1_html_table_with_style}
        
        <h2>PAX5</h2>
        {pysamstats_PAX5_html_table_with_style}
        
        <h2">Coverage IKZF1</h2>
        {pysamstats_coverage_html_table_with_style}
        
        <h1 id="section7">RNASeq-CNV Plot</h1>
        <img src={rnaseqcnv_png_path} alt="RNASeqCNV first plot" width="900" height="600" class="center">
        
        <h2 id="section8">RNASeq-CNV Manual An Table</h2>
        {rna_seq_cnv_manual_an_table_file_html_table}
        
        <h2 id="section9">RNASeq-CNV Log2FC Table</h2>        
        {rna_seq_cnv_log2foldchange_file_html_table}

        <h1 id="section10">ARRIBA Fusions Table</h1>
        {arriba_html_table_with_style}
        
        <h1 id="section11">ARRIBA Fusions Plots</h1>
        <iframe src= {arriba_fusion_pdf_path} name="arriba_fusion_pdf_path" width="100%" height="600" frameborder="0"></iframe>
        
        <h1 id="section12">Fusioncatcher Fusions</h1>
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
    if len(sys.argv) != 14:
        print(sys.argv)
        print(len(sys.argv))
        print("Usage: python generate_report.py <prediction_file> <fusioncatcher_file> ... <output_file>")
        sys.exit(1)

    prediction_file, fusioncatcher_file, arriba_file, rna_seq_cnv_log2foldchange_file, \
        rna_seq_cnv_manual_an_table_file, star_log_final_out_file, \
        comparison_file, pysamstats_files_IKZF1, pysamstats_files_PAX5, \
        pysamstats_files_coverage, hotspots, sample_id, output_file = sys.argv[1:]

    generate_report(prediction_file, fusioncatcher_file, arriba_file,
                    rna_seq_cnv_log2foldchange_file, rna_seq_cnv_manual_an_table_file,
                    star_log_final_out_file,
                    comparison_file, pysamstats_files_IKZF1, pysamstats_files_PAX5,
                    pysamstats_files_coverage, hotspots, sample_id, output_file)
