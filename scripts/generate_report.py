import pandas as pd
import sys
from bs4 import BeautifulSoup

def generate_report(prediction_file, fusioncatcher_file, arriba_file,
                    rna_seq_cnv_log2foldchange_file, rna_seq_cnv_manual_an_table_file,
                    star_log_final_out_file, multiqc_fqc_right, multiqc_fqc_left,
                    aggregated_output, output_file):
    # Read CSV/TSV files
    prediction_data = pd.read_csv(prediction_file, delimiter='\t')  # Example for TSV
    aggregated_output = pd.read_csv(aggregated_output, skiprows=17, delimiter='\t', na_values=['NaN', 'N/A', ''])  # Example for CSV
    aggregated_output = aggregated_output.fillna('')
    prediction_data_subset = prediction_data.iloc[:, :9]
    # Read the HTML file
    with open(multiqc_fqc_right, 'r') as file:
        html_content = file.read()

    # Create a BeautifulSoup object
    soup = BeautifulSoup(html_content, 'html.parser')

    # Find and extract content within a specific <div> element with a class or ID
    general_stats = soup.find('div', id='general_stats')  # Replace with your class name
    # OR
    mqc = soup.find('div', id='fastqc-status-check-heatmap')  # Replace with your ID
    extracted_html = ''
    # Extract the text or content within the specific <div> element
    extracted_general_stats = general_stats.encode_contents().decode() if general_stats else ''
    extracted_mqc = mqc.encode_contents().decode() if mqc else ''


    # Construct HTML output
    html_output = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>Report</title>
         <style>
            /* Define styles for your table class */
            .my-table-class {{
                font-family: Arial, sans-serif;
                border-collapse: collapse;
                width: 100%;
            }}
            .my-table-class th, .my-table-class td {{
                border: 1px solid #ddd;
                padding: 8px;
                text-align: left;
            }}
            .my-table-class th {{
                background-color: #f2f2f2;
                color: #333;
            }}
            .my-table-class tr:nth-child(even) {{
                background-color: #f9f9f9;
            }}
            .my-table-class tr:hover {{
                background-color: #f1f1f1;
            }}
        </style>
    </head>
    <body>
        <h1>Prediction Data</h1>
        {prediction_data_subset.to_html()}


        <h1>MultiQC Results</h1>
        {extracted_general_stats} 
        {extracted_mqc} 
        
        <h1>Aggregated Output</h1>
        {aggregated_output.to_html()}
    </body>
    </html>
    """

    # Write to HTML file
    with open(output_file, 'w') as file:
        file.write(html_output)


# Usage: python generate_report.py <prediction_file> <fusioncatcher_file> ... <output_file>
if __name__ == "__main__":
    if len(sys.argv) != 11:
        print(sys.argv)
        print(len(sys.argv))
        print("Usage: python generate_report.py <prediction_file> <fusioncatcher_file> ... <output_file>")
        sys.exit(1)

    prediction_file, fusioncatcher_file, arriba_file, rna_seq_cnv_log2foldchange_file, \
        rna_seq_cnv_manual_an_table_file, star_log_final_out_file, multiqc_fqc_right, \
        multiqc_fqc_left, comparison_file, output_file = sys.argv[1:]

    generate_report(prediction_file, fusioncatcher_file, arriba_file,
                    rna_seq_cnv_log2foldchange_file, rna_seq_cnv_manual_an_table_file,
                    star_log_final_out_file, multiqc_fqc_right, multiqc_fqc_left,
                    comparison_file, output_file)

