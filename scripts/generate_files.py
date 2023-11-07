import os
import csv


def generate_files(input_csv, output_dir):
    # Create the output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Initialize list for meta.txt data
    meta_data = []
    config_data = [
        "out_dir = 'RNAseqCNV_output'",
        "count_dir = 'data/single_counts'",
        "snv_dir = 'data/vcf_files'"
    ]

    # Read the sample.csv file
    with open(input_csv, 'r') as csvfile:
        csvreader = csv.DictReader(csvfile)
        for row in csvreader:
            sample_id = row['sample_id']
            # Add data to meta.txt
            meta_data.append(f"{sample_id}\t{sample_id}.txt\t{sample_id}.tsv")

    # Write data to meta.txt
    meta_file = os.path.join(output_dir, 'meta.txt')
    with open(meta_file, 'w') as meta:
        meta.write("\n".join(meta_data))

    # Write data to config.txt
    config_file = os.path.join(output_dir, 'config.txt')
    with open(config_file, 'w') as config:
        config.write("\n".join(config_data))

    return
