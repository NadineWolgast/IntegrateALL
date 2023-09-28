import os
import pandas as pd


def merge_reads_per_gene_files(input_dir, out_file):
    # Initialize an empty DataFrame to store the merged data
    merged_df = pd.DataFrame()

    # Iterate through subdirectories (each corresponds to a sample)
    for sample_dir in os.listdir(input_dir):
        sample_path = os.path.join(input_dir, sample_dir)

        # Check if the path is a directory
        if os.path.isdir(sample_path):
            sample_name = os.path.basename(sample_path)
            readsPerGene_file = os.path.join(sample_path, 'ReadsPerGene.out.tab')

            # Check if the ReadsPerGene.out.tab file exists
            if os.path.exists(readsPerGene_file):
                # Read the file, skip the first 4 rows, and select the first and fourth columns
                data = pd.read_csv(readsPerGene_file, sep='\t', header=None, skiprows=4, usecols=[0, 1],
                                   names=['Gene', sample_name])

                # Set the gene name column as the index
                data.set_index('Gene', inplace=True)

                # If merged_df is empty, initialize it with the first data
                if merged_df.empty:
                    merged_df = data
                else:
                    # Merge data with the existing merged_df based on the index (gene names)
                    merged_df = pd.merge(merged_df, data, left_index=True, right_index=True, how='outer')

    # Reset the index to have the gene names as a separate column
    merged_df.reset_index(inplace=True)

    # Save the merged DataFrame to the output file
    merged_df.to_csv(out_file, sep='\t', index=False)
