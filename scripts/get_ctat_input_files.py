import pandas as pd
import os


def get_ctat_input_files(samples_csv):
    df = pd.read_csv(samples_csv, sep=',')

    # Initialize an empty dictionary to store the results
    samples_test = {}

    # Iterate over the rows of the DataFrame
    for index, row in df.iterrows():
        sample_id = row['sample_id']
        left = os.path.basename(row['left'])
        right = os.path.basename(row['right'])

        # Add a backslash in front of the filenames
        left = '/' + left
        right = '/' + right

        # Create a nested dictionary for each sample
        sample_info = {
            "left": left,
            "right": right
        }

        # Add the sample to the samples_test dictionary
        samples_test[sample_id] = sample_info

    # Print the resulting dictionary
    print("samples_test:",samples_test)
    return samples_test
