import os
import csv


def validate_input(samples_csv):

    # Initialize a list to store error messages
    errors = []

    # Check if the samples.csv file exists
    if not os.path.exists(samples_csv):
        errors.append(f"{samples_csv} does not exist")

    # Check if the samples.csv file has the correct format
    try:
        with open(samples_csv, "r") as csv_file:
            csv_reader = csv.reader(csv_file)
            header = next(csv_reader)  # Read the header row

            if len(header) != 3 or header[0] != "sample_id" or header[1] != "left" or header[2] != "right":
                errors.append(f"Invalid header in {samples_csv}. Expected: 'sample_id, left, right'")

            for row in csv_reader:
                sample, left, right = row
                if not os.path.exists(left):
                    errors.append(f"File referenced in 'left' column does not exist: {left}")
                if not os.path.exists(right):
                    errors.append(f"File referenced in 'right' column does not exist: {right}")
    except Exception as e:
        errors.append(f"Error reading {samples_csv}: {str(e)}")
    return errors
