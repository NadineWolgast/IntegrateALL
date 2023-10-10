import csv


def extract_filenames_from_csv(csv_file):
    left_filenames = []
    right_filenames = []

    with open(csv_file, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            left_filenames.append(row['left'].split('/')[-1])  # Extract only the filename
            right_filenames.append(row['right'].split('/')[-1])  # Extract only the filename

    return left_filenames, right_filenames
