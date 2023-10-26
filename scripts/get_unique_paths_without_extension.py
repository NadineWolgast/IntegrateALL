import os


def get_unique_paths_without_extension(dataframe):
    unique_paths = set()

    for fastq_path in dataframe['FASTQ']:
        # Extract the directory portion of the path and remove the file extension
        path_without_extension = os.path.splitext(os.path.dirname(fastq_path))[0]
        unique_paths.add(path_without_extension)
    print("get unique",list(unique_paths))

    return list(unique_paths)
