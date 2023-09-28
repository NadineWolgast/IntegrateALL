import pandas as pd


def create_sample_dataframe(sample_sheet):
    sample_df = pd.read_csv(sample_sheet)
    result_data = []
    for index, row in sample_df.iterrows():
        sample_identifier = row['sample_id']
        left_fastq = row['left']
        right_fastq = row['right']
        result_data.append({'sample_id': f"{sample_identifier}_left", 'FASTQ': left_fastq})
        result_data.append({'sample_id': f"{sample_identifier}_right", 'FASTQ': right_fastq})
    result_df = pd.DataFrame(result_data)
    return result_df
