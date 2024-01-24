import pysam
import pysamstats
import sys
from pysam import Fastafile
import pandas as pd
import csv

"""mutations = pd.DataFrame({
    'Gene': ["ZEB2", "ZEB2", "KRAS", "KRAS", "KRAS", "KRAS", "NRAS", "NRAS", "NRAS", "NRAS"],
    'Hotspot': ["H1038", "Q1072", "G12", "G13", "Q16", "A146", "G12", "G13", "Q61", "A146"],
    'Chromosome': ["2", "2", "12", "12", "12", "12", "1", "1", "1", "1"],
    'Start': [144389982, 144389880, 25245349, 25245346, 25227341, 25225626, 114716125, 114716122, 114713907, 114709581],
    'End': [144389985, 144389883, 25245352, 25245349, 25227344, 25225629, 114716128, 114716125, 114713910, 114709584]
})"""


mutations = pd.DataFrame({
    'Gene': ["ZEB2", "ZEB2", "KRAS", "KRAS", "KRAS", "KRAS", "NRAS", "NRAS", "NRAS", "NRAS", "FLT3", "FLT3", "FLT3", "FLT3", "FLT3", "FLT3", "FLT3", "FLT3", "FLT3",
             "FLT3", "FLT3", "FLT3", "FLT3", "FLT3", "FLT3", "FLT3", "FLT3", "FLT3",
             "FLT3", "FLT3", "FLT3", "FLT3", "FLT3", "FLT3", "FLT3", "FLT3", "FLT3", "PAX5", "IKZF1"],
    'Hotspot': ["H1038", "Q1072", "G12", "G13", "Q16", "A146", "G12", "G13", "Q61", "A146", "P857", "V843", "Y842", "N841", "D839", "M837", "I836", "D835", "R834",
                "A680", "N676", "A627", "K623", "Y599", "R595", "V592", "Y589", "N587",
                "G583", "Q580", "V579", "Q577", "L576", "E573", "Y572", "V491", "S446", "P80R", "N159Y"],
    'Chromosome': ["2", "2", "12", "12", "12", "12", "1", "1", "1", "1", "13", "13", "13", "13", "13", "13", "13", "13", "13",
                   "13", "13", "13", "13", "13", "13", "13", "13", "13",
                   "13", "13", "13", "13", "13", "13", "13", "13", "13", "9", "7"],
    'Start': [144389981, 144389879, 25245348, 25245345, 25227340, 25225625, 114716124, 114716121, 114713906, 114709580,
              28015671, 28018478, 28018481, 28018484, 28018490, 28018496, 28018499, 28018502, 28018505,
              28028190, 28028202, 28033947, 28033959, 28034121, 28034133, 28034142, 28034151, 28034157,
              28034169, 28034178, 28034181, 28034187, 28034190, 28034199, 28034202, 28035618, 28036014, 37015164, 50382592],
    'End': [144389984, 144389882, 25245351, 25245348, 25227343, 25225628, 114716127, 114716124, 114713909, 114709583,
            28015674, 28018481, 28018484, 28018487, 28018493, 28018499, 28018502, 28018505, 28018508,
            28028193, 28028205, 28033950, 28033962, 28034124, 28034136, 28034145, 28034154, 28034160,
            28034172, 28034181, 28034184, 28034190, 28034193, 28034202, 28034205, 28035621, 28036017, 37015167, 50382595]
})


def calculate_pysamstats(bam_file, fasta_file, output_file, start, stop, chromosome):
    mybam = pysam.AlignmentFile(bam_file)
    myfasta = Fastafile(fasta_file)
    b = pysamstats.load_variation(mybam, myfasta, chrom=chromosome, start=start, end=stop, truncate=True)

    print(b)
    with open(output_file, 'w') as file:
        # Write header
        header = "chrom\tpos\tref\treads_all\treads_pp\tmatches\tmatches_pp\tmismatches\tmismatches_pp\tdeletions\tdeletions_pp\tinsertions\tinsertions_pp\tA\tA_pp\tC\tC_pp\tT\tT_pp\tG\tG_pp\tN\tN_pp\n"
        file.write(header)

        for entry in b:
            # Convert bytes to strings and remove 'b'
            entry = [str(val).replace("b'", "").replace("'", "") for val in entry]
            file.write('\t'.join(entry) + '\n')


if __name__ == "__main__":
    print(sys.argv)

    bam_file, fasta_file, sample_id, output_file = sys.argv[1:]

    for index, row in mutations.iterrows():
        hotspot = row['Hotspot']
        gene = row['Gene']
        chromosome = row['Chromosome']
        start = row['Start']
        end = row['End']

        output_file = f"pysamstats_output_dir/{sample_id}/{sample_id}_{gene}_{hotspot}.tsv"
        calculate_pysamstats(bam_file, fasta_file, output_file, start, end, chromosome)

    #calculate_pysamstats(bam_file, fasta_file, sample_id, output_file, start, stop, chromosome, Hotspot)
