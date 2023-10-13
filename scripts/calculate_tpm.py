import pandas as pd
import sys


reads_per_gene_file = sys.argv[1]
total_mapped_reads_file = sys.argv[2]
outfile = sys.argv[3]

# Load the ReadsPerGene.out.tab file
df = pd.read_csv(reads_per_gene_file, sep='\t', skiprows=4, names=["Gene", "Count"], index_col=False)

# Load the total mapped reads
with open(total_mapped_reads_file, 'r') as total_reads_file:
    total_reads = int(total_reads_file.read().strip())

# Calculate TPM
df["TPM"] = (df["Count"] / df["Count"].sum()) * 1e6 / (total_reads / 1e6)

# Write TPM values to the output file
df.to_csv(outfile, sep='\t', index=False)
