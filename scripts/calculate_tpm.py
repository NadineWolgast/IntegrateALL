import pandas as pd

# Load the ReadsPerGene.out.tab file
df = pd.read_csv(snakemake.input.reads_per_gene, sep='\t', skiprows=snakemake.params.skip_rows, names=["Gene", "Count"])

# Load the total mapped reads
with open(snakemake.input.total_mapped_reads, 'r') as total_reads_file:
    total_reads = int(total_reads_file.read().strip())

# Calculate TPM
df["TPM"] = (df["Count"] / df["Count"].sum()) * 1e6 / (total_reads / 1e6)

# Write TPM values to the output file
df.to_csv(snakemake.output.tpm, sep='\t', index=False)