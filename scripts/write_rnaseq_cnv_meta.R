args <- commandArgs(trailingOnly = TRUE)
count_directory <- args[1]
vcf_directory <- args[2]
outputfile <- args[3]

# List files in the count_directory and vcf_directory
count_files <- list.files(path = count_directory, pattern = ".txt")
vcf_files <- list.files(path = vcf_directory, pattern = ".tsv")

# Extract file names without extensions
count_file_names <- sub("\\.txt$", "", count_files)
vcf_file_names <- sub("\\.tsv$", "", vcf_files)

# Create a data frame to store the results
result_df <- data.frame(
  Name = count_file_names,
  CountFile = count_files,
  VCFFile = vcf_files
)

# Write the results to the output file
write.table(result_df, file = outputfile, sep = "\t", row.names = FALSE, col.names= FALSE, quote = FALSE)




