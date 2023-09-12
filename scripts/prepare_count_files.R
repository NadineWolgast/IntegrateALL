library(data.table)

args <- commandArgs(trailingOnly = TRUE)
sample_id <- args[1]
reads_directory <- args[2]

sample <- paste(reads_directory, "/", sample_id, "ReadsPerGene.out.tab", sep="")
print(sample)
reads_per_gene_file <- read.csv(sample, header=FALSE, sep= "\t", skip=4)
counts = reads_per_gene_file[,c(1,4)]
print(head(counts))
out_file_name <- paste(sample_id, ".txt", sep="")
write.table(counts, out_file_name, col.names=F, row.names=F)







