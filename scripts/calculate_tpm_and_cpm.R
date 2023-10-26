library(data.table)

args <- commandArgs(trailingOnly = TRUE)
#print(args)
reads_per_gene_file <- args[1]
outfile <- args[2]
outfile_cpm <- args[3]

counts <- fread(reads_per_gene_file, skip = 4, header = FALSE, select = 1:2, col.names = c("Gene", "Count"))

# Lade die Gesamtanzahl der zugeordneten Reads
total_reads <- fread("annotation/cds_length.tsv", header = TRUE)
total_reads <- as.data.frame(total_reads)
which(is.na(counts$Gene))
ma <- match(counts$Gene, total_reads$ESNG)

total_reads <- total_reads[match(counts$Gene, total_reads$ESNG),]
total_reads$ESNG <- counts$Gene
total_reads$cds_length[which(is.na(total_reads$cds_length))] <- median(total_reads$cds_length, na.rm = T)

table(counts$Gene == total_reads$ESNG)
#total_reads <- na.omit(total_reads)

tpm3 <- function(counts,total_reads) {

  x <- counts$Count/total_reads$cds_length
  return(x*1e6/sum(x, na.rm = T))
}

tpm <- tpm3(counts,total_reads)

# Schreibe TPM-Werte in die Ausgabedatei
write.table(cbind(counts$Gene,tpm), file = outfile, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

countmatrix <- data.frame(counts$Count)
rownames(countmatrix) <- counts$Gene
cpm <- apply(countmatrix,2, function(x) (x/sum(x))*1000000)
colnames(cpm)<- c("CPM")
write.table(cbind(rownames(cpm),cpm), file = outfile_cpm, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
