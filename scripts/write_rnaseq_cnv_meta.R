args <- commandArgs(trailingOnly = TRUE)
sample_id <- args[1]
count_directory <- args[2]
vcf_directory <- args[3]

vcf_files <-  list.files(path = vcf_directory, pattern = "cancer.vcf")
count_files <-  list.files(path = count_directory, pattern = ".tsv")

out <- data.frame()
for(v in vcf_files){
    tmp <- c(sample_id, count_files, v)
    out <- rbind(out,tmp)
}



write.table(out, "meta.txt", row.names = FALSE, sep="\t", col.names = FALSE, quote = FALSE)

