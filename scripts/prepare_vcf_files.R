library(conflicted)
library(data.table)
library(readxl)
library(lubridate)
library(foreach)
library(magrittr)
library(dplyr)
library(readr)

args <- commandArgs(trailingOnly = TRUE)

vcf_directory<- args[1]
output_file <- args[2]
print(args)
sample_file_gz <- paste(vcf_directory, sep="")
sample_file <- R.utils::gunzip(sample_file_gz)

find_CHROM_row <- function(file_path) {
  con <- file(file_path, "r")
  line_num <- 0
  while (TRUE) {
    line <- readLines(con, n = 1)
    if (length(line) == 0) {
      break  # End of file
    }
    line_num <- line_num + 1
    if (grepl("#CHROM", line)) {
      return(line_num)
    }
  }
  return(NULL)  # #CHROM not found
}


# Find the row number where #CHROM occurs
skip <- find_CHROM_row(sample_file)
print("in prepare_vcf_files")
print(skip)


vcf <- read.csv(sample_file, skip=skip-1, header=TRUE, sep="\t")

#print(vcf)

AD <- as.character(sapply(vcf[,10], function(x) strsplit(x, ":", fixed = T)[[1]][2]))
DP <- as.character(sapply(vcf[,10], function(x) strsplit(x, ":", fixed = T)[[1]][3]))

Ref_AD <- as.numeric(sapply(AD, function(x) as.numeric(strsplit(x,",", fixed = T)[[1]][1])))
Alt_AD <- as.numeric(sapply(AD, function(x) as.numeric(strsplit(x,",", fixed = T)[[1]][2])))

vcf$DP_DIY <- DP
vcf$Ref_AD  <- Ref_AD
vcf$Alt_AD  <- Alt_AD
vcf$AF_DIY <- Alt_AD/(Alt_AD+Ref_AD)
colnames(vcf) <- gsub(".X", "",colnames(vcf))
colnames(vcf) <- gsub("X", "",colnames(vcf))
colnames(vcf) <- gsub(".CHROM", "#CHROM",colnames(vcf))
fwrite(vcf, file = paste(vcf_directory,"_clean.vcf", sep=""), sep = "\t",col.names=TRUE, row.names = F, quote = F)

data_file <- paste(vcf_directory,"_clean.vcf", sep="")
the_data <- read_tsv(data_file,
                     col_select = c("#CHROM", "POS", "DP_DIY", "AF_DIY"),
                     col_types  = c("f", "i", "n", "n")) %>% rename(chr = "#CHROM", start = "POS", depth = "DP_DIY", maf = "AF_DIY")

write.table(the_data, file =output_file, row.names = F,sep = "\t")