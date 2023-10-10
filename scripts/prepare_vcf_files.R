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
sample_file <- paste(vcf_directory, sep="")

vcf <- read.csv(sample_file, skip=263, header=TRUE, sep="\t")

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