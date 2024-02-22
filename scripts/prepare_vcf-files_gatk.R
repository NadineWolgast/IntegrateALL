#library(conflicted)
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

file <- read.csv(paste0(vcf_directory,'_sel'),
                 header = F,
                 sep = "\t")
head(file)
table(file$V1)
head(file)
table(file$V9)

# check whether AD is always at pos2
AD_pos <- unlist(lapply(sapply(file$V9, function(x) strsplit(x , ":")), function(x) which(x == "AD")))

if (length(which(AD_pos == 2)) == nrow(file)) {
  
  AD <- as.character(sapply(file$V10, function(x) strsplit(x, ":", fixed = T)[[1]][[2]]))
  REF <- as.numeric(as.character(sapply(AD, function(x) strsplit(x, ",", fixed = T)[[1]][[1]])))
  DP <- sapply(AD, function(x) strsplit(x, ",", fixed = T)[[1]])
  DP <- unlist(lapply(DP, function(x) sum(as.numeric(x))))
  
  tab <- data.frame("chr" = paste0("chr", file$V1),
                    "start" = file$V2,
                    "depth" = DP,
                    "maf" = as.numeric((DP-REF)/(DP)),
                    stringsAsFactors = F)
  write.table(tab, output_file, sep = "\t",row.names = F)
  
}
