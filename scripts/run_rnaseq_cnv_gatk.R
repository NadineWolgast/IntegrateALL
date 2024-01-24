#library(conflicted)
library(RNAseqCNV)
#library(tidyverse)
#library(readxl)
#library(lubridate)
#library(foreach)
#library(doParallel)

# Get input file path from command line arguments
args <- commandArgs(trailingOnly = TRUE)
print(args)
config_file <- args[1]
metadata_file <- args[2]
sample_id <- args[3]
out_dir <- args[4]

new_meta <- paste0("data/",sample_id,"_meta.txt")
text <- paste(sample_id, "\t", paste0(sample_id, ".txt"), "\t", paste0(sample_id, "_Gatk.tsv"), sep = "")
writeLines(text, new_meta)

ou <- paste0('RNAseqCNV_output/gatk/', sample_id, "_gatk")

config_data <- list(
  out_dir = paste0("out_dir = '", ou, "'"),
  count_dir = "count_dir = 'data/single_counts'",
  snv_dir = "snv_dir = 'data/vcf_files/GATK'"
)

new_config <- paste0("data/",sample_id,"new_config.txt")
writeLines(unlist(config_data), new_config)

print(config_data)

RNAseqCNV_wrapper(config = new_config, metadata = new_meta, snv_format = "custom", CNV_matrix = T, arm_lvl = F, batch = F, generate_weights = F)
