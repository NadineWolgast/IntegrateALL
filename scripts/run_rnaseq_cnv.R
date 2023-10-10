library(conflicted)
library(RNAseqCNV)
library(tidyverse)
library(readxl)
library(lubridate)
library(foreach)
library(doParallel)

# Get input file path from command line arguments
args <- commandArgs(trailingOnly = TRUE)
print(args)
config_file <- args[1]
metadata_file <- args[2]

RNAseqCNV_wrapper(config = config_file, metadata = metadata_file, snv_format = "custom", CNV_matrix = T, arm_lvl = F, batch = F, generate_weights = F)

# TODO: How to save in another directory?
# TODO: Test with batch = T for more than 20 samples and set generate_weights = T