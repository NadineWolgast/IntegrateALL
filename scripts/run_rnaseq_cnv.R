library(RNAseqCNV)
library(tidyverse)
library(readxl)
library(lubridate)

# Get input file path from command line arguments
args <- commandArgs(trailingOnly = TRUE)
config_file <- args[1]
metadata_file <- args[2]

RNAseqCNV_wrapper(config = config_file, metadata = metadata_file, snv_format = "custom", CNV_matrix = T, arm_lvl = T, batch = T, generate_weights = T)

# TODO: How to save in another directory?