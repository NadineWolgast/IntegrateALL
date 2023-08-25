library(ALLCatchR)

# Get input file path from command line arguments
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]

allcatch(input_file, ID_class = "symbol", sep = "\t")
# TODO: How to save in another directory?
