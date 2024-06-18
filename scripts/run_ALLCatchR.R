library(ALLCatchRbcrabl1)

# Get input file path from command line arguments
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
out_file <- args[2]

#allcatch(input_file, ID_class = "ensemble_ID", sep = "\t")

allcatch_bcrabl1(input_file, ID_class = "ensemble_ID", sep = "\t", out.file = out_file)
