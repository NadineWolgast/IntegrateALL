args <- commandArgs(trailingOnly = TRUE)
out_directory <- args[1]
count_directory <- args[2]
snv_directory <- args[3]


out = paste0('out_dir = "', out_directory, '"')
count = paste0('count_dir = "', count_directory, '"')
snv = paste0('snv_dir = "', snv_directory, '"')
cat(out, count, snv, sep ="\n", file="config.txt")


