# Performance optimizations for cluster execution
options(mc.cores = parallel::detectCores())
data.table::setDTthreads(0)  # Use all available threads

# Load required libraries
library(data.table)
library(dplyr)
library(parallel)

args <- commandArgs(trailingOnly = TRUE)

if(length(args) < 2) {
  stop("Usage: Rscript prepare_vcf-files_gatk.R <vcf_file> <output_file>")
}

vcf_file <- args[1]
output_file <- args[2]

cat("Processing VCF file:", vcf_file, "\n")
cat("Output file:", output_file, "\n")

# Check if input file exists
if(!file.exists(vcf_file)) {
  stop("VCF file does not exist: ", vcf_file)
}

# Read VCF data with optimized fread for better performance
cat("Reading VCF data...\n")
file <- fread(vcf_file, header = FALSE, sep = "\t", nThread = parallel::detectCores())

if(nrow(file) == 0) {
  cat("No AD entries found in VCF. Creating empty output file.\n")
  fwrite(data.table(chr = character(), start = integer(), depth = integer(), maf = numeric()), 
         output_file, sep = "\t")
  quit(save = "no", status = 0)
}

cat("VCF entries loaded:", nrow(file), "\n")

# Optimized AD position detection using vectorized operations
cat("Checking AD field positions...\n")
format_fields <- file$V9
ad_positions <- mclapply(format_fields, function(x) {
  fields <- strsplit(x, ":", fixed = TRUE)[[1]]
  which(fields == "AD")
}, mc.cores = parallel::detectCores())

ad_positions <- unlist(ad_positions)

# Verify AD is consistently at position 2
if (length(which(ad_positions == 2)) == nrow(file)) {
  cat("AD field found at position 2 in all entries. Processing...\n")
  
  # Optimized parallel processing of AD data
  sample_data <- file$V10
  
  # Extract AD values using parallel processing
  cat("Extracting allele depth information...\n")
  ad_values <- mclapply(sample_data, function(x) {
    strsplit(x, ":", fixed = TRUE)[[1]][2]
  }, mc.cores = parallel::detectCores())
  
  ad_values <- unlist(ad_values)
  
  # Process AD values to get REF and total depth
  cat("Calculating reference and total depths...\n")
  depth_data <- mclapply(ad_values, function(ad) {
    depths <- as.numeric(strsplit(ad, ",", fixed = TRUE)[[1]])
    list(ref = depths[1], total = sum(depths))
  }, mc.cores = parallel::detectCores())
  
  ref_depths <- sapply(depth_data, `[[`, "ref")
  total_depths <- sapply(depth_data, `[[`, "total")
  
  # Create optimized data.table output
  cat("Creating output table...\n")
  result_table <- data.table(
    chr = paste0("chr", file$V1),
    start = file$V2,
    depth = total_depths,
    maf = (total_depths - ref_depths) / total_depths
  )
  
  # Remove invalid entries (zero depth, infinite MAF)
  result_table <- result_table[is.finite(maf) & depth > 0]
  
  cat("Writing results to:", output_file, "\n")
  cat("Final entries:", nrow(result_table), "\n")
  
  # Write output using optimized fwrite
  fwrite(result_table, output_file, sep = "\t", nThread = parallel::detectCores())
  
  cat("VCF processing completed successfully!\n")
  
} else {
  cat("Warning: AD field not consistently at position 2\n")
  cat("AD positions found:", unique(ad_positions), "\n")
  # Create empty output for inconsistent format
  fwrite(data.table(chr = character(), start = integer(), depth = integer(), maf = numeric()), 
         output_file, sep = "\t")
}
