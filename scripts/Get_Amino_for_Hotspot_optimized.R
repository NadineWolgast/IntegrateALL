#!/usr/bin/env Rscript
# Optimized R script for hotspot amino acid analysis with strand-specific handling.
# Processes pysamstats output to determine amino acid changes from sequence variations.

# Suppress warnings and messages for cleaner output
options(warn = -1)
suppressPackageStartupMessages({
  library(Biostrings)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript Get_Amino_for_Hotspot_optimized.R <pysamstats_directory> <gatk_vcf_file> <output_directory>")
}

pysamstats_directory <- args[1]
gatk_file <- args[2]
output_dir <- args[3]

# Optimized mutations data with strand information
mutations <- data.frame(
  Gene = c("ZEB2", "ZEB2", "KRAS", "KRAS", "KRAS", "KRAS", "NRAS", "NRAS", "NRAS", "NRAS", 
           rep("FLT3", 27), "PAX5", "IKZF1"),
  Hotspot = c("H1038", "Q1072", "G12", "G13", "Q16", "A146", "G12", "G13", "Q61", "A146", 
              "P857", "V843", "Y842", "N841", "D839", "M837", "I836", "D835", "R834",
              "A680", "N676", "A627", "K623", "Y599", "R595", "V592", "Y589", "N587",
              "G583", "Q580", "V579", "Q577", "L576", "E573", "Y572", "V491", "S446", 
              "P80R", "N159Y"),
  Chromosome = c("2", "2", "12", "12", "12", "12", "1", "1", "1", "1", 
                 rep("13", 27), "9", "7"),
  Strand = c("-", "-", "-", "-", "-", "-", "-", "-", "-", "-", 
             rep("-", 27), "-", "+"),  # Strand information from GTF analysis
  Start = c(144389982, 144389880, 25245349, 25245346, 25227341, 25225626, 
            114716125, 114716122, 114713907, 114709581, 
            28015672, 28018479, 28018482, 28018485, 28018491, 28018497, 28018500, 28018503, 28018506,
            28028191, 28028203, 28033948, 28033960, 28034122, 28034134, 28034143, 28034152, 28034158,
            28034170, 28034179, 28034182, 28034188, 28034191, 28034200, 28034203, 28035619, 28036015,
            37015167, 50382593),  # 1-based coordinates from original script
  End = c(144389984, 144389882, 25245351, 25245348, 25227343, 25225628, 
          114716127, 114716124, 114713909, 114709583,
          28015674, 28018481, 28018484, 28018487, 28018493, 28018499, 28018502, 28018505, 28018508,
          28028193, 28028205, 28033950, 28033962, 28034124, 28034136, 28034145, 28034154, 28034160,
          28034172, 28034181, 28034184, 28034190, 28034193, 28034202, 28034205, 28035621, 28036017,
          37015169, 50382596),
  stringsAsFactors = FALSE
)

# Function to find GATK VCF header row
find_CHROM_row_header <- function(file_path) {
  con <- file(file_path, "r")
  line_num <- 0
  while (TRUE) {
    line <- readLines(con, n = 1)
    if (length(line) == 0) break
    line_num <- line_num + 1
    if (grepl("#CHROM", line)) {
      close(con)
      return(line_num)
    }
  }
  close(con)
  return(NULL)
}

# Function to search GATK VCF for variants at specific position
find_GATK_variant <- function(file_path, chrom_val, pos_val) {
  skip <- find_CHROM_row_header(file_path)
  if (is.null(skip)) return(NULL)
  
  skip <- skip - 1
  df <- read.table(file_path, skip = skip, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Filter for variants in the region
  filtered_rows <- df[df$X.CHROM == chrom_val & df$POS >= pos_val & df$POS <= pos_val + 2, ]
  
  return(if (nrow(filtered_rows) > 0) filtered_rows else NULL)
}

# Function to determine the most likely alternative base
get_alternative_base <- function(data_row, ref_base) {
  # Column indices for nucleotides (A, C, T, G)
  nuc_cols <- c("A_pp", "C_pp", "T_pp", "G_pp")
  
  # Exclude reference base column
  alt_cols <- nuc_cols[nuc_cols != paste0(ref_base, "_pp")]
  
  # Find the nucleotide with highest count among alternatives
  if (length(alt_cols) > 0) {
    # Extract values for alternative nucleotides
    alt_values <- as.numeric(data_row[alt_cols])
    max_idx <- which.max(alt_values)
    max_col <- alt_cols[max_idx]
    return(substr(max_col, 1, 1))  # Extract nucleotide letter
  }
  
  return(ref_base)  # Fallback to reference
}

# Function to process a single hotspot file with strand-specific logic
process_hotspot_file <- function(file_path, mutations_df, gatk_file, output_dir) {
  file_name <- basename(file_path)
  
  # Skip if not TSV or CSV file
  if (!grepl("\\.(tsv|csv)$", file_name, ignore.case = TRUE)) return(NULL)
  
  # Read the data
  data <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  if (nrow(data) <= 1) return(NULL)  # Skip files with no data
  
  # Extract gene and hotspot from filename
  extracted_info <- sub(".*_([^_]+_[^_]+)\\.(tsv|csv)$", "\\1", file_name)
  gene <- sub("_.*", "", extracted_info)
  hotspot <- sub(".*_", "", extracted_info)
  
  # Get mutation information
  mutation_info <- mutations_df[mutations_df$Gene == gene & mutations_df$Hotspot == hotspot, ]
  if (nrow(mutation_info) == 0) {
    warning(paste("No mutation info found for", gene, hotspot))
    return(NULL)
  }
  
  chromosome <- mutation_info$Chromosome[1]
  strand <- mutation_info$Strand[1]
  start <- mutation_info$Start[1]
  end <- mutation_info$End[1]
  
  # Initialize new_base column
  data$new_base <- ""
  
  # Process each row to identify variants
  variant_found <- FALSE
  for (row in 1:nrow(data)) {
    # Check for significant mismatches (>5% mismatch rate)
    if (data$mismatches_pp[row] > 0 && 
        (data$mismatches_pp[row] / data$matches_pp[row]) >= 0.05) {
      
      ref_base <- data$ref[row]
      alt_base <- get_alternative_base(data[row, ], ref_base)
      data$new_base[row] <- alt_base
      variant_found <- TRUE
    } else {
      data$new_base[row] <- data$ref[row]
    }
  }
  
  if (!variant_found) return(NULL)  # Skip if no variants found
  
  # Collect bases and handle strand-specific processing
  bases_list <- character()
  percentage <- 0
  
  # Process bases based on strand orientation
  if (strand == "+") {
    # Plus strand: process bases in forward direction (IKZF1)
    for (k in 1:nrow(data)) {
      if (data$new_base[k] %in% c("A", "C", "T", "G") && data$reads_pp[k] >= 50) {
        bases_list <- c(bases_list, data$new_base[k])
        # Calculate percentage for the alternative base
        col_name <- paste0(data$new_base[k], "_pp")
        if (col_name %in% colnames(data)) {
          percentage <- data[[col_name]][k] * 100 / data$reads_pp[k]
        }
      } else {
        bases_list <- c(bases_list, data$ref[k])
      }
    }
    
    # Plus strand: use bases directly
    bases_string <- paste(bases_list, collapse = "")
    
  } else {
    # Minus strand: reverse complement processing (all other genes)
    # Reverse the data order first
    data_reversed <- data[nrow(data):1, ]
    
    for (k in 1:nrow(data_reversed)) {
      if (data_reversed$new_base[k] %in% c("A", "C", "T", "G") && data_reversed$reads_pp[k] >= 50) {
        bases_list <- c(bases_list, data_reversed$new_base[k])
        # Calculate percentage for the alternative base
        col_name <- paste0(data_reversed$new_base[k], "_pp")
        if (col_name %in% colnames(data_reversed)) {
          percentage <- data_reversed[[col_name]][k] * 100 / data_reversed$reads_pp[k]
        }
      } else {
        bases_list <- c(bases_list, data_reversed$ref[k])
      }
    }
    
    # Minus strand: reverse complement the sequence
    bases_string <- paste(bases_list, collapse = "")
    bases_string <- as.character(complement(DNAString(bases_string)))
  }
  
  # Translate to amino acid
  if (nchar(bases_string) %% 3 == 0 && nchar(bases_string) > 0) {
    tryCatch({
      new_amino_acid <- as.character(GENETIC_CODE[[bases_string]])
      if (is.na(new_amino_acid)) new_amino_acid <- "Unknown"
    }, error = function(e) {
      new_amino_acid <- "Invalid"
    })
  } else {
    new_amino_acid <- "Invalid_Length"
  }
  
  # Save results if significant variant found
  if (percentage > 0) {
    # Create output files
    output_file <- paste0(extracted_info, "_output_file.csv")
    
    # Save the detailed data
    selected_columns <- data[, c("chrom", "pos", "ref", "reads_pp", "mismatches_pp", 
                                "deletions_pp", "insertions_pp", "A_pp", "C_pp", 
                                "T_pp", "G_pp", "N_pp", "new_base")]
    
    # For minus strand, save in correct order
    if (strand == "-") {
      selected_columns <- selected_columns[nrow(selected_columns):1, ]
    }
    
    write.table(selected_columns, file.path(output_dir, paste0(extracted_info, "_with_bases.tsv")), 
                sep = "\t", row.names = FALSE, quote = FALSE)
    
    # Save summary result
    result <- data.frame(
      Gene = gene,
      Hotspot = hotspot, 
      Chromosome = chromosome,
      Start = start,
      End = end,
      Strand = strand,
      New_Aminoacid = new_amino_acid,
      Percentage = percentage,
      stringsAsFactors = FALSE
    )
    
    write.csv(result, file.path(output_dir, output_file), row.names = FALSE)
    
    # Search for GATK variant information
    gatk_result <- find_GATK_variant(gatk_file, chromosome, start)
    if (!is.null(gatk_result)) {
      write.table(gatk_result, file.path(output_dir, paste0(extracted_info, "_gatk_result.tsv")), 
                  sep = "\t", quote = FALSE, row.names = FALSE)
    } else {
      writeLines(paste("No GATK variant found for chromosome", chromosome, "at position", start),
                file.path(output_dir, paste0(extracted_info, "_gatk_result.tsv")))
    }
    
    cat("✅ Processed", gene, hotspot, paste0("(", strand, " strand)"), 
        "- New AA:", new_amino_acid, "- Percentage:", round(percentage, 2), "%\n")
    
    return(result)
  }
  
  return(NULL)
}

# Main processing
cat("📊 Starting optimized hotspot amino acid analysis\n")
cat("Processing directory:", pysamstats_directory, "\n")
cat("Output directory:", output_dir, "\n")

# Create output directory
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Get list of files to process
file_list <- list.files(pysamstats_directory, pattern = "\\.(tsv|csv)$", full.names = TRUE)
cat("Found", length(file_list), "files to process\n")

# Process files (can be parallelized for better performance)
results_list <- lapply(file_list, function(file) {
  process_hotspot_file(file, mutations, gatk_file, output_dir)
})

# Filter out NULL results
results_list <- results_list[!sapply(results_list, is.null)]

cat("✅ Successfully processed", length(results_list), "hotspot variants\n")

# Combine all results into summary file
if (length(results_list) > 0) {
  all_results <- do.call(rbind, results_list)
  write.csv(all_results, file.path(output_dir, "all_hotspots_summary.csv"), row.names = FALSE)
  cat("📝 Summary saved to all_hotspots_summary.csv\n")
}

cat("🎉 Hotspot analysis completed\n")