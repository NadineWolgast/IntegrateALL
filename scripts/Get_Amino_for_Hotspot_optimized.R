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

# Load hotspots configuration from CSV file
script_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))
if (length(script_dir) == 0) {
  # If script_dir is empty, use current directory (for interactive sessions)
  script_dir <- getwd()
}
config_file <- file.path(dirname(script_dir), "config", "hotspots_config.csv")

if (!file.exists(config_file)) {
  stop(paste("Hotspots config file not found:", config_file))
}

cat("Loading hotspots configuration from", config_file, "\n")
mutations <- read.csv(config_file, stringsAsFactors = FALSE)

# Rename columns to match expected format in the rest of the script
colnames(mutations) <- c("Gene", "Hotspot", "Chromosome", "Start", "End", "Strand")

# Ensure correct data types
mutations$Chromosome <- as.character(mutations$Chromosome)
mutations$Start <- as.integer(mutations$Start)
mutations$End <- as.integer(mutations$End)
mutations$Strand <- as.character(mutations$Strand)

cat("Loaded", nrow(mutations), "hotspot configurations\n")

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

# Function to search GATK VCF for variants in a genomic region
# Note: VCF uses 1-based coordinates, BAM/BED uses 0-based
find_GATK_variant <- function(file_path, chrom_val, start_pos, end_pos) {
  skip <- find_CHROM_row_header(file_path)
  if (is.null(skip)) {
    cat("    ⚠️  Could not find VCF header in file\n")
    return(NULL)
  }

  skip <- skip - 1

  # Try to read VCF with better error handling
  df <- tryCatch({
    read.table(file_path, skip = skip, header = TRUE, sep = "\t",
               stringsAsFactors = FALSE, comment.char = "", quote = "")
  }, error = function(e) {
    cat("    ⚠️  Error reading VCF file:", e$message, "\n")
    return(NULL)
  })

  if (is.null(df) || nrow(df) == 0) {
    cat("    ⚠️  VCF file is empty or could not be read\n")
    return(NULL)
  }

  # Debug: show column names
  cat("    VCF columns:", paste(head(colnames(df), 8), collapse=", "), "...\n")
  cat("    Total variants in VCF:", nrow(df), "\n")

  # Filter for variants in the region
  # Add ±1 buffer to account for potential 0-based vs 1-based coordinate differences
  # VCF is 1-based, pysamstats output might be 0-based
  search_start <- start_pos - 1
  search_end <- end_pos + 1

  cat("    Searching GATK VCF: chr", chrom_val, ":", search_start, "-", search_end, "\n", sep="")

  # Try both with and without "chr" prefix for chromosome matching
  chrom_variants <- c(as.character(chrom_val), paste0("chr", chrom_val))

  # Handle different possible column names for chromosome
  chrom_col <- NULL
  if ("X.CHROM" %in% colnames(df)) {
    chrom_col <- "X.CHROM"
  } else if ("CHROM" %in% colnames(df)) {
    chrom_col <- "CHROM"
  } else if ("#CHROM" %in% colnames(df)) {
    chrom_col <- "#CHROM"
  } else {
    # Take first column as chromosome
    chrom_col <- colnames(df)[1]
    cat("    ⚠️  Using column", chrom_col, "as chromosome\n")
  }

  # Filter for variants
  filtered_rows <- df[df[[chrom_col]] %in% chrom_variants &
                      df$POS >= search_start &
                      df$POS <= search_end, ]

  if (nrow(filtered_rows) > 0) {
    cat("    ✓ Found", nrow(filtered_rows), "GATK variant(s) in region\n")
    # Show the variants found
    for (i in 1:nrow(filtered_rows)) {
      cat("      - Pos:", filtered_rows$POS[i],
          "Ref:", filtered_rows$REF[i],
          "Alt:", filtered_rows$ALT[i], "\n")
    }
  } else {
    cat("    ✗ No GATK variants found in region\n")
    # Debug: show what chromosomes are in the VCF
    unique_chroms <- unique(df[[chrom_col]])
    if (length(unique_chroms) > 0) {
      cat("    VCF chromosomes:", paste(head(unique_chroms, 5), collapse=", "),
          ifelse(length(unique_chroms) > 5, "...", ""), "\n")
    }
    # Show nearby variants for debugging
    nearby <- df[df[[chrom_col]] %in% chrom_variants &
                 df$POS >= (start_pos - 100) &
                 df$POS <= (end_pos + 100), ]
    if (nrow(nearby) > 0) {
      cat("    Nearby variants (±100bp):\n")
      for (i in 1:min(3, nrow(nearby))) {
        cat("      - Pos:", nearby$POS[i], "Ref:", nearby$REF[i], "Alt:", nearby$ALT[i], "\n")
      }
    }
  }

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

  # Read the data with error handling
  data <- tryCatch({
    read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  }, error = function(e) {
    warning(paste("Could not read file", file_name, ":", e$message))
    return(NULL)
  })

  if (is.null(data) || nrow(data) <= 1) return(NULL)  # Skip files with no data

  # Extract gene and hotspot from filename
  extracted_info <- sub(".*_([^_]+_[^_]+)\\.(tsv|csv)$", "\\1", file_name)
  gene <- sub("_.*", "", extracted_info)
  hotspot <- sub(".*_", "", extracted_info)

  cat("Processing:", gene, hotspot, "...\n")
  
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
  variant_position <- NA
  alt_base_found <- ""

  for (row in 1:nrow(data)) {
    # Check for significant mismatches (>5% mismatch rate)
    if (data$mismatches_pp[row] > 0 &&
        (data$mismatches_pp[row] / data$matches_pp[row]) >= 0.05) {

      ref_base <- data$ref[row]
      alt_base <- get_alternative_base(data[row, ], ref_base)
      data$new_base[row] <- alt_base
      variant_found <- TRUE
      variant_position <- data$pos[row]
      alt_base_found <- paste0(ref_base, ">", alt_base)
    } else {
      data$new_base[row] <- data$ref[row]
    }
  }

  if (!variant_found) return(NULL)  # Skip if no variants found

  # CRITICAL: Use the exact codon coordinates from the config file
  # The config contains start and end positions (1-based) that define the codon
  # Extract only rows that match these exact positions

  # For 3-base codon: positions should be start, start+1, start+2 (1-based)
  # But config end is exclusive, so codon positions are: start to end-1
  codon_positions <- seq(start, end - 1)

  # Filter data to only include rows with positions in the codon
  codon_data <- data[data$pos %in% codon_positions, ]

  # Verify we have exactly 3 bases
  if (nrow(codon_data) != 3) {
    cat("⚠️  Warning: Expected 3 bases in codon, found", nrow(codon_data), "\n")
    cat("    Config start:", start, "end:", end, "\n")
    cat("    Expected positions:", paste(codon_positions, collapse=", "), "\n")
    cat("    Found positions:", paste(codon_data$pos, collapse=", "), "\n")
    # Continue anyway, but this might produce wrong results
  }

  cat("    Codon window (1-based positions):", paste(codon_data$pos, collapse=", "), "\n")

  # Collect bases and handle strand-specific processing
  bases_list <- character()
  percentage <- 0

  # Process bases based on strand orientation
  if (strand == "+") {
    # Plus strand: process bases in forward direction (IKZF1)
    for (k in 1:nrow(codon_data)) {
      if (codon_data$new_base[k] %in% c("A", "C", "T", "G") && codon_data$reads_pp[k] >= 50) {
        bases_list <- c(bases_list, codon_data$new_base[k])
        # Calculate VAF ONLY for variant positions (ref != new_base)
        if (codon_data$ref[k] != codon_data$new_base[k]) {
          ref_col <- paste0(codon_data$ref[k], "_pp")
          alt_col <- paste0(codon_data$new_base[k], "_pp")
          if (ref_col %in% colnames(codon_data) && alt_col %in% colnames(codon_data)) {
            ref_count <- codon_data[[ref_col]][k]
            alt_count <- codon_data[[alt_col]][k]
            # VAF = alternative / (reference + alternative) * 100
            percentage <- alt_count * 100 / (ref_count + alt_count)
          }
        }
      } else {
        bases_list <- c(bases_list, codon_data$ref[k])
      }
    }

    # Plus strand: use bases directly
    bases_string <- paste(bases_list, collapse = "")

    # Convert U (RNA) to T (DNA) if present - handle mixed RNA/DNA sequences
    bases_string <- gsub("U", "T", bases_string, ignore.case = TRUE)

  } else {
    # Minus strand: reverse complement processing (all other genes)
    # Reverse the codon data order
    codon_data_reversed <- codon_data[nrow(codon_data):1, ]

    for (k in 1:nrow(codon_data_reversed)) {
      if (codon_data_reversed$new_base[k] %in% c("A", "C", "T", "G") && codon_data_reversed$reads_pp[k] >= 50) {
        bases_list <- c(bases_list, codon_data_reversed$new_base[k])
        # Calculate VAF ONLY for variant positions (ref != new_base)
        if (codon_data_reversed$ref[k] != codon_data_reversed$new_base[k]) {
          ref_col <- paste0(codon_data_reversed$ref[k], "_pp")
          alt_col <- paste0(codon_data_reversed$new_base[k], "_pp")
          if (ref_col %in% colnames(codon_data_reversed) && alt_col %in% colnames(codon_data_reversed)) {
            ref_count <- codon_data_reversed[[ref_col]][k]
            alt_count <- codon_data_reversed[[alt_col]][k]
            # VAF = alternative / (reference + alternative) * 100
            percentage <- alt_count * 100 / (ref_count + alt_count)
          }
        }
      } else {
        bases_list <- c(bases_list, codon_data_reversed$ref[k])
      }
    }

    # Minus strand: reverse complement the sequence
    bases_string <- paste(bases_list, collapse = "")

    # Convert U (RNA) to T (DNA) if present - handle mixed RNA/DNA sequences
    bases_string <- gsub("U", "T", bases_string, ignore.case = TRUE)

    # Error handling for DNAString conversion
    tryCatch({
      bases_string <- as.character(complement(DNAString(bases_string)))
    }, error = function(e) {
      cat("⚠️  Error in complement for", gene, hotspot, ":", e$message, "\n")
      cat("    Sequence:", bases_string, "\n")
      # Try to clean the sequence - keep only valid DNA bases
      bases_string <<- gsub("[^ACGTacgt]", "N", bases_string)
      cat("    Cleaned sequence:", bases_string, "\n")
      bases_string <<- as.character(complement(DNAString(bases_string)))
    })
  }
  
  # Translate to amino acid with improved error handling
  new_amino_acid <- tryCatch({
    seq_len <- nchar(bases_string)

    if (seq_len == 3) {
      # Perfect codon length - translate directly
      aa <- as.character(GENETIC_CODE[[bases_string]])
      if (is.na(aa)) {
        cat("⚠️  Unknown codon:", bases_string, "\n")
        "Unknown"
      } else {
        aa
      }
    } else {
      # Unexpected length - this should not happen with the new codon extraction logic
      cat("⚠️  Unexpected sequence length", seq_len, "bp for", gene, hotspot, "\n")
      cat("    Sequence:", bases_string, "\n")

      if (seq_len > 3) {
        # Try to extract 3 bases from the middle/start
        codon <- substr(bases_string, 1, 3)
        aa <- as.character(GENETIC_CODE[[codon]])
        if (!is.na(aa)) {
          cat("    Using first 3 bases:", codon, "->", aa, "\n")
          bases_string <<- codon  # Update for output
          aa
        } else {
          paste0("Invalid_Length_", seq_len)
        }
      } else {
        paste0("TooShort_", seq_len, "bp")
      }
    }
  }, error = function(e) {
    cat("⚠️  Translation error for", gene, hotspot, ":", e$message, "\n")
    cat("    Sequence:", bases_string, "(length:", nchar(bases_string), ")\n")
    "Translation_Error"
  })
  
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
    
    # Save summary result with enhanced information
    result <- data.frame(
      Gene = gene,
      Hotspot = hotspot,
      Chromosome = chromosome,
      Start = start,
      End = end,
      Strand = strand,
      Variant = alt_base_found,
      Variant_Position = variant_position,
      Sequence = bases_string,
      Sequence_Length = nchar(bases_string),
      New_Aminoacid = new_amino_acid,
      Percentage = percentage,
      stringsAsFactors = FALSE
    )
    
    write.csv(result, file.path(output_dir, output_file), row.names = FALSE)
    
    # Search for GATK variant information in the entire hotspot region
    cat("  Searching GATK VCF for", gene, hotspot, "in region chr", chromosome, ":", start, "-", end, "\n", sep="")
    gatk_result <- find_GATK_variant(gatk_file, chromosome, start, end)
    if (!is.null(gatk_result)) {
      write.table(gatk_result, file.path(output_dir, paste0(extracted_info, "_gatk_result.tsv")),
                  sep = "\t", quote = FALSE, row.names = FALSE)
    } else {
      writeLines(paste("No GATK variant found for chromosome", chromosome, "in region", start, "-", end),
                file.path(output_dir, paste0(extracted_info, "_gatk_result.tsv")))
    }
    
    # Enhanced output with sequence and variant information
    cat("✅ Processed", gene, hotspot, paste0("(", strand, " strand)"), "\n")
    cat("    Variant:", alt_base_found, "at position", variant_position, "\n")
    cat("    Sequence:", bases_string, "(length:", nchar(bases_string), "bp)\n")
    cat("    New AA:", new_amino_acid, "- Percentage:", round(percentage, 2), "%\n")
    
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

# Process files with error handling (can be parallelized for better performance)
results_list <- lapply(file_list, function(file) {
  tryCatch({
    process_hotspot_file(file, mutations, gatk_file, output_dir)
  }, error = function(e) {
    cat("❌ Error processing", basename(file), ":", e$message, "\n")
    return(NULL)
  })
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