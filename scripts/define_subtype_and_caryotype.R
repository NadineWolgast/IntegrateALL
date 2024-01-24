library("readr")
library("utils")

args <- commandArgs(trailingOnly = TRUE)

#print(args)

prediction_file <- args[1]
rna_seq_cnv_estimation_file <- args[2]
fusioncatcher_file <- args[3]
arriba_file <- args[4]
chromosome_counts_karyotype_file <- args[5]
anno_gene_fusions_file <- args[6]
outfile <- args[7]
#print(arriba_file)

subtype_data <- read.table(anno_gene_fusions_file, header = TRUE, sep = "")

prediction_data <- read.table(prediction_file, header=TRUE, sep="\t")

# Extract sample-specific information
sample_info <- prediction_data$sample
subtype <- prediction_data$Prediction
confidence <- prediction_data$Confidence
karyotype_data <- read.table(rna_seq_cnv_estimation_file, header=TRUE, sep="\t")
chromosome_number <- karyotype_data$chrom_n


# Define the dataframe
chromosome_counts_to_karyotype <- data.frame(
  Chromosome_count = c("24-30", "31-39", "51-58", "59-78"),
  Karyotype = c("near haploid", "low hypodiploid", "high hyperdiploid", "near triploid"),
  Subtype = c("Near haploid", "Low hypodiploid", "Hyperdiploid", "Low hypodiploid")
)

# Define the function to check Subtype and Chromosome_number
check_subtype_and_chromosome <- function(Subtype, Chromosome_number,arriba_file, fusioncatcher_file) {

  # Check if Subtype is in the dataframe
  if (Subtype %in% chromosome_counts_to_karyotype$Subtype) {

    # Find the corresponding row for the Subtype
    matching_row <- chromosome_counts_to_karyotype[chromosome_counts_to_karyotype$Subtype == Subtype, ]

    # Extract the Chromosome_count range for the Subtype
    chromosome_range <- matching_row$Chromosome_count

    # Check if Chromosome_number falls within the range
    if (Chromosome_number %in% unlist(strsplit(chromosome_range, "-"))) {
      return(paste("Bei", Subtype, ": Analysis of the virtual karyotype established a chromosome count of",
                    Chromosome_number, "confirming the subtype allocation."))
    } else {
      return(paste("Bei", Subtype, ": Analysis of the virtual karyotype established a chromosome count of",
                    Chromosome_number, "not confirming the subtype allocation."))
    }
  } else {
        #print(subtype_data$Subtype)
         if (Subtype %in% subtype_data$Subtype) {
            matching_row_index <- which(subtype_data$Subtype == Subtype)
             if (length(matching_row_index) > 0) {
                # Extract the matching row from the dataframe
                matching_row <- subtype_data[matching_row_index, ]
                print("matching_row")
                print(matching_row)

                # Access specific columns
                X5_end_partner <- matching_row$X5_end_partner
                print("X5_end_partner")
                print(X5_end_partner)
                X3_end_partner <- matching_row$X3_end_partner
                print("X3_end_partner")
                print(X3_end_partner)

            } else {
                cat("No matching row found for Subtype.\n")
            }

            # Check if X5_end_partner and X3_end_partner are present in both arriba and fusioncatcher files

            arriba_data <- read.csv(arriba_file, header = TRUE, sep = "\t")


            fusioncatcher_data <- read.table(fusioncatcher_file, header = TRUE, sep = "\t")

            fusioncatcher_found <- list()
            arriba_found <- list()
            k <- 1
            l <- 1
            # Loop through each row of pairs of X5_end_partner and X3_end_partner
            for (i in 1:nrow(matching_row)) {
                row <- matching_row[i, ]

                X5_end_partner <- row$X5_end_partner
                X3_end_partner <- row$X3_end_partner

                # Check if the pair is present in FusionCatcher data
                matching_fc_rows <- which(
                    fusioncatcher_data$Gene_1_symbol.5end_fusion_partner. %in% c(X5_end_partner) &
                    fusioncatcher_data$Gene_2_symbol.3end_fusion_partner. %in% c(X3_end_partner)
                )

                if (length(matching_fc_rows) > 0) {
                    print("matching_fc_rows")
                    print(matching_fc_rows)
                    found_fusions <- fusioncatcher_data[matching_fc_rows, c("Gene_1_symbol.5end_fusion_partner.", "Gene_2_symbol.3end_fusion_partner.")]
                    fusioncatcher_found[[k]] <- found_fusions
                    k <- k + 1
                } else {
                    invisible()
                }


                # Check if the pair is present in Arriba data
                matching_arriba_rows <- which(
                    arriba_data$X.gene1 %in% c(X5_end_partner) &
                    arriba_data$gene2 %in% c(X3_end_partner)
                )

                if (length(matching_arriba_rows) > 0) {
                    arriba_found[[l]] <- arriba_data[matching_arriba_rows, c("X.gene1", "gene2")]
                    l <- l + 1
                } else {
                    invisible()
                }
            }


            # Check if any fusions were found in FusionCatcher or Arriba data
            if (length(fusioncatcher_found) > 0 || length(arriba_found) > 0) {
                result_text <- "Gene fusion calling identified"

                if (length(fusioncatcher_found) > 0) {
                    result_text <- paste(result_text, "fusions confirmed by FusionCatcher:")
                }

                if (length(arriba_found) > 0) {
                    result_text <- paste(result_text, "fusions confirmed by Arriba:")
                }

                if (length(fusioncatcher_found) > 0) {
                    print("fusioncatcher_found")
                    print(outfile)
                    print(fusioncatcher_found)
                    fusion_pairs <- unique(do.call(rbind, fusioncatcher_found))
                    print("fusion_pairs")
                    print(fusion_pairs)
                    #unique_fusion_pairs <- unique(apply(fusion_pairs, 1, function(x) paste(sort(x), collapse = " ")))
                    fusion_text <- paste(fusion_pairs, collapse = " ")
                    #fusion_text <- paste(apply(fusion_pairs, 1, function(x) paste(sort(x), collapse = " ")), collapse = "; ")
                    fusion_text <- gsub('["c()]', '', fusion_text)
                    fusion_text <- gsub(', ', ':', fusion_text)
                    fusion_text <- gsub(' ', '; ', fusion_text)
                }



                if (length(arriba_found) > 0) {
                    #fusion_text <- c(fusion_text, "Fusion (Arriba):")
                    for (i in 1:length(arriba_found)) {
                        fusion_text <- c(fusion_text, paste(arriba_found[[i]]$X.gene1, arriba_found[[i]]$gene2, collapse = " "))
                    }
                }

                return(paste(result_text, fusion_text, collapse = "\n"))
            } else {
                return("Gene fusion calling identified fusions not confirming the subtype allocation.")
            }


         } else {
            return("Subtype not found in the subtype data.")
         }
    }
}


result <- check_subtype_and_chromosome(subtype, chromosome_number, arriba_file,fusioncatcher_file)
cat(result, "\n")


# Create the output text
output_text <- paste("Based on the gene expression profile, the sample", sample_info,
                     "was allocated to the", subtype, "subtype with", confidence ,".\n" ,result)


# Write the output text to a CSV file
#output_file <- "out.csv"
write.csv(data.frame(Output=output_text), file=outfile, row.names=FALSE, quote=FALSE)

# Print a message to indicate that the process is complete
cat("Output file 'out.csv' has been generated.\n")
