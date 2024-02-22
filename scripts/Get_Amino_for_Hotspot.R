#options(warn=-1)
suppressMessages(library(Biostrings))



args <- commandArgs(trailingOnly = TRUE)
#print(args)
pysamstats_directory <- args[1]
ctat_file <- args[2]
output_dir <- args[3]

mutations <- data.frame(
  Gene = c("ZEB2", "ZEB2", "KRAS", "KRAS", "KRAS", "KRAS", "NRAS", "NRAS", "NRAS", "NRAS", "FLT3", "FLT3", "FLT3", "FLT3", "FLT3", "FLT3", "FLT3", "FLT3", "FLT3",
           "FLT3", "FLT3", "FLT3", "FLT3", "FLT3", "FLT3", "FLT3", "FLT3", "FLT3",
           "FLT3", "FLT3", "FLT3", "FLT3", "FLT3", "FLT3", "FLT3", "FLT3", "FLT3", "PAX5"),
  Hotspot = c("H1038", "Q1072", "G12", "G13", "Q16", "A146", "G12", "G13", "Q61", "A146", "P857", "V843", "Y842", "N841", "D839", "M837", "I836", "D835", "R834",
              "A680", "N676", "A627", "K623", "Y599", "R595", "V592", "Y589", "N587",
              "G583", "Q580", "V579", "Q577", "L576", "E573", "Y572", "V491", "S446", "P80R"),
  Chromosome = c("2", "2", "12", "12", "12", "12", "1", "1", "1", "1", "13", "13", "13", "13", "13", "13", "13", "13", "13",
                 "13", "13", "13", "13", "13", "13", "13", "13", "13",
                 "13", "13", "13", "13", "13", "13", "13", "13", "13", "9"),
  Start = c("144389982", "144389880" , "25245349", "25245346", "25227341", "25225626", "114716125", "114716122", "114713907", "114709581", "28015672", "28018479", "28018482", "28018485", "28018491", "28018497", "28018500", "28018503", "28018506",
            "28028191", "28028203", "28033948", "28033960", "28034122", "28034134", "28034143", "28034152", "28034158",
            "28034170", "28034179", "28034182", "28034188", "28034191", "28034200", "28034203", "28035619", "28036015", "37015167"),
  End = c("144389984", "144389882", "25245351", "25245348", "25227343", "25225628", "114716127", "114716124", "114713909", "114709583", "28015674", "28018481", "28018484", "28018487", "28018493", "28018499", "28018502", "28018505", "28018508",
          "28028193", "28028205", "28033950", "28033962", "28034124", "28034136", "28034145", "28034154", "28034160",
          "28034172", "28034181", "28034184", "28034190", "28034193", "28034202", "28034205", "28035621", "28036017", "37015169")
)

mutations$Start <- as.integer(mutations$Start)
mutations$End <- as.integer(mutations$End)

find_CHROM_row_header <- function(file_path) {
  con <- file(file_path, "r")
  line_num <- 0
  while (TRUE) {
    line <- readLines(con, n = 1)
    if (length(line) == 0) {
      break  # End of file
    }
    line_num <- line_num + 1
    if (grepl("#CHROM", line)) {
      print("line_num")
      print(line_num)
      return(line_num)
    }
  }
  return(NULL)  # #CHROM not found
}

# Find the row number where #CHROM occurs
#skip <- find_CHROM_row_header(file_path)


#vcf <- read.csv(ctat_file, skip=skip-1, header=TRUE, sep="\t")



find_CHROM_row <- function(file_path, chrom_val, pos_val) {
  skip <- find_CHROM_row_header(file_path)
  skip <- skip - 1
  
  # Einlesen des Dateiinhalts als DataFrame
  df <- read.csv(file_path, header = FALSE, skip = skip, stringsAsFactors = FALSE, sep = "\t")
  
  # Setzen der Spaltennamen
  headers <- c("X.CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "Sample_Values")
  colnames(df) <- headers
  print(names(df))
  print("vor subset")
  # Subset basierend auf den Bedingungen
  filtered_rows <- subset(df, X.CHROM == chrom_val & POS >= pos_val & POS <= pos_val + 2)
  
  # Überprüfung, ob Zeilen gefunden wurden
  if (nrow(filtered_rows) > 0) {
    #found_info <- filtered_rows$INFO  # Annahme: Es wird nur der erste gefundene Wert zurückgegeben
    found_info <- filtered_rows
    print(found_info)
    return(found_info)
  } else {
    print("Keine passenden Zeilen gefunden")
  }
}



# TODO:
# 1. Mutations Liste mit allen benötigten mutationen ergänzen. Evtl. Trennung der liste für sense und antisense?



file_list <- list.files(path = pysamstats_directory, full.names = TRUE)

for (file in file_list) {
  if (grepl("\\.tsv$", file, ignore.case = TRUE)) {
    #print("not IKZF1")
    data <- read.csv(file, sep="\t", header = TRUE)
    if(NROW(data) > 2) {
      
      #print("rows")
      #print(NROW(data))
      file_name <- basename(file)
      #print(paste("Extracted info from filename:", file_name))
      pattern <- "([^_]+_[^_]+)\\.tsv$" # Pattern to match the text between the first and second underscore before .tsv
      extracted_info <- sub(".*_([^_]+_[^_]+)\\.tsv$", "\\1", file_name)
      #print(paste("Extracted info from filename:", extracted_info))
      output_file <- paste0(extracted_info, "_output_file.csv")
      gene <- substr(extracted_info, 1, regexpr("_", extracted_info) - 1)
      hotspot <- substr(extracted_info, regexpr("_", extracted_info) + 1, nchar(extracted_info))
      # Get the corresponding information from mutations table
      entry <- subset(mutations, Gene == gene & Hotspot == hotspot)
      chromosome <- entry$Chromosome
      start <- as.numeric(entry$Start)
      end <- as.numeric(entry$End)
      
      # Find and save "new bases" in data
      data$new_base <- ""
      for(row in 1:NROW(data)){
        # Iterate data rows, and check if mismatches > 0 and mismatches_pp/matches_pp > 0.05
        if(data$mismatches_pp[row] > 0 & (data$mismatches_pp[row] / data$matches_pp[row]) >= 0.05){
          # Get ref and switch case based on Ref base and maximum Character_pp, append result to data
          ref <- data$ref[row]
          res <- switch(
            ref,
            "A" = {if(colnames(data[row, c(17,19,21,23)])[apply(data[row, c(17,19,21,23)], 1, which.max)] == "C_pp"){"C"} else if(colnames(data[row, c(17,19,21,23)])[apply(data[row, c(17,19,21,23)], 1, which.max)] == "T_pp"){"T"} else if(colnames(data[row, c(17,19,21,23)])[apply(data[row, c(17,19,21,23)], 1, which.max)] == "G_pp"){"G"} else if(colnames(data[row, c(17,19,21,23)])[apply(data[row, c(17,19,21,23)], 1, which.max)] == "N_pp"){"N"} else{"A"}},
            "C" = {if(colnames(data[row, c(15,19,21,23)])[apply(data[row, c(15,19,21,23)], 1, which.max)] == "A_pp"){"A"} else if(colnames(data[row, c(15,19,21,23)])[apply(data[row, c(15,19,21,23)], 1, which.max)] == "T_pp"){"T"} else if(colnames(data[row, c(15,19,21,23)])[apply(data[row, c(15,19,21,23)], 1, which.max)] == "G_pp"){"G"} else if(colnames(data[row, c(15,19,21,23)])[apply(data[row, c(15,19,21,23)], 1, which.max)] == "N_pp"){"N"} else{"C"}},
            "T" = {if(colnames(data[row, c(15,17,21,23)])[apply(data[row, c(15,17,21,23)], 1, which.max)] == "C_pp"){"C"} else if(colnames(data[row, c(15,17,21,23)])[apply(data[row, c(15,17,21,23)], 1, which.max)] == "A_pp"){"A"} else if(colnames(data[row, c(15,17,21,23)])[apply(data[row, c(15,17,21,23)], 1, which.max)] == "G_pp"){"G"} else if(colnames(data[row, c(15,17,21,23)])[apply(data[row, c(15,17,21,23)], 1, which.max)] == "N_pp"){"N"} else{"T"}},
            "G" = {if(colnames(data[row, c(15,17,19,23)])[apply(data[row, c(15,17,19,23)], 1, which.max)] == "C_pp"){"C"} else if(colnames(data[row, c(15,17,19,23)])[apply(data[row, c(15,17,19,23)], 1, which.max)] == "T_pp"){"T"} else if(colnames(data[row, c(15,17,19,23)])[apply(data[row, c(15,17,19,23)], 1, which.max)] == "A_pp"){"A"} else if(colnames(data[row, c(15,17,19,23)])[apply(data[row, c(15,17,19,23)], 1, which.max)] == "N_pp"){"N"} else{"G"}},
            "N" = {if(colnames(data[row, c(15,17,19,21)])[apply(data[row, c(15,17,19,21)], 1, which.max)] == "C_pp"){"C"} else if(colnames(data[row, c(15,17,19,21)])[apply(data[row, c(15,17,19,21)], 1, which.max)] == "T_pp"){"T"} else if(colnames(data[row, c(15,17,19,21)])[apply(data[row, c(15,17,19,21)], 1, which.max)] == "A_pp"){"A"} else if(colnames(data[row, c(15,17,19,21)])[apply(data[row, c(15,17,19,21)], 1, which.max)] == "G_pp"){"G"} else{"N"}},
          )
          data[row, 24] <- res
          percentage <- 0
          v <- list("A", "C", "T", "G")
          #Reverse data
          current_block <- data[nrow(data):1, ]
          bases <- list()
          
          
          
          # Check whether new bases are available and calculate percentage:
          for(k in 1:nrow(current_block)){
            if(current_block$new_base[k] %in% v & current_block$reads_pp[k] >= 50){
              if (current_block$new_base[k] == TRUE){
                print(" TRUE")
                bases <- append(bases, "T")
              } else{
                bases <- append(bases, current_block$new_base[k])
              }
              
              col <- sprintf("%s_pp", current_block$new_base[k])
              percentage <- switch(
                col,
                "A_pp" = {current_block$A_pp[k] * 100 / current_block$reads_pp[k]},
                "C_pp" = {current_block$C_pp[k] * 100 / current_block$reads_pp[k]},
                "T_pp" = {current_block$T_pp[k] * 100 / current_block$reads_pp[k]},
                "G_pp" = {current_block$G_pp[k] * 100 / current_block$reads_pp[k]},
              )
              
            }else{
              if (current_block$ref[k] == TRUE){
                bases <- append(bases, "T")
              }else {
                bases <- append(bases, current_block$ref[k])
              }
              
            }
          }
          print("bases")
          print(bases)
          bases <- paste(unlist(bases), collapse="")
          
          reverse_bases <- toString(complement(DNAString(bases)))
          
          newAminoAcid <- GENETIC_CODE[[reverse_bases]]
          
          
          if(percentage > 0){
            write.csv(data.frame(Gene = character(), Hotspot = character(), Chromosome = character(), Start = numeric(), End = numeric(), New_Aminoacid = character(), Percentage = numeric()), file = paste0(output_dir,"/", output_file), row.names = FALSE)
            selected_columns <- current_block[c("chrom", "pos", "ref", "reads_pp", "mismatches_pp", "deletions_pp", "insertions_pp", "A_pp", "C_pp", "T_pp", "G_pp", "N_pp", "new_base")]
            selected_columns_right_order <- selected_columns[nrow(selected_columns):1, ]
            write.table(selected_columns_right_order, file = paste0(output_dir,"/",extracted_info, "_with_bases.tsv"), row.names = FALSE)
            result <- data.frame(Gene = gene, Hotspot = hotspot, Chromosome = chromosome, Start = start, End = end, New_Aminoacid = newAminoAcid, Percentage = percentage)
            write.table(result, paste0(output_dir,"/", output_file), append = TRUE, sep = ",", col.names = FALSE, row.names = FALSE)
            
            result_row <- find_CHROM_row(ctat_file, chromosome, start)
            if(!is.null(result_row)){
              print("result_row:")
              print(result_row)
              result_df <- as.data.frame(result_row, stringsAsFactors = FALSE)
              col_names <- c("X.CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "Sample")
              #result_df <- data.frame(t(result_row), stringsAsFactors = FALSE)
              #colnames(result_df) <- col_names
              write.table(result_df, file = paste0(output_dir, "/", extracted_info, "_ctat_result.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, append = TRUE)
              write.table(result_df, file = paste0(output_dir, "/", extracted_info, "_ctat_result_1.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, append = TRUE)
              
              #output <- paste("CTAT INFO data for chromosome", chromosome, "and", start, ":", result_row)
              #writeLines(output, paste0(output_dir,"/",extracted_info, "_ctat_result.tsv"))
              #cat("CTAT INFO data for chromosome", chromosome, "at position", start, ":\n", file = paste0(output_dir, "/", extracted_info, "_ctat_result.tsv"))
              #cat(paste(result_row, collapse = "\n"), file = paste0(output_dir, "/", extracted_info, "_ctat_result.tsv"), append = TRUE)
            } else {
              print("result_row empty")
              #output <- paste("Zeile nicht gefunden", chromosome, "and", start)
              #writeLines(output, paste0(output_dir,"/",extracted_info, "_ctat_result.tsv"))
              cat("Empty result for chromosome", chromosome, "at position", start, "\n", file = paste0(output_dir, "/", extracted_info, "_ctat_result.tsv"))
            }
            percentage <- 0
          }else{
            # Nothing
          }
          
        }else {
          # Nothing
        }
      }
    }
  }
  
}

for (file in file_list) {
  if (grepl("\\.csv$", file, ignore.case = TRUE)) {
    print("IKZF1")
    data <- read.csv(file, sep="\t", header = TRUE)
    file_name <- basename(file)
    pattern <- "([^_]+_[^_]+)\\.tsv$" # Pattern to match the text between the first and second underscore before .tsv
    extracted_info <- sub(".*_([^_]+_[^_]+)\\.tsv$", "\\1", file_name)
    print(paste("Extracted info from filename:", extracted_info))
    output_file <- paste0(extracted_info, "_output_file.csv")
    gene <- substr(extracted_info, 1, regexpr("_", extracted_info) - 1)
    hotspot <- substr(extracted_info, regexpr("_", extracted_info) + 1, nchar(extracted_info))
    # Get the corresponding information from mutations table
    #entry <- subset(mutations, Gene == gene & Hotspot == hotspot)
    chromosome <- "chr7"
    start <- 50382593
    end <- 50382596
    gene <- "IKZF1"
    hotspot <- "Y159N"
    # Find and save "new bases" in data
    data$new_base <- ""
    for(row in 1:NROW(data)){
      # Iterate data rows, and check if mismatches > 0 and mismatches_pp/matches_pp > 0.05
      if(data$mismatches_pp[row] > 0 & (data$mismatches_pp[row] / data$matches_pp[row]) >= 0.05){
        # Get ref and switch case based on Ref base and maximum Character_pp, append result to data
        ref <- data$ref[row]
        res <- switch(
          ref,
          "A" = {if(colnames(data[row, c(17,19,21,23)])[apply(data[row, c(17,19,21,23)], 1, which.max)] == "C_pp"){"C"} else if(colnames(data[row, c(17,19,21,23)])[apply(data[row, c(17,19,21,23)], 1, which.max)] == "T_pp"){"T"} else if(colnames(data[row, c(17,19,21,23)])[apply(data[row, c(17,19,21,23)], 1, which.max)] == "G_pp"){"G"} else if(colnames(data[row, c(17,19,21,23)])[apply(data[row, c(17,19,21,23)], 1, which.max)] == "N_pp"){"N"} else{"A"}},
          "C" = {if(colnames(data[row, c(15,19,21,23)])[apply(data[row, c(15,19,21,23)], 1, which.max)] == "A_pp"){"A"} else if(colnames(data[row, c(15,19,21,23)])[apply(data[row, c(15,19,21,23)], 1, which.max)] == "T_pp"){"T"} else if(colnames(data[row, c(15,19,21,23)])[apply(data[row, c(15,19,21,23)], 1, which.max)] == "G_pp"){"G"} else if(colnames(data[row, c(15,19,21,23)])[apply(data[row, c(15,19,21,23)], 1, which.max)] == "N_pp"){"N"} else{"C"}},
          "T" = {if(colnames(data[row, c(15,17,21,23)])[apply(data[row, c(15,17,21,23)], 1, which.max)] == "C_pp"){"C"} else if(colnames(data[row, c(15,17,21,23)])[apply(data[row, c(15,17,21,23)], 1, which.max)] == "A_pp"){"A"} else if(colnames(data[row, c(15,17,21,23)])[apply(data[row, c(15,17,21,23)], 1, which.max)] == "G_pp"){"G"} else if(colnames(data[row, c(15,17,21,23)])[apply(data[row, c(15,17,21,23)], 1, which.max)] == "N_pp"){"N"} else{"T"}},
          "G" = {if(colnames(data[row, c(15,17,19,23)])[apply(data[row, c(15,17,19,23)], 1, which.max)] == "C_pp"){"C"} else if(colnames(data[row, c(15,17,19,23)])[apply(data[row, c(15,17,19,23)], 1, which.max)] == "T_pp"){"T"} else if(colnames(data[row, c(15,17,19,23)])[apply(data[row, c(15,17,19,23)], 1, which.max)] == "A_pp"){"A"} else if(colnames(data[row, c(15,17,19,23)])[apply(data[row, c(15,17,19,23)], 1, which.max)] == "N_pp"){"N"} else{"G"}},
          "N" = {if(colnames(data[row, c(15,17,19,21)])[apply(data[row, c(15,17,19,21)], 1, which.max)] == "C_pp"){"C"} else if(colnames(data[row, c(15,17,19,21)])[apply(data[row, c(15,17,19,21)], 1, which.max)] == "T_pp"){"T"} else if(colnames(data[row, c(15,17,19,21)])[apply(data[row, c(15,17,19,21)], 1, which.max)] == "A_pp"){"A"} else if(colnames(data[row, c(15,17,19,21)])[apply(data[row, c(15,17,19,21)], 1, which.max)] == "G_pp"){"G"} else{"N"}},
        )
        data[row, 24] <- res
        percentage <- 0
        v <- list("A", "C", "T", "G")
        #Reverse data
        current_block <- data
        bases <- list()
        
        
        
        # Check whether new bases are available and calculate percentage:
        for(k in 1:nrow(current_block)){
          if(current_block$new_base[k] %in% v & current_block$reads_pp[k] >= 50){
            bases <- append(bases, current_block$new_base[k])
            col <- sprintf("%s_pp", current_block$new_base[k])
            percentage <- switch(
              col,
              "A_pp" = {current_block$A_pp[k] * 100 / current_block$reads_pp[k]},
              "C_pp" = {current_block$C_pp[k] * 100 / current_block$reads_pp[k]},
              "T_pp" = {current_block$T_pp[k] * 100 / current_block$reads_pp[k]},
              "G_pp" = {current_block$G_pp[k] * 100 / current_block$reads_pp[k]},
            )
            
          }else{
            bases <- append(bases, current_block$ref[k])
          }
        }
        
        bases <- paste(unlist(bases), collapse="")
        bases <- toString((DNAString(bases)))
        print("bases")
        print(bases)
        #reverse_bases <- toString(complement(DNAString(bases)))
        #newAminoAcid <- GENETIC_CODE[[reverse_bases]]
        newAminoAcid <- GENETIC_CODE[[bases]]
        print("newAminoAcid")
        print(newAminoAcid)
        
        
        
        if(percentage > 0){
          write.csv(data.frame(Gene = character(), Hotspot = character(), Chromosome = character(), Start = numeric(), End = numeric(), New_Aminoacid = character(), Percentage = numeric()), file = paste0(output_dir,"/", output_file), row.names = FALSE)
          selected_columns <- current_block[c("chrom", "pos", "ref", "reads_pp", "mismatches_pp", "deletions_pp", "insertions_pp", "A_pp", "C_pp", "T_pp", "G_pp", "N_pp", "new_base")]
          selected_columns_right_order <- selected_columns[nrow(selected_columns):1, ]
          write.table(selected_columns_right_order, file = paste0(output_dir,"/",extracted_info, "_with_bases.tsv"), row.names = FALSE)
          chromosome <- "chr7"
          start <- 50382593
          end <- 50382596
          gene <- "IKZF1"
          hotspot <- "Y159N"
          result <- data.frame(Gene = gene, Hotspot = hotspot, Chromosome = chromosome, Start = start, End = end, New_Aminoacid = newAminoAcid, Percentage = percentage)
          write.table(result, paste0(output_dir,"/", output_file), append = TRUE, sep = ",", col.names = FALSE, row.names = FALSE)
          
          result_row <- find_CHROM_row(ctat_file, chromosome, start)
          if(!is.null(result_row)){
            #info_column_index <- which(headers =="INFO")
            #info_data <- result_row[info_column_index]
            output <- paste("CTAT INFO data for", chromosome, "and", start, ":", result_row)
            writeLines(output, paste0(output_dir,"/",extracted_info, "_ctat_result.tsv"))
          } else {
            output <- paste("Zeile nicht gefunden IKZF1", chromosome, "and", start)
            writeLines(output, paste0(output_dir,"/",extracted_info, "_ctat_result.tsv"))
          }
          percentage <- 0
        }else{
          # Nothing
        }
        
      }else {
        # Nothing
      }
    }
    
  }
  
}
