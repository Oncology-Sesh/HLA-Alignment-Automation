# Fetch_Alignment.R
# This file is the main alignment generation file.
# It takes input from the IMGT-IPD webpage.
# It aligns these sequence's using the DECIPHER package.
# Finally the alighnment is then printed as per the IPD-IMGT format. 

# ----------------------------- Libraries -----------------------------
library(httr2)
library(readr)
library(jsonlite)
library(DECIPHER)
library(tidyverse)
library(ggplot2)
library(Biostrings)
library(DECIPHER)
library(readxl)

# ----------------------------- Set Path -----------------------------
setwd("C:/Users/wells/Documents/Github/HLA-Alignment-Automation")

# -------------------- Data Aquisition and Wrangling ----------------

# Get Data from patient file
# Summary
# This code:
#         Starts with wide-format HLA data (A1, A2, etc.).
#         Converts it into long format (pivot_longer) to tidy the data.
#         Cleans up the Locus names (A1 → A).
#         Assigns an allele number (Allele-1, Allele-2).
#         Converts it back to wide, now in a clean tidy structure for each locus.

patient_dat_orig <- read_xlsx(path = "../data/HLA6d Andrew5.xlsx") #get csv of patient data
patient_dat <- as.data.frame(patient_dat_orig) %>% 
  pivot_longer(cols = 3:8, names_to = "Locus", values_to = "Allele_ID") %>% # Pivot longer all columns except the metadata
  mutate(Allele = rep(c("Allele-1", "Allele-2"), length.out = n()), # Create columns that assign Allele-1 and Allele-2 
         Locus = str_replace(Locus, "([ABC])[12]", "\\1"))%>% # replace A,B,C and 1&2 
  pivot_wider(.,names_from = Allele,values_from = Allele_ID)



# Fetch Allele information for testing 
locus = "C"
allele1 = "01:02:01"
allele2 = "07:01:01:01"
allele1 = paste0(locus,"*",allele1)
allele2 = paste0(locus,"*",allele2)



#################### Fetch Sequence and Perform Alignment ##########################################
# get IDs
print("Finding IDs...")
Allelelist_390 <- read_csv("data/Allelelist.390.txt",comment = "#", show_col_types = FALSE)
ID1 = Allelelist_390$AlleleID[Allelelist_390$Allele==allele1]
ID2 = Allelelist_390$AlleleID[Allelelist_390$Allele==allele2]

#ID2 = "HLA14780" #This was not present in original Allele reference for the input file being created. 

# reference allele logic from Mike's code

if (locus == "A"){
  IDref = "HLA00037"
}
if (locus == "B") {
  IDref = "HLA00132"
}
if (locus == "C"){
  IDref = "HLA00434"
}

allele_ref <- Allelelist_390$Allele[Allelelist_390$AlleleID==IDref] # Allele reference file

# Query HLA function
query_HLA <- function(ID){
  print(paste0("Querying HLA database for reference coding: ",ID)) # print output
  res = request(paste0("https://www.ebi.ac.uk/cgi-bin/ipd/api/allele/",ID))%>%req_perform() # scrape HLA data from IPD-IMGT
  json = fromJSON(rawToChar(res$body)) # Convert to JSON
  coding = json$sequence$coding # Get the Coding sequence
  return(coding) # Return the coding sequence
}

coding0 <- query_HLA(IDref) # Reference
coding1 <- query_HLA(ID1) # Allele-1
coding2 <- query_HLA(ID2) # Allel-2
codes = list(coding0,coding1,coding2) # List of all threee 


## Using DECIPHER to align sequences
seqs <- DNAStringSet(c(coding0,coding1,coding2)) # Create DNAStrings objct
names(seqs) <- c(IDref,ID1, ID2) #Provide names to the sequence #Assign proper names
aligned_seq <- AlignSeqs(seqs) #Align the sequences 

#################### Format Alignment files ##########################################

# Step 1: Convert to character matrix
char_mat <- as.matrix(aligned_seq)

# Step 2: Choose reference sequence (first sequence)
ref_seq <- char_mat[1, ]
ref_name <- rownames(char_mat)[1]

# Step 3: Annotate sequences: '-' for matches, base for mismatches, '*' gaps preserved or Truncated/Indels...etc
annot_mat <- char_mat
for (i in 2:nrow(char_mat)) {
  for (j in 1:ncol(char_mat)) {
    base <- char_mat[i, j]
    if (base == ref_seq[j] && base != "-") {
      annot_mat[i, j] <- "-"  # match
    } else if (base == "-") {
      annot_mat[i, j] <- "*"  # gap
    } else if (base == ref_seq[j] && base == "-") {
      annot_mat[i,j] <- "."
    } else {
      annot_mat[i, j] <- base  # mismatch
    }
  }
}

annot_mat[1, ] <- ref_seq  # reference unchanged
rownames(annot_mat) <-c(allele_ref,allele1,allele2)

#################### Insert Codon Space Function ########################
#' Insert space after every codon (3 nucleotide) in a sequence alignment block
#'
#' This function takes a character vector (typically a sequence of nucleotides)
#' and inserts a space after every `codon_size` elements, simulating codon boundaries.
#'
#' @param vec A character vector. Usually a vector of single-letter nucleotides (e.g., "A", "T", "G", ...).
#' @param global_start Integer (default = 1). The global start position for the first nucleotide.
#'        Used to maintain alignment if the input is part of a larger sequence.
#' @param codon_size Integer (default = 3). The number of characters per codon (usually 3 for DNA/RNA).
#'
#' @return A character vector the same as `vec`, but with `" "` (a space) inserted after every codon.
#'         Useful for formatting sequences for visual alignment or printing.
#'
#' @examples
#' insert_codon_space(c("A", "T", "G", "C", "G", "A"))
#' # Output: "A" "T" "G" " " "C" "G" "A"
#' 
insert_codon_space <- function(vec, global_start = 1, codon_size = 3) {
  # Create an empty character vector
  out <- character(0)
  
  
  #Loop over the sequence vector (each row of the annot_matrix)
  for (i in seq_along(vec)) {
    #Add i (each nucleotide) to the empty vector
    out <- c(out, vec[i])
    
    # If it's the end of a codon but not after the last char
    if (((global_start + i - 1) %% codon_size) == 0 && i != length(vec)) {
      
    # Insert a space
      out <- c(out, " ")
    }
  }
  out
}



#################### Codon Label Function ########################
#' Generate a ruler line for codon numbering in a sequence alignment block
#'
#' This function produces a character vector representing a single horizontal line
#' that labels every 5th codon within a given alignment block. The labels are left-aligned
#' within a fixed-width field (e.g. 4 characters wide), which can be stacked above
#' aligned DNA or protein sequences for display.
#'
#' @param block_start Integer. The **1-based global position** where this alignment block starts.
#'        For example, 1 for the first block, 76 for the second block (if 75-width), etc.
#'
#' @param block_end Integer. The **global end position** of this alignment block.
#'        For example, 75, 150, 225 depending on how many residues are shown per line.
#'
#' @param codon_size Integer. Number of characters per codon. Usually 3 (for DNA/RNA),
#'        but could be different for grouped amino acids, etc.
#'
#' @param codon_start_num Integer. The **codon number corresponding to position 1**
#'        in the global alignment. Could be negative or 1 depending on how alignment is anchored.
#'        E.g. -25 if position 1 corresponds to codon -25.
#'
#' @param codon_label_width Integer. The number of character positions to reserve for each label
#'        (default = 4). Wider labels allow space for triple-digit codon numbers like "100 ".
#'
#' @return A character vector of length (block_end - block_start + 1) containing the label line
#'         with codon numbers inserted at the appropriate positions (every 5th codon).
#'

get_alt_codon_label_line <- function(block_start,block_end,codon_size,codon_start_num, codon_label_width = 4) {
  
  # width of one alighnment block
  block_width   <- block_end - block_start + 1
  
  # create empty space vector
  codon_indices <- rep(" ", block_width)
  
  # Concat the codon_indices vector 
  cat(codon_indices)
  
  ## Integer offset: how many complete codons before this block?
  codon_offset  <- (block_start - 1) %/% codon_size
  
  ## DEBUG (SANITY CHECK): print the first codon number for the block as well as the boundries of a typical block. 
  message(sprintf("Block %4d–%-4d  first codon = %d",block_start, block_end,codon_start_num + codon_offset))
  
  # Iterate through the alignment block every 'codon_size' steps (i.e., one codon per loop)
  for (i in seq(1, block_width, by = codon_size)) {
    
    # Calculate codon number for this position (based on global codon_start_num)
    codon_num <- codon_start_num + ((i - 1) %/% codon_size)
    
    # Only label every 5th codon AND make sure the full label will fit in the block
    if (codon_num %% 5 == 0 && (i + codon_size - 1) <= block_width) {
      label <- sprintf("%-*s", codon_label_width, codon_num)
      codon_indices[i:(i + codon_label_width - 1)] <- strsplit(label, "")[[1]]
    }
  }
  codon_indices
}

# Not using this function anymore, remove it in final version.
# format_codon_labels_with_spacing <- function(codon_labels, codon_size = 3) {
#   out <- character(0)
#   for (i in seq_along(codon_labels)) {
#     out <- c(out, codon_labels[i])
#     if ((i %% codon_size) == 0 && i != length(codon_labels)) {
#       out <- c(out, " ")
#     }
#   }
#   paste0(out, collapse = "")
# }

#################### Generate Output Alignment files ##########################################
# Step 8: Print annotated alignment block-by-block

output_file <- paste0("results/Alignment_",locus,".txt")
sink(output_file)
cat("CODON-ALIGNED ANNOTATED ALIGNMENT vs Reference: HLA00434", "\n")
codon_size <- 3
base_offset<- -25
block_size <- 75
seq_names <- rownames(annot_mat)
seq_length <- ncol(annot_mat)
name_width <- 15

for (block_start in seq(1, seq_length, by = block_size)) {
  block_end <- min(block_start + block_size - 1, seq_length)
  block_width <- block_end - block_start + 1

  # Codon label line (based on global codon number)
  
  # Add spacing to codon labels to match sequence spacing
  
  ruler <- insert_codon_space(as.character((1:block_width) %% 10), global_start = block_start)
  cat(format("Ruler", width=name_width), paste(ruler, collapse = ""), "\n")
  
  # Print aligned with sequence labels
  # Codon number line (labels every 5th codon)
  codon_line <- get_alt_codon_label_line(
    block_start = block_start,
    block_end = block_end,
    codon_size = codon_size,
    codon_start_num = base_offset
  )
  #cat(format("Codon", width = name_width), paste(codon_line, collapse = ""))
  codon_line_spaced <- insert_codon_space(codon_line, global_start = block_start)
  cat(format("AA Codon", width = name_width), paste(codon_line_spaced, collapse = ""), "\n")
  
  for (i in 1:nrow(annot_mat)) {
    name <- format(seq_names[i], width = name_width)
    chars <- annot_mat[i, block_start:block_end]
    chars_barred <- insert_codon_space(chars, global_start = block_start)
    cat(name, paste(chars_barred, collapse = ""), "\n")
  }
  cat("\n")
  base_offset <- base_offset + (block_width %/% codon_size)
}

sink()
## 



