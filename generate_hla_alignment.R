
# ------------------------ Generate HLA Alignment File ------------------------
#'
#' This function generates a patient specific alignment file by querying IMGT/HLA Database.
#'
#'
#' @param locus str: A, B, or C
#' @param allele1 str: 6-8 Digit code for patients allele
#' @param allele2 str: 6-8 Digit code for patients allele
#' @param file_path str: Output file name
#' @param reference str (optional): Custom reference allele. If not specified, uses default
#'        references (A*03:01:01:01, B*07:02:01:01, or C*07:02:01:01).
#'
#' @return Returns nothing but saves the alignment file as a .txt to the file path location
#'
#'
#' @examples
#' # Using default reference
#' generate_hla_alignment(locus = "A", allele1 = "01:02:01", allele2 = "07:01:01:01",
#'                        file_path = "pat1_alignment_HLA_")
#'
#' # Using custom reference
#' generate_hla_alignment(locus = "A", allele1 = "01:02:01", allele2 = "07:01:01:01",
#'                        file_path = "pat1_alignment_HLA_", reference = "A*02:01:01:01")
#'

source("Fetch_Alignment.R")

# Example 1: Using default reference allele
# -----------------------------------------
locus = "C"
allele1 = "01:02:01"
allele2 = "07:01:01:01"
allele1 = paste0(locus,"*",allele1)
allele2 = paste0(locus,"*",allele2)
output_path <- "results/test_"

# Call with default reference (C*07:02:01:01 for locus C)
generate_hla_alighnment(allele1, allele2, locus, filepath = output_path)


# Example 2: Using custom reference allele
# -----------------------------------------
# Uncomment to run with a custom reference:
#
# locus = "A"
# allele1 = "A*01:01:01:01"
# allele2 = "A*26:01:01:01"
# custom_ref = "A*02:01:01:01"
# output_path <- "results/custom_ref_"
#
# generate_hla_alighnment(allele1, allele2, locus, filepath = output_path, reference = custom_ref)
