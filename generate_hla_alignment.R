
# ------------------------ Generate HLA Alignment File ------------------------ 
#' 
#' This function generates a patient specific alignment file by querying IMGT/HLA Database.
#' 
#'
#' @param locus str: A, B, or C
#' @param allele1 str: 6-8 Digit code for patients allele 
#' @param allele2 str: 6-8 Digit code for patients allele 
#' @param file_path str: Output file name
#'
#' @return Returns nothing but saves the alignment file as a .txt to the file path location 
#'         
#'
#' @examples
#' generate_hla_alignment(locus = "A", allele1 = "01:02:01", allele2 = "07:01:01:01", file_path = "pat1_alignment_HLA_C.txt")
#' 
#' 
generate_hla_alignment = function(locus, allele1, allele2, file_path){
  
}

