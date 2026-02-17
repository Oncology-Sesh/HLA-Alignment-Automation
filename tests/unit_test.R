library(readr)

# ------------ Set Dir ------------
setwd("C:/Users/wells/Documents/Github/HLA-Alignment-Automation")
source("./Fetch_Alignment.R")

# ------------ HLA A -------------
allele_ref = "A*03:01:01:01"
allele1 = "A*01:01:01:01"
allele2 = "A*26:01:01"

locus = "A"
output_path <- "tests/test/NIH_pt000419_HLA_"
generate_hla_alighnment(allele1,allele2,locus, filepath = output_path,reference = allele_ref)


# ------------ HLA B -------------
allele_ref = "B*07:02:01:01"
allele1 = "B*08:01:01"
allele2 = "B*27:05:02"
locus = "B"
output_path <- "tests/test/NIH_pt000419_HLA_"

generate_hla_alighnment(allele1,allele2,locus, filepath = output_path, reference = allele_ref)

# ------------- Comp File ---------
comp_file_path = paste0(getwd(),"/tests/comp/NIH_pt000419_HLA_B.txt")
comp_alignment <- read_csv(comp_file_path, col_names = FALSE,show_col_types = FALSE)
comp_alignment <- as.matrix(comp_alignment)
print(comp_alignment)

# ------------- Test File ---------
test_file_path = paste0(getwd(),"/tests/test/NIH_pt000419_HLA_B.txt")
test_alignment <- read_csv(test_file_path, col_names = FALSE,show_col_types = FALSE)
test_alignment <- as.matrix(test_alignment)
print(test_alignment)

