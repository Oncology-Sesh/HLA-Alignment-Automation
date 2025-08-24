library(readr)

# ------------ Set Dir ------------
setwd("C:/Users/wells/Documents/Github/HLA-Alignment-Automation")

# ------------- Comp File ---------
comp_file_path = paste0(getwd(),"/tests/comp/NIH_pt000419_HLA_C.txt")
comp_alignment <- read_csv(comp_file_path, col_names = FALSE,show_col_types = FALSE)
comp_alignment <- as.matrix(comp_alignment)
print(comp_alignment)

# ------------- Test File ---------
test_file_path = paste0(getwd(),"/tests/test/Alignment_C.txt")
test_alignment <- read_csv(test_file_path, col_names = FALSE,show_col_types = FALSE)
test_alignment <- as.matrix(test_alignment)
print(test_alignment)

