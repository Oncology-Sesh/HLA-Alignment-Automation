library(SeshDevTools)
library(readr)

#------- Generate Alignments -----------
# HN_20_A
allele_ref = "A*03:01:01:01"
allele1= "A*02:01:01:01"
allele2= "A*11:01:01:01"

locus = "A"
output_path <- "tests/reload-test-data/alignments-auto/HN_20_"
generate_hla_alignment(allele1,allele2,locus, filepath = output_path,reference = allele_ref)



#------- Run Baseline Reload -----------
alignment_path = "tests/reload-test-data/alignments/HN_20_A.txt"
plasma_path = "tests/reload-test-data/counts/HN_20_A_readcounts_plasma.wig"
normal_path = "tests/reload-test-data/counts/HN_20_A_readcounts_normal.wig"

allele_alignment_file <- as.matrix(read_csv(alignment_path, col_names = F)) # this is to get the allele names.
ref_allele <- strsplit(allele_alignment_file[2],'      ')[[1]][1]
first_allele <- strsplit(allele_alignment_file[3],'      ')[[1]][1]
second_allele <- strsplit(allele_alignment_file[4],'      ')[[1]][1]

RELOAD_result_baseline = RELOAD(IMGTHLA_alignment_file_with_path = alignment_path,reference_allele = ref_allele,allele_1 = first_allele,allele_2 = second_allele,plasma_wig_file_with_path = plasma_path,normal_wig_file_with_path = normal_path,HLA_gene = locus )


#------- Run Comp Reload -----------

alignment_path_test = "tests/reload-test-data/alignments-auto/HN_20_A.txt"
plasma_path = "tests/reload-test-data/counts/HN_20_A_readcounts_plasma.wig"
normal_path = "tests/reload-test-data/counts/HN_20_A_readcounts_normal.wig"

allele_alignment_file <- as.matrix(read_csv(alignment_path_test, col_names = F)) # this is to get the allele names.
ref_allele_comp <- strsplit(allele_alignment_file[2],'   ')[[1]][1]
first_allele_comp <- strsplit(allele_alignment_file[3],'   ')[[1]][1]
second_allele_comp <- strsplit(allele_alignment_file[4],'   ')[[1]][1]

RELOAD_result_comp = RELOAD(IMGTHLA_alignment_file_with_path = alignment_path_test,reference_allele = ref_allele_comp,allele_1 = first_allele_comp,allele_2 = second_allele_comp,plasma_wig_file_with_path = plasma_path,normal_wig_file_with_path = normal_path,HLA_gene = locus )



