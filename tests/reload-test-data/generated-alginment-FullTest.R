# generated-alignment-FullTest.R
# Comprehensive test script for HLA alignment generation and RELOAD analysis
# Tests all samples in the alignments and counts folders

library(SeshDevTools)
library(readr)
library(dplyr)

# Set working directory to project root
# setwd("/Users/mwe11s/Documents/Projects/HLA-Alignment-Automation")

# Define paths
base_path <- "tests/reload-test-data"
alignments_path <- file.path(base_path, "alignments")
counts_path <- file.path(base_path, "counts")
auto_alignments_path <- file.path(base_path, "alignments-auto")

# Create auto alignments directory if it doesn't exist
if (!dir.exists(auto_alignments_path)) {
  dir.create(auto_alignments_path, recursive = TRUE)
}

#' Extract allele information from an alignment file
#'
#' @param alignment_file Path to the alignment file
#' @return A list with reference_allele, allele_1, allele_2
extract_alleles_from_alignment <- function(alignment_file) {
  lines <- readLines(alignment_file, n = 5)

  # Extract allele names from lines 2, 3, 4 (reference, allele1, allele2)
  # Format: " A*03:01:01:01      ATG GCC..."
  extract_allele_name <- function(line) {
    # Remove leading space and split by multiple spaces
    trimmed <- trimws(line)
    # Get the allele name (first part before multiple spaces)
    parts <- strsplit(trimmed, "\\s{2,}")[[1]]
    return(trimws(parts[1]))
  }

  ref_allele <- extract_allele_name(lines[2])
  allele_1 <- extract_allele_name(lines[3])
  allele_2 <- extract_allele_name(lines[4])

  return(list(
    reference_allele = ref_allele,
    allele_1 = allele_1,
    allele_2 = allele_2
  ))
}

#' Run RELOAD and extract results safely
#'
#' @param alignment_path Path to alignment file
#' @param plasma_path Path to plasma wig file
#' @param normal_path Path to normal wig file
#' @param locus HLA gene (A, B, or C)
#' @return List with pval and diff_mean for both alleles, or NA on error
run_reload_safe <- function(alignment_path, plasma_path, normal_path, locus) {
  tryCatch({
    # Extract alleles from alignment file
    alleles <- extract_alleles_from_alignment(alignment_path)

    # Run RELOAD
    result <- RELOAD(
      IMGTHLA_alignment_file_with_path = alignment_path,
      reference_allele = alleles$reference_allele,
      allele_1 = alleles$allele_1,
      allele_2 = alleles$allele_2,
      plasma_wig_file_with_path = plasma_path,
      normal_wig_file_with_path = normal_path,
      HLA_gene = locus
    )

    # Extract p-values and mean differences
    return(list(
      allele_1_pval = result[[1]]$p.value,
      allele_1_diff_mean = result[[1]]$estimate,
      allele_2_pval = result[[2]]$p.value,
      allele_2_diff_mean = result[[2]]$estimate,
      success = TRUE
    ))
  }, error = function(e) {
    message(paste("RELOAD error:", e$message))
    return(list(
      allele_1_pval = NA,
      allele_1_diff_mean = NA,
      allele_2_pval = NA,
      allele_2_diff_mean = NA,
      success = FALSE
    ))
  })
}

#' Generate alignment file safely
#'
#' @param allele_1 First allele
#' @param allele_2 Second allele
#' @param locus HLA gene (A, B, or C)
#' @param output_prefix Output file path prefix
#' @param reference Optional reference allele
#' @return TRUE if successful, FALSE otherwise
generate_alignment_safe <- function(allele_1, allele_2, locus, output_prefix, reference = NULL) {
  tryCatch({
    generate_hla_alignment(
      allele1 = allele_1,
      allele2 = allele_2,
      locus = locus,
      filepath = output_prefix,
      reference = reference
    )
    return(TRUE)
  }, error = function(e) {
    message(paste("Alignment generation error:", e$message))
    return(FALSE)
  })
}

# Get list of all alignment files
alignment_files <- list.files(alignments_path, pattern = "\\.txt$", full.names = FALSE)
# Exclude list files
alignment_files <- alignment_files[!grepl("list", alignment_files, ignore.case = TRUE)]

# Get list of all count files (unique sample-gene combinations)
count_files <- list.files(counts_path, pattern = "_readcounts_plasma\\.wig$", full.names = FALSE)
count_samples <- gsub("_readcounts_plasma\\.wig$", "", count_files)

# Find samples that have both alignment and count files
# Parse sample info from filenames
parse_sample_info <- function(filename) {
  # Handle different naming patterns:
  # MBC_191_A.txt -> sample_id = MBC_191, gene = A
  # PC_17_HLA_A.txt -> sample_id = PC_17, gene = A
  # HN_20_A.txt -> sample_id = HN_20, gene = A

  # Remove .txt extension
  base <- gsub("\\.txt$", "", filename)

  # Check if it has HLA_ pattern
  if (grepl("_HLA_", base)) {
    parts <- strsplit(base, "_HLA_")[[1]]
    sample_id <- parts[1]
    gene <- parts[2]
  } else {
    # Last part after _ is the gene
    parts <- strsplit(base, "_")[[1]]
    gene <- parts[length(parts)]
    sample_id <- paste(parts[-length(parts)], collapse = "_")
  }

  return(list(sample_id = sample_id, gene = gene))
}

# Initialize results dataframe
results <- data.frame(
  sample_id = character(),
  gene = character(),
  allele_num = integer(),
  allele_code = character(),
  alignment_file = logical(),
  reference_pval = numeric(),
  generated_pval = numeric(),
  reference_diff_mean = numeric(),
  generated_diff_mean = numeric(),
  match = character(),
  stringsAsFactors = FALSE
)

# Process each alignment file
message("Starting full test of all samples...")
message(paste("Found", length(alignment_files), "alignment files"))

for (align_file in alignment_files) {
  # Parse sample info
  info <- parse_sample_info(align_file)
  sample_id <- info$sample_id
  gene <- info$gene

  message(paste("\n=== Processing:", sample_id, gene, "==="))

  # Check if count files exist for this sample
  # Handle different naming conventions:
  # - MBC/HN samples: MBC_191_A_readcounts_plasma.wig / MBC_191_A_readcounts_normal.wig
  # - PC samples: PC_17_A_readcounts_cf.wig / PC_17_A_readcounts_normal.wig (no HLA_, uses cf instead of plasma)

  # Try plasma/normal naming first (MBC/HN samples)
  plasma_file <- paste0(sample_id, "_", gene, "_readcounts_plasma.wig")
  normal_file <- paste0(sample_id, "_", gene, "_readcounts_normal.wig")

  plasma_path <- file.path(counts_path, plasma_file)
  normal_path <- file.path(counts_path, normal_file)

  # If not found, try cf/normal naming (PC samples)
  if (!file.exists(plasma_path)) {
    plasma_file <- paste0(sample_id, "_", gene, "_readcounts_cf.wig")
    plasma_path <- file.path(counts_path, plasma_file)
  }

  if (!file.exists(plasma_path) || !file.exists(normal_path)) {
    message(paste("  Skipping - count files not found"))
    message(paste("    Tried:", plasma_file, "and", normal_file))
    next
  }

  # Get reference alignment path
  ref_alignment_path <- file.path(alignments_path, align_file)

  # Extract alleles from reference alignment
  alleles <- tryCatch({
    extract_alleles_from_alignment(ref_alignment_path)
  }, error = function(e) {
    message(paste("  Error extracting alleles:", e$message))
    NULL
  })

  if (is.null(alleles)) next

  message(paste("  Reference:", alleles$reference_allele))
  message(paste("  Allele 1:", alleles$allele_1))
  message(paste("  Allele 2:", alleles$allele_2))

  # Generate new alignment
  output_prefix <- file.path(auto_alignments_path, paste0(sample_id, "_"))
  gen_success <- generate_alignment_safe(
    allele_1 = alleles$allele_1,
    allele_2 = alleles$allele_2,
    locus = gene,
    output_prefix = output_prefix,
    reference = alleles$reference_allele
  )

  gen_alignment_path <- paste0(output_prefix, gene, ".txt")
  alignment_file_created <- gen_success && file.exists(gen_alignment_path)

  message(paste("  Alignment generated:", alignment_file_created))

  # Run RELOAD on reference alignment
  message("  Running RELOAD on reference alignment...")
  ref_results <- run_reload_safe(ref_alignment_path, plasma_path, normal_path, gene)

  # Run RELOAD on generated alignment (if it exists)
  gen_results <- list(
    allele_1_pval = NA,
    allele_1_diff_mean = NA,
    allele_2_pval = NA,
    allele_2_diff_mean = NA,
    success = FALSE
  )

  if (alignment_file_created) {
    message("  Running RELOAD on generated alignment...")
    gen_results <- run_reload_safe(gen_alignment_path, plasma_path, normal_path, gene)
  }

  # Compare results and determine match status
  compare_values <- function(ref_val, gen_val, tolerance = 1e-6) {
    if (is.na(ref_val) || is.na(gen_val)) return(FALSE)
    return(abs(ref_val - gen_val) < tolerance)
  }

  # Allele 1 results
  allele_1_pval_match <- compare_values(ref_results$allele_1_pval, gen_results$allele_1_pval)
  allele_1_mean_match <- compare_values(ref_results$allele_1_diff_mean, gen_results$allele_1_diff_mean)
  allele_1_match <- if (allele_1_pval_match && allele_1_mean_match) "PASS" else "FAIL"

  # Allele 2 results
  allele_2_pval_match <- compare_values(ref_results$allele_2_pval, gen_results$allele_2_pval)
  allele_2_mean_match <- compare_values(ref_results$allele_2_diff_mean, gen_results$allele_2_diff_mean)
  allele_2_match <- if (allele_2_pval_match && allele_2_mean_match) "PASS" else "FAIL"

  # Add allele 1 row
  results <- rbind(results, data.frame(
    sample_id = sample_id,
    gene = gene,
    allele_num = 1,
    allele_code = alleles$allele_1,
    alignment_file = alignment_file_created,
    reference_pval = ref_results$allele_1_pval,
    generated_pval = gen_results$allele_1_pval,
    reference_diff_mean = as.numeric(ref_results$allele_1_diff_mean),
    generated_diff_mean = as.numeric(gen_results$allele_1_diff_mean),
    match = allele_1_match,
    stringsAsFactors = FALSE
  ))

  # Add allele 2 row
  results <- rbind(results, data.frame(
    sample_id = sample_id,
    gene = gene,
    allele_num = 2,
    allele_code = alleles$allele_2,
    alignment_file = alignment_file_created,
    reference_pval = ref_results$allele_2_pval,
    generated_pval = gen_results$allele_2_pval,
    reference_diff_mean = as.numeric(ref_results$allele_2_diff_mean),
    generated_diff_mean = as.numeric(gen_results$allele_2_diff_mean),
    match = allele_2_match,
    stringsAsFactors = FALSE
  ))

  message(paste("  Allele 1 match:", allele_1_match))
  message(paste("  Allele 2 match:", allele_2_match))
}

# Print summary
message("\n========== TEST SUMMARY ==========")
message(paste("Total tests:", nrow(results)))
message(paste("Alignments created:", sum(results$alignment_file)))
message(paste("PASS:", sum(results$match == "PASS", na.rm = TRUE)))
message(paste("FAIL:", sum(results$match == "FAIL", na.rm = TRUE)))

# Save results to CSV
output_csv <- file.path(base_path, "test_results_summary.csv")
write_csv(results, output_csv)
message(paste("\nResults saved to:", output_csv))

# Print the results dataframe
print(results)

# Return summary statistics
summary_stats <- list(
  total_tests = nrow(results),
  alignments_created = sum(results$alignment_file),
  pass_count = sum(results$match == "PASS", na.rm = TRUE),
  fail_count = sum(results$match == "FAIL", na.rm = TRUE),
  pass_rate = sum(results$match == "PASS", na.rm = TRUE) / nrow(results) * 100
)

message(paste("\nPass rate:", round(summary_stats$pass_rate, 2), "%"))
