# HLA-Alignment-Automation

Automation script for generating codon-aligned HLA alignment files for the RELOAD algorithm.

## Overview

This R-based tool automates the generation of codon-aligned Human Leukocyte Antigen (HLA) alignment files by:

1. Fetching HLA sequences from the IMGT/IPD (Immuno Polymorphism Database) API
2. Performing multiple sequence alignment using DECIPHER
3. Generating formatted alignment output with reference-based annotation

The tool is designed for immunology and oncology research where precise HLA typing and alignment is critical for patient-specific profiling and therapeutic target identification.

## Prerequisites

- R >= 3.6.0 (recommended: R 4.0+)
- Internet connection (for IMGT/IPD API access)

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/Oncology-Sesh/HLA-Alignment-Automation.git
   cd HLA-Alignment-Automation
   ```

2. Install required R packages:
   ```R
   install.packages(c("httr2", "readr", "jsonlite", "tidyverse", "ggplot2", "readxl"))

   # Bioconductor packages
   if (!require("BiocManager", quietly = TRUE))
       install.packages("BiocManager")
   BiocManager::install(c("DECIPHER", "Biostrings"))
   ```

## Usage

### Basic Usage

```R
# Load the module
source("Fetch_Alignment.R")

# Generate alignment for two HLA alleles at a locus (using default reference)
generate_hla_alighnment(
  allele1 = "A*01:01:01",       # First allele (with locus prefix)
  allele2 = "A*26:01:01:01",    # Second allele (with locus prefix)
  locus = "A",                   # Locus: A, B, or C
  filepath = "results/patient_" # Output path prefix
)
# Creates: results/patient_A.txt
```

### Using a Custom Reference Allele

```R
source("Fetch_Alignment.R")

# Specify a custom reference allele instead of the default
generate_hla_alighnment(
  allele1 = "A*01:01:01",
  allele2 = "A*26:01:01:01",
  locus = "A",
  filepath = "results/patient_",
  reference = "A*02:01:01:01"    # Custom reference allele
)
```

### Complete Patient Workflow

```R
source("Fetch_Alignment.R")

# Generate alignments for all three HLA loci (using default references)
generate_hla_alighnment("A*01:01:01", "A*26:01:01:01", "A", "results/patient_123_HLA_")
generate_hla_alighnment("B*08:01:01", "B*27:05:02", "B", "results/patient_123_HLA_")
generate_hla_alighnment("C*01:02:01", "C*07:01:01:01", "C", "results/patient_123_HLA_")

# Or with custom references
generate_hla_alighnment("A*01:01:01", "A*26:01:01:01", "A", "results/patient_123_HLA_",
                        reference = "A*02:01:01:01")
```

### Running the Example Script

```R
source("generate_hla_alignment.R")
```

## Project Structure

```
HLA-Alignment-Automation/
├── Fetch_Alignment.R          # Main module with core functions
├── generate_hla_alignment.R   # Example/wrapper script
├── data/
│   └── Allelelist.390.txt     # Allele ID database (IPD-IMGT/HLA 3.9.0)
├── results/                   # Generated alignment outputs
├── tests/
│   ├── comp/                  # Reference comparison files
│   ├── test/                  # Test output files
│   └── unit_test.R            # Unit test script
└── README.md
```

## Core Functions

### `query_HLA(ID)`
Queries the EBI IPD-IMGT/HLA API to retrieve the coding sequence for a given allele ID.

### `insert_codon_space(vec, global_start, codon_size)`
Formats sequence output by inserting spaces at codon boundaries (every 3 nucleotides).

### `get_alt_codon_label_line(block_start, block_end, codon_size, codon_start_num, codon_label_width)`
Generates ruler lines labeling codon positions at regular intervals.

### `generate_hla_alighnment(allele1, allele2, locus, filepath, reference = NULL)`
Main function that orchestrates the alignment workflow:
- Looks up allele IDs from the database
- Fetches sequences from the API
- Performs multiple sequence alignment
- Annotates and formats the output

**Parameters:**
- `allele1`: First patient allele (e.g., "A*01:01:01")
- `allele2`: Second patient allele (e.g., "A*26:01:01:01")
- `locus`: HLA locus - "A", "B", or "C"
- `filepath`: Output file path prefix
- `reference`: (Optional) Custom reference allele. If NULL, uses default references per locus

## Output Format

The tool generates codon-aligned alignment files in IPD-IMGT format:

```
CODON-ALIGNED ANNOTATED ALIGNMENT vs Reference: HLA00434

Ruler           123 456 789 012 345 678 901 234 567 890 ...
AA Codon        -25                 -20                 -15 ...
C*07:02:01:01   ATG CGG GTC ATG GCG CCC CGA GCC CTC CTC ...
C*01:02:01      --- --- --- --- --- --- --- A-- --- A-- ...
C*07:01:01:01   --- --- --- --- --- --- --- --- --- --- ...
```

**Annotation Key:**
- `-` = nucleotide matches reference
- `A/T/G/C` = mismatch (shows actual nucleotide)
- `*` = gap/deletion
- `|` = exon/intron boundary

## Default Reference Alleles

When no custom reference is specified, the following default reference alleles are used:

| Locus | Allele ID | Allele Name |
|-------|-----------|-------------|
| HLA-A | HLA00037 | A*03:01:01:01 |
| HLA-B | HLA00132 | B*07:02:01:01 |
| HLA-C | HLA00434 | C*07:02:01:01 |

You can override these defaults by passing the `reference` parameter to `generate_hla_alighnment()`. The reference allele must exist in the Allelelist.390.txt database.

## Dependencies

| Package | Purpose |
|---------|---------|
| httr2 | HTTP client for API requests |
| readr | File I/O |
| jsonlite | JSON parsing |
| DECIPHER | Multiple sequence alignment |
| Biostrings | Biological sequence handling |
| tidyverse | Data manipulation |

## Data Source

Allele data is sourced from IPD-IMGT/HLA Database version 3.9.0:
- API: https://www.ebi.ac.uk/ipd/imgt/hla/
- GitHub: https://github.com/ANHIG/IMGTHLA

## Running Tests

```R
source("tests/unit_test.R")
```

Tests generate alignments for HLA-A, B, and C loci and compare against reference files in `tests/comp/`.

## License

This project is maintained by Oncology-Sesh.
