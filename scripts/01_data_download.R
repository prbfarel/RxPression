# ============================================
# Script 01: Load and Process GDSC Data
# Deskripsi: Script untuk mengekstrak dan menyelaraskan data ekspresi gen 
#            dan respons obat dari PharmacoSet GDSC
# ============================================

# Purpose: Load GDSC data and create aligned matrices
# Author: prbfarel
# Date: 2025-12-04
# ============================================

# Load required libraries
library(PharmacoGx)
library(SummarizedExperiment)
library(tidyverse)
library(reshape2)

# Set output directories
dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE,showWarnings = FALSE)


# ==============================
# Step 1: Load GDSC PharmacoSet
# ==============================
cat(" STEP 1: Loading GDSC PharmacoSet...\n")

cat("Step 1: Loading GDSC PharmacoSet...\n")
gdsc <- readRDS("data/raw/GDSC_2020(v2-8.2).rds")
cat("✓ PharmacoSet loaded successfully\n\n")

# ====================================
# Step 2: Extract Gene Expression Data
# ====================================
cat(" STEP 2: Extracting gene expression data...\n")

# Extract RNA-seq data from molecularProfiles
# use "RNA" data contains 17.611 genes × 1,018 samples
rna_se <- gdsc@molecularProfiles[["rna"]]

# Extract expression matrix
expr_raw <- assay(rna_se, "exprs")   # Menggunakan assay "exprs"
cat(sprintf("  - Raw expression matrix: %d genes × %d samples\n", 
            nrow(expr_raw), ncol(expr_raw)))

# Transpose to samples × genes format
expr_matrix <- as.data.frame(t(expr_raw))
cat(sprintf("  - Transposed to: %d samples × %d genes\n", 
            nrow(expr_matrix), ncol(expr_matrix)))

# Extract metadata for CEL file to cell line mapping
rna_meta <- colData(rna_se)
cat(sprintf("  - Metadata extracted: %d entries\n", nrow(rna_meta)))

# Create CEL file ID to cell line ID mapping
# colnames(expr_raw) berisi CEL file IDs (e.g., "5500994157493061613625_A01.cel")
# rna_meta$sampleid berisi cell line names (e.g., "COLO 205")
cell_to_cell <- setNames(rna_meta$sampleid, colnames(expr_raw))

# Map rownames dari CEL IDs ke cell line names
expr_matrix$cell_line <- cell_to_cell[rownames(expr_matrix)]
cat("✓ Gene expression data extracted and mapped (cell lines added as column)\n\n")


# ====================================
# Step 3: Extract Drug Response Data
# ====================================

cat("Step 3: Extracting drug response data...\n")

# Extract treatmentResponse components
sens_info <- gdsc@treatmentResponse$info
sens_profiles <- gdsc@treatmentResponse$profiles

cat(sprintf("  - Treatment info: %d experiments\n", nrow(sens_info)))
cat(sprintf("  - Sensitivity profiles: %d measurements\n", nrow(sens_profiles)))

# Create long-format drug response data
# sampleid: cell line ID
# treatmentid: drug ID
# value: aac_recomputed (Area Above Curve)
drug_long <- data.frame(
  cellid = sens_info$sampleid,
  drugid = sens_info$treatmentid,
  value = sens_profiles[, "aac_recomputed"]
)

cat(sprintf("  - Initial drug response entries: %d\n", nrow(drug_long)))

# Remove missing values
drug_long <- drug_long[!is.na(drug_long$value), ]
cat(sprintf("  - After removing NAs: %d entries\n", nrow(drug_long)))

# Convert long format to wide format (samples × drugs)
# fun.aggregate = mean for duplicate handle
drug_matrix <- dcast(drug_long,
                     cellid ~ drugid,
                     value.var = "value",
                     fun.aggregate = mean)

# Set cell line names as rownames
rownames(drug_matrix) <- drug_matrix$cellid
drug_matrix$cellid <- NULL

cat(sprintf("  - Wide format drug matrix: %d samples × %d drugs\n", 
            nrow(drug_matrix), ncol(drug_matrix)))
cat("✓ Drug response data extracted and reshaped\n\n")

# ====================================
# Step 4: Align Data by Common Samples
# ====================================

cat("Step 4: Aligning expression and drug response data...\n")

# Find common cell lines (using the cell_line column for expr_matrix)
common_cell_lines <- intersect(unique(expr_matrix$cell_line), rownames(drug_matrix))
cat(sprintf("  - Unique cell lines in expression: %d\n", length(unique(expr_matrix$cell_line))))
cat(sprintf("  - Cell lines in drug response: %d\n", nrow(drug_matrix)))
cat(sprintf("  - Common cell lines: %d\n", length(common_cell_lines)))

if(length(common_cell_lines) == 0) {
  stop("ERROR: No common cell lines found between expression and drug response data!")
}

# For expr_matrix, select one sample per common cell line (e.g., the first one)
expr_aligned <- expr_matrix[expr_matrix$cell_line %in% common_cell_lines, ]
expr_aligned <- expr_aligned[!duplicated(expr_aligned$cell_line), ]
rownames(expr_aligned) <- expr_aligned$cell_line
expr_aligned$cell_line <- NULL

# Align drug_matrix to common cell lines
drug_aligned <- drug_matrix[common_cell_lines, ]

# Verify alignment
if(!all(rownames(expr_aligned) == rownames(drug_aligned))) {
  stop("ERROR: Cell line alignment failed!")
}

cat("✓ Data alignment successful\n\n")

# ====================================
# Step 5: Extract Metadata
# ====================================

cat("Step 5: Extracting metadata...\n")

# Extract drug information
drug_info <- as.data.frame(gdsc@treatment)
cat(sprintf("  - Drug metadata: %d drugs\n", nrow(drug_info)))

# Extract cell line information
cell_info <- as.data.frame(gdsc@sample)
cat(sprintf("  - Cell line metadata: %d cell lines\n", nrow(cell_info)))

cat("✓ Metadata extracted\n\n")

# ====================================
# Step 6: Save Processed Data
# ====================================

cat("Step 6: Saving processed data...\n")

# Save aligned expression data
saveRDS(expr_aligned, "data/raw/gene_expression.rds")
cat(sprintf("  ✓ Saved: gene_expression.rds (%d samples × %d genes)\n", 
            nrow(expr_aligned), ncol(expr_aligned)))

# Save aligned drug response data
saveRDS(drug_aligned, "data/raw/drug_response.rds")
cat(sprintf("  ✓ Saved: gene_expression.rds (%d samples × %d genes)\n", 
            nrow(expr_aligned), ncol(expr_aligned)))

# Save drug metadata
saveRDS(drug_info, "data/raw/drug_info.rds")
cat(sprintf("  ✓ Saved: drug_info.rds (%d drugs)\n", nrow(drug_info)))

# Save cell line metadata
saveRDS(cell_info, "data/raw/cell_info.rds")
cat(sprintf("  ✓ Saved: cell_info.rds (%d cell lines)\n", nrow(cell_info)))

# ====================================
# Step 7: Data Quality Summary
# ====================================

cat("\n=============================================================================\n")
cat("Data Quality Summary\n")
cat("=============================================================================\n\n")

# Expression data summary
cat("Gene Expression Data:\n")
cat(sprintf("  - Dimensions: %d samples × %d genes\n", 
            nrow(expr_aligned), ncol(expr_aligned)))
cat(sprintf("  - Missing values: %d (%.2f%%)\n", 
            sum(is.na(expr_aligned)), 
            100 * sum(is.na(expr_aligned)) / prod(dim(expr_aligned))))
cat(sprintf("  - Value range: [%.2f, %.2f]\n", 
            min(expr_aligned, na.rm = TRUE), 
            max(expr_aligned, na.rm = TRUE)))

# Drug response data summary
cat("\nDrug Response Data:\n")
cat(sprintf("  - Dimensions: %d samples × %d drugs\n", 
            nrow(drug_aligned), ncol(drug_aligned)))
cat(sprintf("  - Missing values: %d (%.2f%%)\n", 
            sum(is.na(drug_aligned)), 
            100 * sum(is.na(drug_aligned)) / prod(dim(drug_aligned))))
cat(sprintf("  - Value range: [%.2f, %.2f]\n", 
            min(drug_aligned, na.rm = TRUE), 
            max(drug_aligned, na.rm = TRUE)))

# Sample overlap
cat("\nSample Alignment:\n")
cat(sprintf("  - All sample names match: %s\n", 
            ifelse(all(rownames(expr_aligned) == rownames(drug_aligned)), "YES", "NO")))
cat(sprintf("  - Sample order identical: %s\n", 
            ifelse(identical(rownames(expr_aligned), rownames(drug_aligned)), "YES", "NO")))

cat("\n=============================================================================\n")
cat("Data Extraction Completed Succesfully!\n")
cat("=============================================================================\n\n")

cat("Next steps:\n")
cat("  - Run 02_eda.R for exploratory data analysis\n")
cat("  - Check data distributions and relationships\n")
cat("  - Identify potential issues before preprocessing\n\n")

# =============================================================================
# End Of Script
# =============================================================================