# ============================================
# Script 03: Data Preprocessing
# Deskripsi: Script untuk membersihkan, menormalisasi, dan mempersiapkan data untuk machine learning modeling
# ============================================

# Input: 
#   - gene_expression.rds (802 samples × 17,611 genes)
#   - drug_response.rds (802 samples × 190 drugs)
#   - quality_flags.rds
# Output: 
#   - expr_train.rds, expr_test.rds
#   - drug_train.rds, drug_test.rds
#   - feature_info.rds (selected genes info)
#   - preprocessing_report.txt
# Author: prbfarel
# Date: 2025-12-18
# ============================================

# Load required packages
library(tidyverse)
library(caret)
library(VIM)    # For missing value imputation

# Set random seed for reproducibility
set.seed(123)

# Create output directories
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("models", recursive = TRUE, showWarnings = FALSE)
dir.create("reports", recursive = TRUE, showWarnings = FALSE)

# =================================
# Step 1: Load Data & Quality Flags
# =================================

cat("Step 1: Loading data and quality flags...\n")

expr_aligned <- readRDS("data/raw/gene_expression.rds")
drug_aligned <- readRDS("data/raw/drug_response.rds")
quality_flags <- readRDS("data/processed/quality_flags.rds")

cat(sprintf("  - Gene expression: %d samples × %d genes\n", 
            nrow(expr_aligned), ncol(expr_aligned)))
cat(sprintf("  - Drug response: %d samples × %d drugs\n", 
            nrow(drug_aligned), ncol(drug_aligned)))

# Check quality flags
cat("\nQuality flags detected:\n")
if(!is.null(quality_flags$abnormal_samples)) {
  cat(sprintf("  - %d abnormal samples\n", 
              length(quality_flags$abnormal_samples)))
}
if(!is.null(quality_flags$zero_variance_genes)) {
  cat(sprintf("  - %d zero variance genes\n", 
              length(quality_flags$zero_variance_genes)))
}
if(!is.null(quality_flags$high_missing_drugs)) {
  cat(sprintf("  - %d high-missing drugs (>80%%)\n", 
              length(quality_flags$high_missing_drugs)))
}
cat("\n")

# ==================================
# Step 2: Remove Problematic Samples
# ==================================

cat("Step 2: Removing problematic samples...\n")

initial_samples <- nrow(expr_aligned)

# Remove abnormal samples if any
if(!is.null(quality_flags$abnormal_samples)) {
  abnormal_idx <- which(rownames(expr_aligned) %in% quality_flags$abnormal_samples)
  
  if(length(abnormal_idx) > 0) {
    cat(sprintf("  Removing %d abnormal samples:\n", length(abnormal_idx)))
    for(sample_name in quality_flags$abnormal_samples) {
      cat(sprintf("    - %s\n", sample_name))
    }
    
    expr_aligned <- expr_aligned[-abnormal_idx, ]
    drug_aligned <- drug_aligned[-abnormal_idx, ]
  }
}

# Verify alignment after removal
if(!identical(rownames(expr_aligned), rownames(drug_aligned))) {
  stop("ERROR: Sample alignment broken after removal!")
}

cat(sprintf("  ✓ Removed: %d samples\n", 
            initial_samples - nrow(expr_aligned)))
cat(sprintf("  ✓ Remaining: %d samples\n\n", nrow(expr_aligned)))

# ==================================
# Step 3: Remove Low-Variance Genes
# ==================================

cat("Step 3: Removing low-variance genes...\n")

initial_genes <- ncol(expr_aligned)

# Calculate variance for all genes
gene_vars <- apply(expr_aligned, 2, var, na.rm = TRUE)

# Remove bottom 10% lowest variance genes
variance_threshold <- quantile(gene_vars, 0.1)
low_var_genes <- gene_vars < variance_threshold

cat(sprintf("  Variance threshold (10th percentile): %.4f\n", 
            variance_threshold))
cat(sprintf("  Low-variance genes to remove: %d\n", sum(low_var_genes)))

# Remove low-variance genes
expr_aligned <- expr_aligned[, !low_var_genes]

# Also remove zero variance genes if any
if(!is.null(quality_flags$zero_variance_genes)) {
  zero_var_in_remaining <- quality_flags$zero_variance_genes[
    quality_flags$zero_variance_genes %in% colnames(expr_aligned)
  ]
  
  if(length(zero_var_in_remaining) > 0) {
    expr_aligned <- expr_aligned[, 
                                 !colnames(expr_aligned) %in% zero_var_in_remaining]
    cat(sprintf("  Additional zero-variance genes removed: %d\n", 
                length(zero_var_in_remaining)))
  }
}

cat(sprintf("  ✓ Removed: %d genes\n", 
            initial_genes - ncol(expr_aligned)))
cat(sprintf("  ✓ Remaining: %d genes\n\n", ncol(expr_aligned)))

# ==================================
# Step 4: Remove High-Missing Drugs
# ==================================

cat("Step 4: Removing drugs with >80% missing values...\n")

initial_drugs <- ncol(drug_aligned)

# Calculate missing percentage per drug
drug_missing_pct <- 100 * colSums(is.na(drug_aligned)) / nrow(drug_aligned)

# Identify high-missing drugs (>80%)
high_missing_drugs <- names(drug_missing_pct)[drug_missing_pct > 80]

if(length(high_missing_drugs) > 0) {
  cat(sprintf("  Removing %d drugs with >80%% missing:\n", 
              length(high_missing_drugs)))
  
  # Show first 10 if many
  to_show <- min(10, length(high_missing_drugs))
  for(i in 1:to_show) {
    cat(sprintf("    - %s (%.1f%% missing)\n", 
                high_missing_drugs[i], 
                drug_missing_pct[high_missing_drugs[i]]))
  }
  if(length(high_missing_drugs) > 10) {
    cat(sprintf("    ... and %d more\n", length(high_missing_drugs) - 10))
  }
  
  # Remove high-missing drugs
  drug_aligned <- drug_aligned[, !colnames(drug_aligned) %in% high_missing_drugs]
}

cat(sprintf("  ✓ Removed: %d drugs\n", 
            initial_drugs - ncol(drug_aligned)))
cat(sprintf("  ✓ Remaining: %d drugs\n\n", ncol(drug_aligned)))


# ==============================================
# Step 5: Handle Missing Values in Drug Response
# ==============================================

cat("Step 5: Handling missing values in drug response...\n")

# Calculate missing stats
total_values <- prod(dim(drug_aligned))
missing_values <- sum(is.na(drug_aligned))
missing_pct <- 100 * missing_values / total_values

cat(sprintf("  Current missing: %d values (%.2f%%)\n", 
            missing_values, missing_pct))

# Strategy: KNN imputation for drugs with <50% missing
drug_missing_pct_updated <- 100 * colSums(is.na(drug_aligned)) / nrow(drug_aligned)

drugs_to_impute <- names(drug_missing_pct_updated)[drug_missing_pct_updated > 0 & 
                                                     drug_missing_pct_updated < 50]
drugs_complete <- names(drug_missing_pct_updated)[drug_missing_pct_updated == 0]

cat(sprintf("  - Drugs with 0%% missing: %d (no imputation needed)\n", 
            length(drugs_complete)))
cat(sprintf("  - Drugs with 0-50%% missing: %d (will use KNN imputation)\n", 
            length(drugs_to_impute)))

# KNN Imputation
if(length(drugs_to_impute) > 0) {
  cat("\n  Performing KNN imputation (k=5)...\n")
  
  # KNN imputation using VIM package
  # Note: VIM::kNN works on data frames
  drug_aligned_df <- as.data.frame(drug_aligned)
  
  # Impute (this may take a few minutes)
  drug_imputed <- kNN(drug_aligned_df, k = 5, imp_var = FALSE)
  
  # Convert back to matrix format
  drug_aligned <- as.data.frame(drug_imputed)
  rownames(drug_aligned) <- rownames(drug_aligned_df)
  
  # Verify imputation
  remaining_missing <- sum(is.na(drug_aligned))
  cat(sprintf("  ✓ Imputation complete\n"))
  cat(sprintf("  ✓ Remaining missing: %d values (%.2f%%)\n\n", 
              remaining_missing, 
              100 * remaining_missing / prod(dim(drug_aligned))))
} else {
  cat("  ✓ No imputation needed (all drugs either complete or already removed)\n\n")
}

# ==========================================
# Step 6: Feature Selection (Variance-Based)
# ==========================================

cat("Step 6: Feature selection (variance-based filtering)...\n")

# Recalculate variance after cleaning
gene_vars_clean <- apply(expr_aligned, 2, var, na.rm = TRUE)

# Select top N most variable genes
n_features_to_keep <- 5000   # Keep top 5000 genes

if(ncol(expr_aligned) > n_features_to_keep) {
  
  top_var_genes <- names(sort(gene_vars_clean, decreasing = TRUE)[1:n_features_to_keep])
  
  cat(sprintf("  Selecting top %d most variable genes\n", n_features_to_keep))
  cat(sprintf("  Variance range of selected genes: [%.4f, %.4f]\n", 
              min(gene_vars_clean[top_var_genes]), 
              max(gene_vars_clean[top_var_genes])))
  
  # Keep only selected genes
  expr_aligned <- expr_aligned[, top_var_genes]
  
  # Save variance info
  gene_variance_info <- data.frame(
    gene = names(gene_vars_clean),
    variance = gene_vars_clean,
    selected = names(gene_vars_clean) %in% top_var_genes
  ) %>% arrange(desc(variance))
  
  saveRDS(gene_variance_info, "data/processed/gene_variance_info.rds")
  
  cat(sprintf("  ✓ Selected: %d genes\n", ncol(expr_aligned)))
  cat(sprintf("  ✓ Gene variance info saved\n\n"))
  
} else {
  cat(sprintf("  Current genes (%d) <= target (%d), keeping all\n\n", 
              ncol(expr_aligned), n_features_to_keep))
}

# ===============================
# Step 7: Normalization & Scaling
# ===============================

cat("Step 7: Normalization and scaling...\n")

# 7.1 Gene Expression: Z-score normalization per sample
cat("  7.1 Gene expression normalization:\n")
cat("      - Method: Z-score (per sample)\n")

expr_normalized <- t(scale(t(expr_aligned)))   # Scale per sample (row-wise)

# Check for any NA introduced by scaling (constant samples)
if(any(is.na(expr_normalized))) {
  cat("      WARNING: NAs introduced by scaling, replacing with 0\n")
  expr_normalized[is.na(expr_normalized)] <- 0
}


cat(sprintf("      - Before: mean=%.3f, sd=%.3f\n", 
            mean(expr_aligned), sd(as.matrix(expr_aligned))))
cat(sprintf("      - After:  mean=%.3f, sd=%.3f\n", 
            mean(expr_normalized), sd(expr_normalized)))

# Drug Response: Already in [0,1] scale, keep as it is
cat("\n  7.2 Drug response normalization:\n")
cat("      - AAC already in [0,1] scale\n")
cat("      - No additional scaling needed\n")
cat(sprintf("      - Range: [%.3f, %.3f]\n\n", 
            min(drug_aligned, na.rm = TRUE), 
            max(drug_aligned, na.rm = TRUE)))

# Convert to data frame
expr_normalized <- as.data.frame(expr_normalized)
drug_aligned <- as.data.frame(drug_aligned)

# ===============================
# Step 8: Train/Test split
# ===============================

cat("Step 8: Creating train/test split (80/20)...\n")

# Create stratified split id possible (by tissue type)
# For now, random split
n_samples <- nrow(expr_normalized)
n_train <- floor(0.8 * n_samples)
n_test <- n_samples - n_train

# Random sampling for train/test
train_idx <- sample(1:n_samples, n_train)
test_idx <- setdiff(1:n_samples, train_idx)

# Split expression data
expr_train <- expr_normalized[train_idx, ]
expr_test <- expr_normalized[test_idx, ]

# Split drug response data
drug_train <- drug_aligned[train_idx, ]
drug_test <- drug_aligned[test_idx, ]

cat(sprintf("  Train set: %d samples (%.1f%%)\n", 
            nrow(expr_train), 100 * nrow(expr_train) / n_samples))
cat(sprintf("  Test set:  %d samples (%.1f%%)\n", 
            nrow(expr_test), 100 * nrow(expr_test) / n_samples))

# Verify split
cat("\n  Verification:\n")
cat(sprintf("    - Train + Test = %d samples ✓\n", 
            nrow(expr_train) + nrow(expr_test)))
cat(sprintf("    - No overlap: %s ✓\n", 
            ifelse(length(intersect(rownames(expr_train), 
                                    rownames(expr_test))) == 0, "TRUE", "FALSE")))
cat(sprintf("    - Same genes: %s ✓\n", 
            ifelse(identical(colnames(expr_train), colnames(expr_test)), 
                   "TRUE", "FALSE")))
cat(sprintf("    - Same drugs: %s ✓\n\n", 
            ifelse(identical(colnames(drug_train), colnames(drug_test)), 
                   "TRUE", "FALSE")))

# ===============================
# Step 9: Save Processed Data
# ===============================

cat("Step 9: Saving processed data...\n")

# Save train  data
saveRDS(expr_train, "data/processed/expr_train.rds")
saveRDS(drug_train, "data/processed/drug_train.rds")

cat(sprintf("  ✓ expr_train.rds: %d samples × %d genes\n", 
            nrow(expr_train), ncol(expr_train)))
cat(sprintf("  ✓ drug_train.rds: %d samples × %d drugs\n", 
            nrow(drug_train), ncol(drug_train)))

# Save test data
saveRDS(expr_test, "data/processed/expr_test.rds")
saveRDS(drug_test, "data/processed/drug_test.rds")

cat(sprintf("  ✓ expr_test.rds:  %d samples × %d genes\n", 
            nrow(expr_test), ncol(expr_test)))
cat(sprintf("  ✓ drug_test.rds:  %d samples × %d drugs\n\n", 
            nrow(expr_test), ncol(drug_test)))

# Save train/test indices
split_info <- list(
  train_idx = train_idx,
  test_idx = test_idx,
  train_samples = rownames(expr_train),
  test_samples = rownames(expr_test)
)
saveRDS(split_info, "data/processed/train_test_split.rds")
cat("  ✓ train_test_split.rds saved\n\n")

# ======================================
# Step 10: Generate Preprocessing Report
# ======================================

cat("Step 10: Generating preprocessing report...\n")

report <- c(
  "=============================================================================",
  "GDSC Data Preprocessing Report",
  "=============================================================================",
  "",
  paste("Date:", Sys.time()),
  paste("Random seed:", 123),
  "",
  "Input Data:",
  sprintf("  - Original samples: %d", 802),
  sprintf("  - Original genes: %d", 17611),
  sprintf("  - Original drugs: %d", 190),
  "",
  "Preprocessing Steps:",
  "",
  "1. Sample Cleaning:",
  sprintf("   - Abnormal samples removed: %d", 
          ifelse(is.null(quality_flags$abnormal_samples), 0, 
                 length(quality_flags$abnormal_samples))),
  sprintf("   - Remaining samples: %d", n_samples),
  "",
  "2. Gene Filtering:",
  sprintf("   - Low-variance genes removed: %d", 17611 - ncol(expr_aligned)),
  sprintf("   - Variance threshold: %.4f", variance_threshold),
  sprintf("   - Selected top variable genes: %d", ncol(expr_aligned)),
  sprintf("   - Final genes: %d", ncol(expr_aligned)),
  "",
  "3. Drug Filtering:",
  sprintf("   - High-missing drugs removed (>80%%): %d", 
          ifelse(is.null(high_missing_drugs), 0, length(high_missing_drugs))),
  sprintf("   - Final drugs: %d", ncol(drug_aligned)),
  "",
  "4. Missing Value Imputation:",
  sprintf("   - Method: KNN (k=5)"),
  sprintf("   - Drugs imputed: %d", length(drugs_to_impute)),
  sprintf("   - Final missing rate: %.2f%%", 
          100 * sum(is.na(drug_aligned)) / prod(dim(drug_aligned))),
  "",
  "5. Normalization:",
  sprintf("   - Expression: Z-score per sample"),
  sprintf("   - Drug response: No additional scaling (already 0-1)"),
  "",
  "6. Train/Test Split:",
  sprintf("   - Train samples: %d (%.1f%%)", 
          n_train, 100 * n_train / n_samples),
  sprintf("   - Test samples: %d (%.1f%%)", 
          n_test, 100 * n_test / n_samples),
  sprintf("   - Split method: Random sampling"),
  "",
  "Output Files:",
  sprintf("   - expr_train.rds: %d × %d", nrow(expr_train), ncol(expr_train)),
  sprintf("   - drug_train.rds: %d × %d", nrow(drug_train), ncol(drug_train)),
  sprintf("   - expr_test.rds: %d × %d", nrow(expr_test), ncol(expr_test)),
  sprintf("   - drug_test.rds: %d × %d", nrow(drug_test), ncol(drug_test)),
  sprintf("   - train_test_split.rds: split indices"),
  sprintf("   - gene_variance_info.rds: gene selection info"),
  "",
  "Data Characteristics:",
  "",
  "Expression (train):",
  sprintf("   - Mean: %.3f", mean(as.matrix(expr_train))),
  sprintf("   - SD: %.3f", sd(as.matrix(expr_train))),
  sprintf("   - Range: [%.3f, %.3f]", 
          min(expr_train), max(expr_train)),
  "",
  "Drug Response (train):",
  sprintf("   - Mean AAC: %.3f", mean(as.matrix(drug_train), na.rm = TRUE)),
  sprintf("   - Median AAC: %.3f", median(as.matrix(drug_train), na.rm = TRUE)),
  sprintf("   - Range: [%.3f, %.3f]", 
          min(drug_train, na.rm = TRUE), max(drug_train, na.rm = TRUE)),
  sprintf("   - Missing: %.2f%%", 
          100 * sum(is.na(drug_train)) / prod(dim(drug_train))),
  "",
  "=============================================================================",
  "Preprocessing Completed Succesfully",
  "=============================================================================",
  "",
  "Next Steps:",
  "  1. Review preprocessing report",
  "  2. Check data distributions",
  "  3. Run 04_feature_selection.R (if needed)",
  "  4. Run 05_train.R for model training",
  ""
)

# Save report
writeLines(report, "reports/preprocessing_report.txt")
cat("  ✓ Preprocessing report saved to reports/preprocessing_report.txt\n\n")

# Print report summary
cat("=============================================================================\n")
cat("Preprocessing Summary\n")
cat("=============================================================================\n\n")

cat("Before → After:\n")
cat(sprintf("  Samples: %d → %d (removed %d)\n", 
            802, n_samples, 802 - n_samples))
cat(sprintf("  Genes:   %d → %d (removed %d)\n", 
            17611, ncol(expr_aligned), 17611 - ncol(expr_aligned)))
cat(sprintf("  Drugs:   %d → %d (removed %d)\n", 
            190, ncol(drug_aligned), 190 - ncol(drug_aligned)))

cat("\nTrain/Test Split:\n")
cat(sprintf("  Train: %d samples × %d genes × %d drugs\n", 
            nrow(expr_train), ncol(expr_train), ncol(drug_train)))
cat(sprintf("  Test:  %d samples × %d genes × %d drugs\n", 
            nrow(expr_test), ncol(expr_test), ncol(drug_test)))

cat("\nData Quality:\n")
cat(sprintf("  Expression normalized: YES (Z-score)\n"))
cat(sprintf("  Missing values handled: YES (KNN imputation)\n"))
cat(sprintf("  Low-variance genes removed: YES\n"))
cat(sprintf("  High-missing drugs removed: YES\n"))

cat("\n=============================================================================\n")
cat("PREPROCESSING COMPLETED SUCCESSFULLY!\n")
cat("=============================================================================\n\n")

cat("Generated files (in data/processed/):\n")
cat("  1. expr_train.rds - Training expression data\n")
cat("  2. drug_train.rds - Training drug response data\n")
cat("  3. expr_test.rds - Test expression data\n")
cat("  4. drug_test.rds - Test drug response data\n")
cat("  5. train_test_split.rds - Split information\n")
cat("  6. gene_variance_info.rds - Gene selection details\n")
cat("\nReport:\n")
cat("  - reports/preprocessing_report.txt\n\n")

cat("Next steps:\n")
cat("  - Review preprocessing_report.txt\n")
cat("  - Optionally: Run additional feature selection\n")
cat("  - Run 05_train.R to start model training\n\n")

# =============================================================================
# End of Script
# =============================================================================