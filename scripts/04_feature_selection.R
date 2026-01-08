# ============================================
# Script 04: Advance Feature Selection
# Deskripsi: Script untuk melakukan feature selection lanjutan menggunakan
#            multiple methods untuk mengidentifikasi genes paling predictive
# ============================================

# Input: 
#   - expr_train.rds (640 samples × 5,000 genes)
#   - drug_train.rds (640 samples × 173 drugs)
# Output: 
#   - feature_selection_results.rds (scores per drug)
#   - selected_features_per_drug.rds (top genes per drug)
#   - feature_selection_report.txt
#   - Feature importance plots
# Author: prbfarel
# Date: 2025-12-27
# ============================================

# Load required packages
library(tidyverse)
library(caret)
library(glmnet)      # For LASSO
library(randomForest) # For RF feature importance
library(infotheo)    # For mutual information
library(parallel)    # For parallel processing
library(ggplot2)

# Set random seed
set.seed(123)

# Create output directories
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/feature_selection", recursive = TRUE, showWarnings = FALSE)
dir.create("reports", recursive = TRUE, showWarnings = FALSE)

# =================================
# Step 1: Load Preprocessed Data
# =================================

cat("Step 1: Loading preprocessed data...\n")

expr_train <- readRDS("data/processed/expr_train.rds")
drug_train <- readRDS("data/processed/drug_train.rds")

cat(sprintf("  - Expression: %d samples × %d genes\n", 
            nrow(expr_train), ncol(expr_train)))
cat(sprintf("  - Drug response: %d samples × %d drugs\n\n", 
            nrow(drug_train), ncol(drug_train)))

# Verify alignment
if(!identical(rownames(expr_train), rownames(drug_train))) {
  stop("ERROR: Sample names do not match between expression and drug data!")
}

# ======================================================
# Step 2: Correlation-Based Feature Selection (per Drug)
# ======================================================

cat("Step 2: Computing gene-drug correlations...\n")

# Function to compute correlation for one drug
compute_correlations <- function(drug_name, expr_data, drug_data) {
  drug_response <- drug_data[, drug_name]
  
  # Remove samples with NA (if any remaining)
  valid_idx <- !is.na(drug_response)
  
  if(sum(valid_idx) < 10) {
    # Not enough samples
    return(NULL)
  }
  
  # Compute correlation for each gene
  cors <- sapply(colnames(expr_data), function(gene) {
    cor(expr_data[valid_idx, gene], 
        drug_response[valid_idx], 
        method = "pearson")
  })
  
  # Return as data frame
  data.frame(
    drug = drug_name,
    gene = names(cors),
    correlation = cors,
    abs_correlation = abs(cors)
  )
}

# Compute for all drugs
cat("  Computing correlations for all drugs...\n")
correlation_results <- lapply(colnames(drug_train), function(drug) {
  compute_correlations(drug, expr_train, drug_train)
})

# Remove NULL results
correlation_results <- correlation_results[!sapply(correlation_results, is.null)]

# Combine into single data frame
correlation_df <- bind_rows(correlation_results)

cat(sprintf("  ✓ Computed correlations for %d drugs\n", 
            length(unique(correlation_df$drug))))
cat(sprintf("  ✓ Total correlation values: %s\n\n", 
            format(nrow(correlation_df), big.mark = ",")))

# ==================================================
# Step 3: Mutual Information-Based Feature Selection
# ==================================================

cat("Step 3: Computing mutual information scores...\n")

# Function to compute MI for one drug
compute_mi <- function(drug_name, expr_data, drug_data, n_bins = 10) {
  drug_response <- drug_data[, drug_name]
  
  valid_idx <- !is.na(drug_response)
  
  if(sum(valid_idx) < 10) {
    return(NULL)
  }
  
  # Discretize drug response for MI calculation
  drug_discrete <- discretize(drug_response[valid_idx], nbins = n_bins)
  
  # Compute MI for each gene
  mi_scores <- sapply(colnames(expr_data), function(gene) {
    gene_discrete <- discretize(expr_data[valid_idx, gene], nbins = n_bins)
    mutinformation(gene_discrete, drug_discrete)
  })
  
  data.frame(
    drug = drug_name,
    gene = names(mi_scores),
    mi_score = mi_scores
  )
}

# Compute for first 20 drugs as example (MI is slow)
cat("  Computing MI for subset of drugs (this may take a while)...\n")
example_drugs <- colnames(drug_train)[1:20]

mi_results <- lapply(example_drugs, function(drug) {
  if((which(example_drugs == drug) %% 5) == 0) {
    cat(sprintf("    Processing drug %d/%d\n", 
                which(example_drugs == drug), length(example_drugs)))
  }
  compute_mi(drug, expr_train, drug_train)
})

mi_results <- mi_results[!sapply(mi_results, is.null)]
mi_df <- bind_rows(mi_results)

cat(sprintf("  ✓ Computed MI for %d drugs\n\n", 
            length(unique(mi_df$drug))))

# =====================================
# Step 4: Lasso-Based Feature Selection
# =====================================

cat("Step 4: LASSO feature selection...\n")

# Function to run LASSO for one drug
lasso_selection <- function(drug_name, expr_data, drug_data) {
  drug_response <- drug_data[, drug_name]
  
  valid_idx <- !is.na(drug_response)
  
  if(sum(valid_idx) < 20) {
    return(NULL)
  }
  
  # Prepare data
  X <- as.matrix(expr_data[valid_idx, ])
  y <- drug_response[valid_idx]
  
  # Run cross-validated LASSO
  cv_lasso <- tryCatch({
    cv.glmnet(X, y, alpha = 1, nfolds = 5, standardize = TRUE)
  }, error = function(e) NULL)
  
  if(is.null(cv_lasso)) {
    return(NULL)
  }
  
  # Get coefficients at lambda.min
  coefs <- coef(cv_lasso, s = "lambda.min")
  coefs <- coefs[-1, ]  # Remove intercept
  
  # Extract non-zero coefficients
  nonzero_idx <- which(coefs != 0)
  
  if(length(nonzero_idx) == 0) {
    return(NULL)
  }
  
  data.frame(
    drug = drug_name,
    gene = names(coefs)[nonzero_idx],
    lasso_coef = coefs[nonzero_idx],
    abs_lasso_coef = abs(coefs[nonzero_idx])
  )
}

# Run LASSO for first 10 drugs as example
cat("  Running LASSO for subset of drugs...\n")
lasso_drugs <- colnames(drug_train)[1:10]

lasso_results <- lapply(lasso_drugs, function(drug) {
  if((which(lasso_drugs == drug) %% 3) == 0) {
    cat(sprintf("    Processing drug %d/%d\n", 
                which(lasso_drugs == drug), length(lasso_drugs)))
  }
  lasso_selection(drug, expr_train, drug_train)
})

lasso_results <- lasso_results[!sapply(lasso_results, is.null)]
lasso_df <- bind_rows(lasso_results)

cat(sprintf("  ✓ LASSO completed for %d drugs\n", 
            length(unique(lasso_df$drug))))
cat(sprintf("  ✓ Average features selected: %.1f\n\n", 
            mean(table(lasso_df$drug))))

# =======================================
# Step 5: Random Forest Feature Selection
# =======================================

cat("Step 5: Random Forest feature importance...\n")

# Function to compute RF importance for one drug
rf_importance <- function(drug_name, expr_data, drug_data, n_trees = 100) {
  drug_response <- drug_data[, drug_name]
  
  valid_idx <- !is.na(drug_response)
  
  if(sum(valid_idx) < 20) {
    return(NULL)
  }
  
  # Prepare data (sample features to speed up)
  set.seed(123)
  sample_genes <- sample(colnames(expr_data), min(500, ncol(expr_data)))
  X <- expr_data[valid_idx, sample_genes]
  y <- drug_response[valid_idx]
  
  # Train RF
  rf_model <- tryCatch({
    randomForest(x = X, y = y, ntree = n_trees, importance = TRUE)
  }, error = function(e) NULL)
  
  if(is.null(rf_model)) {
    return(NULL)
  }
  
  # Extract importance
  importance_scores <- importance(rf_model)[, "%IncMSE"]
  
  data.frame(
    drug = drug_name,
    gene = names(importance_scores),
    rf_importance = importance_scores
  )
}

# Run RF for first 5 drugs as example (RF is very slow)
cat("  Training Random Forest for subset of drugs...\n")
rf_drugs <- colnames(drug_train)[1:5]

rf_results <- lapply(rf_drugs, function(drug) {
  cat(sprintf("    Processing drug %d/%d (%s)\n", 
              which(rf_drugs == drug), length(rf_drugs), drug))
  rf_importance(drug, expr_train, drug_train)
})

rf_results <- rf_results[!sapply(rf_results, is.null)]
rf_df <- bind_rows(rf_results)

cat(sprintf("  ✓ RF completed for %d drugs\n\n", 
            length(unique(rf_df$drug))))

# =========================================
# Step 6: Combine Feature Selection Methods
# =========================================

cat("Step 6: Combining feature selection methods...\n")

# Select top N features from each method for each drug
top_n <- 100  # Top 100 genes per method

# Function to get top features for one drug
combine_methods <- function(drug_name) {
  
  # Top by correlation
  top_cor <- correlation_df %>%
    filter(drug == drug_name) %>%
    arrange(desc(abs_correlation)) %>%
    head(top_n) %>%
    pull(gene)
  
  # Top by MI (if available)
  top_mi <- if(drug_name %in% mi_df$drug) {
    mi_df %>%
      filter(drug == drug_name) %>%
      arrange(desc(mi_score)) %>%
      head(top_n) %>%
      pull(gene)
  } else {
    character(0)
  }
  
  # LASSO selected (if available)
  lasso_genes <- if(drug_name %in% lasso_df$drug) {
    lasso_df %>%
      filter(drug == drug_name) %>%
      pull(gene)
  } else {
    character(0)
  }
  
  # RF importance (if available)
  top_rf <- if(drug_name %in% rf_df$drug) {
    rf_df %>%
      filter(drug == drug_name) %>%
      arrange(desc(rf_importance)) %>%
      head(top_n) %>%
      pull(gene)
  } else {
    character(0)
  }
  
  # Combine: genes that appear in multiple methods get higher scores
  all_genes <- c(top_cor, top_mi, lasso_genes, top_rf)
  gene_counts <- table(all_genes)
  
  # Create ensemble score
  ensemble_genes <- data.frame(
    gene = names(gene_counts),
    vote_count = as.numeric(gene_counts),
    in_correlation = names(gene_counts) %in% top_cor,
    in_mi = names(gene_counts) %in% top_mi,
    in_lasso = names(gene_counts) %in% lasso_genes,
    in_rf = names(gene_counts) %in% top_rf
  ) %>%
    arrange(desc(vote_count))
  
  # Get correlation values for ranking
  cor_values <- correlation_df %>%
    filter(drug == drug_name, gene %in% ensemble_genes$gene) %>%
    select(gene, abs_correlation)
  
  ensemble_genes <- ensemble_genes %>%
    left_join(cor_values, by = "gene") %>%
    arrange(desc(vote_count), desc(abs_correlation))
  
  list(
    drug = drug_name,
    top_features = head(ensemble_genes, 200),  # Top 200 genes
    all_features = ensemble_genes
  )
}

# Combine for all drugs
cat("  Combining methods for all drugs...\n")
drugs_with_results <- unique(correlation_df$drug)

combined_results <- lapply(drugs_with_results, combine_methods)
names(combined_results) <- drugs_with_results

cat(sprintf("  ✓ Combined results for %d drugs\n\n", 
            length(combined_results)))

# =================================
# Step 7: Create Final Feature Sets
# =================================

cat("Step 7: Creating final feature sets per drug...\n")

# Extract top 200 genes per drug
selected_features_per_drug <- lapply(combined_results, function(x) {
  x$top_features$gene
})

# Count how many drugs each gene is selected for
gene_selection_counts <- table(unlist(selected_features_per_drug))
gene_selection_counts <- sort(gene_selection_counts, decreasing = TRUE)

cat(sprintf("  Top 10 most frequently selected genes:\n"))
print(head(gene_selection_counts, 10))

cat(sprintf("\n  Average genes selected per drug: %.1f\n", 
            mean(sapply(selected_features_per_drug, length))))
cat(sprintf("  Median genes selected per drug: %.1f\n\n", 
            median(sapply(selected_features_per_drug, length))))

# ===========================================
# Step 8: Visualize Feature Selection Results
# ===========================================

cat("Step 8: Creating visualizations...\n")

# 8.1 Gene selection frequency
png("figures/feature_selection/01_gene_selection_frequency.png",
    width = 12, height = 6, units = "in", res = 300)

top_genes_df <- data.frame(
  gene = names(head(gene_selection_counts, 30)),
  count = as.numeric(head(gene_selection_counts, 30))
)

ggplot(top_genes_df, aes(x = reorder(gene, count), y = count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Top 30 Most Frequently Selected Genes",
       subtitle = "Number of drugs for which each gene was selected (top 200)",
       x = "Gene", y = "Number of Drugs") +
  theme_bw(base_size = 12)

dev.off()
cat("  ✓ Saved: 01_gene_selection_frequency.png\n")

# 8.2 Feature selection method comparison (for example drugs)
example_drugs_viz <- names(combined_results)[1:min(5, length(combined_results))]

png("figures/feature_selection/02_method_comparison.png",
    width = 14, height = 10, units = "in", res = 300)

par(mfrow = c(3, 2))

for(drug in example_drugs_viz) {
  features <- combined_results[[drug]]$top_features
  
  method_counts <- c(
    Correlation = sum(features$in_correlation),
    MI = sum(features$in_mi),
    LASSO = sum(features$in_lasso),
    RF = sum(features$in_rf)
  )
  
  barplot(method_counts, 
          main = paste("Feature Selection Methods:", drug),
          ylab = "Number of Features",
          col = c("skyblue", "lightgreen", "coral", "plum"),
          las = 2,
          cex.main = 0.9)
}

dev.off()
cat("  ✓ Saved: 02_method_comparison.png\n")

# 8.3 Correlation distribution for example drug
example_drug <- names(combined_results)[1]
example_cor <- correlation_df %>% filter(drug == example_drug)
example_selected <- combined_results[[example_drug]]$top_features$gene

png("figures/feature_selection/03_correlation_distribution.png",
    width = 12, height = 6, units = "in", res = 300)

ggplot(example_cor, aes(x = correlation)) +
  geom_histogram(aes(fill = gene %in% example_selected), 
                 bins = 50, alpha = 0.7) +
  scale_fill_manual(values = c("gray70", "steelblue"),
                    labels = c("Not Selected", "Selected (Top 200)"),
                    name = "") +
  labs(title = paste("Gene-Drug Correlation Distribution:", example_drug),
       x = "Pearson Correlation",
       y = "Number of Genes") +
  theme_bw(base_size = 12) +
  theme(legend.position = "top")

dev.off()
cat("  ✓ Saved: 03_correlation_distribution.png\n\n")

# =====================
# Step 9: Save Results
# =====================

cat("Step 9: Saving feature selection results...\n")

# Save all results
saveRDS(correlation_df, "data/processed/correlation_results.rds")
saveRDS(mi_df, "data/processed/mi_results.rds")
saveRDS(lasso_df, "data/processed/lasso_results.rds")
saveRDS(rf_df, "data/processed/rf_results.rds")
saveRDS(combined_results, "data/processed/feature_selection_results.rds")
saveRDS(selected_features_per_drug, "data/processed/selected_features_per_drug.rds")

cat("  ✓ correlation_results.rds saved\n")
cat("  ✓ mi_results.rds saved\n")
cat("  ✓ lasso_results.rds saved\n")
cat("  ✓ rf_results.rds saved\n")
cat("  ✓ feature_selection_results.rds saved\n")
cat("  ✓ selected_features_per_drug.rds saved\n\n")

# ========================
# Step 10: Generate Report
# ========================

cat("Step 10: Generating feature selection report...\n")

report <- c(
  "=============================================================================",
  "GDSC Feature Selection Report",
  "=============================================================================",
  "",
  paste("Date:", Sys.time()),
  paste("Random seed:", 123),
  "",
  "Input Data:",
  sprintf("  - Samples: %d", nrow(expr_train)),
  sprintf("  - Initial genes: %d", ncol(expr_train)),
  sprintf("  - Drugs: %d", ncol(drug_train)),
  "",
  "Feature Selection Mehtods:",
  "",
  "1. Correlation Analysis:",
  sprintf("   - Drugs analyzed: %d", length(unique(correlation_df$drug))),
  sprintf("   - Correlation values computed: %s", 
          format(nrow(correlation_df), big.mark = ",")),
  sprintf("   - Method: Pearson correlation"),
  "",
  "2. Mutual Information:",
  sprintf("   - Drugs analyzed: %d", length(unique(mi_df$drug))),
  sprintf("   - Method: Discretization + MI"),
  sprintf("   - Bins: 10"),
  "",
  "3. LASSO Regression:",
  sprintf("   - Drugs analyzed: %d", length(unique(lasso_df$drug))),
  sprintf("   - Method: Cross-validated LASSO (alpha=1)"),
  sprintf("   - Average features per drug: %.1f", mean(table(lasso_df$drug))),
  "",
  "4. Random Forest:",
  sprintf("   - Drugs analyzed: %d", length(unique(rf_df$drug))),
  sprintf("   - Method: RF with 100 trees"),
  sprintf("   - Importance metric: %%IncMSE"),
  "",
  "Ensemble Feature Selection:",
  "",
  sprintf("  - Top N per method: %d", top_n),
  sprintf("  - Final features per drug: 200"),
  sprintf("  - Drugs with feature sets: %d", length(selected_features_per_drug)),
  "",
  "Selected Features Summary:",
  "",
  sprintf("  - Total unique genes selected: %d", 
          length(unique(unlist(selected_features_per_drug)))),
  sprintf("  - Most frequent gene: %s (in %d drugs)", 
          names(gene_selection_counts)[1], 
          gene_selection_counts[1]),
  sprintf("  - Average genes per drug: %.1f", 
          mean(sapply(selected_features_per_drug, length))),
  sprintf("  - Median genes per drug: %.1f", 
          median(sapply(selected_features_per_drug, length))),
  "",
  "Top 10 Most Frequently Selected Genes:",
  "",
  paste(sprintf("  %d. %s: %d drugs", 
                1:10, 
                names(head(gene_selection_counts, 10)),
                head(gene_selection_counts, 10)), 
        collapse = "\n"),
  "",
  "",
  "Output FIles:",
  "",
  "  Data files (in data/processed/):",
  "    - correlation_results.rds",
  "    - mi_results.rds",
  "    - lasso_results.rds",
  "    - rf_results.rds",
  "    - feature_selection_results.rds",
  "    - selected_features_per_drug.rds",
  "",
  "  Figures (in figures/feature_selection/):",
  "    - 01_gene_selection_frequency.png",
  "    - 02_method_comparison.png",
  "    - 03_correlation_distribution.png",
  "",
  "=============================================================================",
  "FEATURE SELECTION COMPLETED SUCCESSFULLY",
  "=============================================================================",
  "",
  "Next Steps:",
  "  1. Review feature selection visualizations",
  "  2. Examine top selected genes for biological relevance",
  "  3. Run 05_train.R to train models with selected features",
  ""
)

writeLines(report, "reports/feature_selection_report.txt")
cat("  ✓ Report saved to reports/feature_selection_report.txt\n\n")

# =============================================================================
# Summaries
# =============================================================================

cat("=============================================================================\n")
cat("Feature Selection Summary\n")
cat("=============================================================================\n\n")

cat("Methods Applied:\n")
cat(sprintf("  1. Correlation:  %d drugs\n", length(unique(correlation_df$drug))))
cat(sprintf("  2. Mutual Info:  %d drugs\n", length(unique(mi_df$drug))))
cat(sprintf("  3. LASSO:        %d drugs\n", length(unique(lasso_df$drug))))
cat(sprintf("  4. Random Forest: %d drugs\n", length(unique(rf_df$drug))))

cat("\nFeatures Sets Created:\n")
cat(sprintf("  - Per-drug feature sets: %d drugs\n", 
            length(selected_features_per_drug)))
cat(sprintf("  - Features per drug: 200 genes\n"))
cat(sprintf("  - Total unique genes: %d\n", 
            length(unique(unlist(selected_features_per_drug)))))

cat("\nMost Important Genes (top 5):\n")
for(i in 1:5) {
  cat(sprintf("  %d. %s: selected for %d drugs\n", 
              i, 
              names(gene_selection_counts)[i],
              gene_selection_counts[i]))
}

cat("\n=============================================================================\n")
cat("FEATURE SELECTION COMPLETED SUCCESSFULLY!\n")
cat("=============================================================================\n\n")

cat("Generated files:\n")
cat("  Data:\n")
cat("    - correlation_results.rds\n")
cat("    - mi_results.rds\n")
cat("    - lasso_results.rds\n")
cat("    - rf_results.rds\n")
cat("    - feature_selection_results.rds\n")
cat("    - selected_features_per_drug.rds\n")
cat("\n  Figures:\n")
cat("    - 01_gene_selection_frequency.png\n")
cat("    - 02_method_comparison.png\n")
cat("    - 03_correlation_distribution.png\n")
cat("\n  Report:\n")
cat("    - feature_selection_report.txt\n\n")

cat("Next steps:\n")
cat("  - Review feature selection visualizations\n")
cat("  - Check biological relevance of top genes\n")
cat("  - Run 05_train.R to train models\n\n")

# =============================================================================
# End of Script
# =============================================================================