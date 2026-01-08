# ============================================
# Script 02: Exploratory Data Analysis (EDA)
# Deskripsi: Script untuk melakukan analisis eksplorasi data ekspresi gen dan respons obat dari data GDSC
# ============================================

# Input: 
#   - gene_expression.rds (802 samples × 17,611 genes)
#   - drug_response.rds (802 samples × 190 drugs)
#   - drug_info.rds
#   - cell_info.rds
# Output: 
#   - EDA plots (PNG/PDF)
#   - Summary statistics
#   - Quality report
# Author: prbfarel
# Date: 2025-12-09
# ============================================

# Load required packages
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(corrplot)
library(gridExtra)
library(RColorBrewer)
library(scales)

# Set plot theme
theme_set(theme_bw(base_size = 12))

# Create output directories
dir.create("figures/eda", recursive = TRUE, showWarnings = FALSE)
dir.create("reports", recursive = TRUE, showWarnings = FALSE)

# ==============================
# Step 1: Load Data
# ==============================

cat("Step 1: Loading processed data...\n")

expr_aligned <- readRDS("data/raw/gene_expression.rds")
drug_aligned <- readRDS("data/raw/drug_response.rds")
drug_info <- readRDS("data/raw/drug_info.rds")
cell_info <- readRDS("data/raw/cell_info.rds")

cat(sprintf("  ✓ Gene expression: %d samples × %d genes\n", 
            nrow(expr_aligned), ncol(expr_aligned)))
cat(sprintf("  ✓ Drug response: %d samples × %d drugs\n", 
            nrow(drug_aligned), ncol(drug_aligned)))
cat(sprintf("  ✓ Drug metadata: %d drugs\n", nrow(drug_info)))
cat(sprintf("  ✓ Cell line metadata: %d cell lines\n\n", nrow(cell_info)))


# ==============================
# Step 2: Basic Statistics
# ==============================

cat("Step 2: Computing basic statistics...\n")

expr_stats <- data.frame(
  Metric = c("Mean", "Median", "SD", "Min", "Max",
             "Q1", "Q3", "Missing Values"),
  Value = c(
    mean(as.matrix(expr_aligned), na.rm = TRUE),
    median(as.matrix(expr_aligned), na.rm = TRUE),
    sd(as.matrix(expr_aligned), na.rm = TRUE),
    min(expr_aligned, na.rm = TRUE),
    max(expr_aligned, na.rm = TRUE),
    quantile(as.matrix(expr_aligned), 0.25, na.rm = TRUE),
    quantile(as.matrix(expr_aligned), 0.75, na.rm = TRUE),
    sum(is.na(expr_aligned))
  )
)

cat("\nGene Expression Statistics:\n")
print(expr_stats, row.names = FALSE)

# Drug Response Statistics
drug_stats <- data.frame(
  Metric = c("Mean", "Median", "SD", "Min", "Max", 
             "Q1", "Q3", "Missing Values", "Missing %"),
  Value = c(
    mean(as.matrix(drug_aligned), na.rm = TRUE),
    median(as.matrix(drug_aligned), na.rm = TRUE),
    sd(as.matrix(drug_aligned), na.rm = TRUE),
    min(drug_aligned, na.rm = TRUE),
    max(drug_aligned, na.rm = TRUE),
    quantile(as.matrix(drug_aligned), 0.25, na.rm = TRUE),
    quantile(as.matrix(drug_aligned), 0.75, na.rm = TRUE),
    sum(is.na(drug_aligned)),
    100 * sum(is.na(drug_aligned)) / prod(dim(drug_aligned))
  )
)

cat("\nDrug Response (AAC) Statistics:\n")
print(drug_stats, row.names = FALSE)
cat("\n")

# ================================
# Step 3: Gene Expression Analysis
# ================================

cat("Step 3: Analyzing gene expression distribution...\n")

# 3.1 Overall distribution
expr_vec <- as.vector(as.matrix(expr_aligned))

png("figures/eda/01_expression_distribution.png",
    width = 12, height = 8, units = "in", res = 300)
par(mfrow = c(2,2))

# Histogram
hist(expr_vec, breaks = 100, col = "skyblue", border = "white",
     main = "Gene Expression Distribution",
     xlab = "Expression Level (log2)", ylab = "Frequency")

# Density plot
plot(density(expr_vec, na.rm = TRUE), main = "Expression Density",
     xlab = "Expression Level (log2)", lwd = 2, col = "darkblue")

# Q-Q plot
qqnorm(sample(expr_vec, 10000), main = "Q-Q Plot (Sample 10k points)")
qqline(sample(expr_vec, 10000), col = "red", lwd = 2)

# Boxplot
boxplot(expr_vec, horizontal = TRUE, col = "lightgreen",
        main = "Expression Boxplot", xlab = "Expression Level (log2)")

dev.off()

cat("  ✓ Saved: 01_expression_distribution.png\n")

# 3.2 Per-sample expression
sample_means <- rowMeans(expr_aligned, na.rm = TRUE)
sample_sds <- apply(expr_aligned, 1, sd, na.rm = TRUE)

png("figures/eda/02_sample_quality.png", 
    width = 14, height = 6, units = "in", res = 300)
par(mfrow = c(1, 2))

# Mean expression per sample
plot(sample_means, pch = 20, col = "blue",
     main = "Mean Expression per Sample",
     xlab = "Sample Index", ylab = "Mean Expression",
     ylim = c(min(sample_means) - 0.5, max(sample_means) + 0.5))
abline(h = median(sample_means), col = "red", lwd = 2, lty = 2)

# SD per sample
plot(sample_sds, pch = 20, col = "darkgreen",
     main = "Expression Variability per Sample",
     xlab = "Sample Index", ylab = "Standard Deviation")
abline(h = median(sample_sds), col = "red", lwd = 2, lty = 2)

dev.off()

cat("  ✓ Saved: 02_sample_quality.png\n")

# 3.3 Per-gene expression
gene_means <- colMeans(expr_aligned, na.rm = TRUE)
gene_vars <- apply(expr_aligned, 2, var, na.rm = TRUE)

png("figures/eda/03_gene_statistics.png", 
    width = 14, height = 6, units = "in", res = 300)
par(mfrow = c(1, 2))

# Mean expression per gene
hist(gene_means, breaks = 100, col = "coral", border = "white",
     main = "Mean Expression per Gene",
     xlab = "Mean Expression", ylab = "Number of Genes")

# Variance per gene
hist(log10(gene_vars +1), breaks = 100, col = "purple", border = "white",
     main = "Gene Expression Variance",
     xlab = "log10(Variance + 1", ylab = "Number of Genes")

dev.off()

cat("  ✓ Saved: 03_gene_statistics.png\n")

# 3.4 Low variance genes
low_var_threshold <- quantile(gene_vars, 0.1)
n_low_var <- sum(gene_vars < low_var_threshold)

cat(sprintf("\n  - Low variance genes (bottom 10%%): %d genes\n", n_low_var))
cat(sprintf("  - Variance threshold: %.4f\n", low_var_threshold))

# ================================
# Step 4: Drug Response Analysis
# ================================

# 4.1 Overall AAC distribution
drug_vec <- as.vector(as.matrix(drug_aligned))
drug_vec_clean <- drug_vec[!is.na(drug_vec)]

png("figures/eda/04_drug_response_distribution.png", 
    width = 12, height = 8, units = "in", res = 300)
par(mfrow = c(2, 2))

# Histogram
hist(drug_vec_clean, breaks = 50, col = "lightblue", border = "white",
     main = "Drug Response (AAC) Distribution",
     xlab = "AAC Value", ylab = "Frequency")

# Density plot
plot(density(drug_vec_clean), main = "AAC Density",
     xlab = "AAC Value", lwd = 2, col = "darkblue")
abline(v = median(drug_vec_clean), col = "red", lwd = 2, lty = 2)

# Q-Q plot
qqnorm(sample(drug_vec_clean, 10000), main = "Q-Q Plot")
qqline(sample(drug_vec_clean, 10000), col = "red", lwd = 2)

# Boxplot
boxplot(drug_vec_clean, horizontal = TRUE, col = "lightcoral",
        main = "AAC Boxplot", xlab = "AAC Value")

dev.off()

cat("  ✓ Saved: 04_drug_response_distribution.png\n")

# 4.2 Missing values per drug
drug_na_count <- colSums(is.na(drug_aligned))
drug_na_percent <- 100 * drug_na_count / nrow(drug_aligned)

# Sort by missing percentage
drug_na_df <- data.frame(
  Drug = names(drug_na_percent),
  Missing_Count = drug_na_count,
  Missing_Percent = drug_na_percent
) %>% arrange(desc(Missing_Percent))

# Top 10 drugs with most missing values
cat("\nTop 10 drugs with most missing values:\n")
print(head(drug_na_df, 10), row.names = FALSE)

# Plot missing values
png("figures/eda/05_missing_values.png", 
    width = 14, height = 8, units = "in", res = 300)
par(mfrow = c(2, 1))

# Histogram of missing percentages
hist(drug_na_percent, breaks = 50, col = "orange", border = "white",
     main = "Distribution of Missing Values per Drug",
     xlab = "Missing Values (%)", ylab = "Number of Drugs")
abline(v = median(drug_na_percent), col = "red", lwd = 2, lty = 2)

# Barplot of top 20 drugs with most missing
top20_missing <- head(drug_na_df, 20)
barplot(top20_missing$Missing_Percent, 
        names.arg = top20_missing$Drug,
        las = 2, col = "salmon",
        main = "Top 20 Drugs with Most Missing Values",
        ylab = "Missing Values (%)",
        cex.names = 0.7)
abline(h = 50, col = "red", lwd = 2, lty = 2)

dev.off()

cat("  ✓ Saved: 05_missing_values.png\n")

# 4.3 Drug response per drug
drug_means <- colMeans(drug_aligned, na.rm = TRUE)
drug_sds <- apply(drug_aligned, 2, sd, na.rm = TRUE)

png("figures/eda/06_drug_statistics.png", 
    width = 14, height = 6, units = "in", res = 300)
par(mfrow = c(1, 2))

# Mean AAC per drug
hist(drug_means, breaks = 30, col = "lightgreen", border = "white",
     main = "Mean AAC per Drug",
     xlab = "Mean AAC", ylab = "Number of Drugs")
abline(v = median(drug_means, na.rm = TRUE), col = "red", lwd = 2, lty = 2)

# SD per drug
hist(drug_sds, breaks = 30, col = "plum", border = "white",
     main = "AAC Variability per Drug",
     xlab = "Standard Deviation", ylab = "Number of Drugs")

dev.off()

cat("  ✓ Saved: 06_drug_statistics.png\n")

# ================================
# Step 5: Sample Clustering & PCA
# ================================

cat("\nStep 5: Performing sample clustering and PCA...\n")

# 5.1 PCA on gene expression (using subset of high variance genes)
# Select top 2000 most variable genes for PCA
top_var_genes <- names(sort(gene_vars, decreasing = TRUE)[1:2000])
expr_for_pca <- expr_aligned[, top_var_genes]

# Perform PCA
pca_result <- prcomp(expr_for_pca, scale. = TRUE, center = TRUE)

# Calculate variance explained
var_explained <- (pca_result$sdev^2 / sum(pca_result$sdev^2)) * 100

png("figures/eda/07_pca_analysis.png", 
    width = 14, height = 6, units = "in", res = 300)
par(mfrow = c(1, 2))

# Scree plot
plot(var_explained[1:20], type = "b", pch = 19, col = "blue",
     main = "PCA = Variance Explained",
     xlab = "Principal Componenent",
     ylab = "Variance Explained (%)")

# PC1 vs PC2
plot(pca_result$x[, 1], pca_result$x[, 2],
     pch = 20, col = alpha("darkblue", 0.5),
     main = "PCA - PC1 vs PC2",
     xlab = sprintf("PC1 (%.2f%% variance)", var_explained[1]),
     ylab = sprintf("PC2 (%.2f%% variance)", var_explained[2]))
grid()

dev.off()

cat(sprintf("  ✓ PC1 explains %.2f%% of variance\n", var_explained[1]))
cat(sprintf("  ✓ PC2 explains %.2f%% of variance\n", var_explained[2]))
cat(sprintf("  ✓ Top 10 PCs explain %.2f%% of variance\n", sum(var_explained[1:10])))
cat("  ✓ Saved: 07_pca_analysis.png\n")


# ================================
# Step 6: Correlation Analysis
# ================================

cat("\nStep 6: Analyzing correlations...\n")

# 6.1 Sample correlation (using subset for speed)
set.seed(42)
sample_subset <- sample(1:ncol(expr_aligned), min(500, ncol(expr_aligned)))
expr_subset <- expr_aligned[, sample_subset]

sample_cor <- cor(t(expr_subset))

png("figures/eda/08_sample_correlation.png", 
    width = 10, height = 10, units = "in", res = 300)
pheatmap(sample_cor, 
         color = colorRampPalette(c("blue", "white", "red"))(100),
         show_rownames = FALSE, 
         show_colnames = FALSE,
         main = "Sample-Sample Correlation (500 genes subset)",
         border_color = NA)
dev.off()

cat("  ✓ Saved: 08_sample_correlation.png\n")

# 6.2 Drug-drug correlation
# Remove drugs with >50% missing values for correlation
drugs_for_cor <- drug_aligned[, drug_na_percent < 50]
drug_cor <- cor(drugs_for_cor, use = "pairwise.complete.obs")

png("figures/eda/09_drug_correlation.png", 
    width = 12, height = 12, units = "in", res = 300)
pheatmap(drug_cor, 
         color = colorRampPalette(c("blue", "white", "red"))(100),
         show_rownames = FALSE, 
         show_colnames = FALSE,
         main = sprintf("Drug-Drug Correlation (n=%d drugs)", ncol(drugs_for_cor)),
         border_color = NA)
dev.off()

cat(sprintf("  ✓ Drugs used for correlation: %d (<%% 50 missing)\n", 
            ncol(drugs_for_cor)))
cat("  ✓ Saved: 09_drug_correlation.png\n")

# =========================================
# Step 7: Expression-Response Relationship
# =========================================

cat("\nStep 7: Exploring expression-response relationships...\n")

# Select a few example drugs with low missing values
example_drugs <- names(sort(drug_na_percent)[1:6])

# For each drug, find top correlated genes
cat("\nExample: Top genes correlated with drug response:\n")

for(i in 1:3) {
  drug <- example_drugs[i]
  drug_values <- drug_aligned[, drug]
  
  #Calculate correlation with all genes
  gene_cors <- apply(expr_aligned, 2, function(gene_expr) {
    valid_idx <- !is.na(drug_values)
    cor(gene_expr[valid_idx], drug_values[valid_idx],
        use = "complete.obs")
  })
  
  # Top 5 correlated genes
  top_genes <- names(sort(abs(gene_cors), decreasing = TRUE)[1:5])
  
  cat(sprintf("\n Drug: %s\n", drug))
  for(gene in top_genes) {
    cat(sprintf("    - %s: r = %.3f\n", gene, gene_cors[gene]))
  }
}

# Create scatter plots for first example
drug1 <- example_drugs[1]
drug1_values <- drug_aligned[, drug1]
valid_idx <- !is.na(drug1_values)

gene_cors1 <- apply(expr_aligned, 2, function(gene_expr) {
  cor(gene_expr[valid_idx], drug1_values[valid_idx], use = "complete.obs")
})

top_genes1 <- names(sort(abs(gene_cors1), decreasing = TRUE)[1:6])

png("figures/eda/10_expression_response_examples.png", 
    width = 14, height = 10, units = "in", res = 300)
par(mfrow = c(2, 3))

for(gene in top_genes1) {
  plot(expr_aligned[valid_idx, gene], drug1_values[valid_idx],
       pch = 20, col = alpha("darkblue", 0.5),
       main = sprintf("%s\nr = %.3f", gene, gene_cors1[gene]),
       xlab = "Gene Expression",
       ylab = sprintf("AAC (%s)", drug1),
       cex.main = 0.9)
  abline(lm(drug1_values[valid_idx] ~ expr_aligned[valid_idx, gene]), 
         col = "red", lwd = 2)
}

dev.off()

cat("\n  ✓ Saved: 10_expression_response_examples.png\n")

# ==========================
# Step 8: Data Quality Flags
# ==========================

cat("\nStep 8: Generating data quality flags...\n")

# Identify potential issues
quality_flags <- list()

# Check 1: Samples with abnormal expression
sample_mean_outliers <- which(abs(scale(sample_means)) > 3)
if(length(sample_mean_outliers) > 0) {
  quality_flags$abnormal_samples <- rownames(expr_aligned)[sample_mean_outliers]
  cat(sprintf("  ⚠ Warning: %d samples with abnormal mean expression\n", 
              length(sample_mean_outliers)))
}

# Check 2: Genes with zero variance
zero_var_genes <- which(gene_vars == 0)
if(length(zero_var_genes) > 0) {
  quality_flags$zero_variance_genes <- names(zero_var_genes)
  cat(sprintf("  ⚠ Warning: %d genes with zero variance\n", 
              length(zero_var_genes)))
}

# Check 3: Drugs with >80% missing
high_missing_drugs <- drug_na_df$Drug[drug_na_df$Missing_Percent > 80]
if(length(high_missing_drugs) > 0) {
  quality_flags$high_missing_drugs <- high_missing_drugs
  cat(sprintf("  ⚠ Warning: %d drugs with >80%% missing values\n", 
              length(high_missing_drugs)))
}

# Save quality flags
saveRDS(quality_flags, "data/processed/quality_flags.rds")
cat("  ✓ Saved: quality_flags.rds\n")


# ==========================
# Step 9: Summary Report
# ==========================

cat("\n=============================================================================\n")
cat("EDA Summary Report\n")
cat("=============================================================================\n\n")

cat("Dataset Overview:\n")
cat(sprintf("  - Samples: %d cell lines\n", nrow(expr_aligned)))
cat(sprintf("  - Genes: %d features\n", ncol(expr_aligned)))
cat(sprintf("  - Drugs: %d compounds\n", ncol(drug_aligned)))

cat("\nGene Expression:\n")
cat(sprintf("  - Mean: %.2f (log2 scale)\n", mean(expr_vec, na.rm = TRUE)))
cat(sprintf("  - Range: [%.2f, %.2f]\n", 
            min(expr_aligned, na.rm = TRUE), 
            max(expr_aligned, na.rm = TRUE)))
cat(sprintf("  - Missing: %.2f%%\n", 
            100 * sum(is.na(expr_aligned)) / prod(dim(expr_aligned))))

cat("\nDrug Response (AAC):\n")
cat(sprintf("  - Mean: %.3f\n", mean(drug_vec_clean)))
cat(sprintf("  - Median: %.3f\n", median(drug_vec_clean)))
cat(sprintf("  - Range: [%.3f, %.3f]\n", 
            min(drug_vec_clean), max(drug_vec_clean)))
cat(sprintf("  - Missing: %.2f%%\n", 
            100 * sum(is.na(drug_aligned)) / prod(dim(drug_aligned))))

cat("\nData Quality:\n")
cat(sprintf("  - Low variance genes: %d\n", n_low_var))
cat(sprintf("  - Drugs with <20%% missing: %d\n", 
            sum(drug_na_percent < 20)))
cat(sprintf("  - Drugs with >80%% missing: %d\n", 
            length(high_missing_drugs)))

cat("\nPCA Results:\n")
cat(sprintf("  - PC1 variance: %.2f%%\n", var_explained[1]))
cat(sprintf("  - PC2 variance: %.2f%%\n", var_explained[2]))
cat(sprintf("  - Top 10 PCs variance: %.2f%%\n", sum(var_explained[1:10])))

cat("\n=============================================================================\n")
cat("EDA Completed Successfully!\n")
cat("=============================================================================\n\n")

cat("Generated files:\n")
cat("  Figures (in figures/eda/):\n")
cat("    1. 01_expression_distribution.png\n")
cat("    2. 02_sample_quality.png\n")
cat("    3. 03_gene_statistics.png\n")
cat("    4. 04_drug_response_distribution.png\n")
cat("    5. 05_missing_values.png\n")
cat("    6. 06_drug_statistics.png\n")
cat("    7. 07_pca_analysis.png\n")
cat("    8. 08_sample_correlation.png\n")
cat("    9. 09_drug_correlation.png\n")
cat("    10. 10_expression_response_examples.png\n")
cat("\n  Data (in data/processed/):\n")
cat("    - quality_flags.rds\n\n")

cat("Next steps:\n")
cat("  - Review all generated plots\n")
cat("  - Check quality flags for issues\n")
cat("  - Run 03_preprocessing.R to clean and normalize data\n\n")

# =============================================================================
# End Of Script
# =============================================================================