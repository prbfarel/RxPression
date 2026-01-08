# ============================================
# Script 06: Model Evaluation & Visualization
# Deskripsi: Script untuk evaluasi model pada test set dan membuat visualisasi 
#            publication-quality untuk hasil machine learning
# ============================================
# Input: 
#   - Trained models (from 05_train.R)
#   - Test data (expr_test.rds, drug_test.rds)
# Output: 
#   - Test set predictions
#   - Performance metrics (train vs test)
#   - Publication quality figures
#   - Comprehensive report
# Author: prbfarel
# Date: 2026-01-06
# ============================================

# Load required packages
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(pheatmap)
library(RColorBrewer)
library(scales)
library(caret)

# Set random seed
set.seed(123)

# Create output directories
dir.create("figures/evaluation", recursive = TRUE, showWarnings = FALSE)
dir.create("results/predictions", recursive = TRUE, showWarnings = FALSE)
dir.create("reports", recursive = TRUE, showWarnings = FALSE)

# =================================
# Step 1: Load Data & Models
# =================================

cat("Step 1: Loading data and models...\n")

# Load test data
expr_test <- readRDS("data/processed/expr_test.rds")
drug_test <- readRDS("data/processed/drug_test.rds")

cat(sprintf("  Test data: %d samples × %d genes\n", 
            nrow(expr_test), ncol(expr_test)))
cat(sprintf("  Test drugs: %d samples × %d drugs\n", 
            nrow(drug_test), ncol(drug_test)))

# Load training results
training_log <- readRDS("models/performance/training_log.rds")

cat(sprintf("  Training log: %d drugs trained\n", nrow(training_log)))

# Load feature selection (if used)
if(file.exists("data/processed/selected_features_per_drug.rds")) {
  selected_features <- readRDS("data/processed/selected_features_per_drug.rds")
  use_features <- TRUE
  cat("  ✓ Using selected features\n")
} else {
  use_features <- FALSE
  cat("  ⚠ No feature selection file found\n")
}

cat("\n")

# =================================
# Step 2: Predict on Test Set
# =================================

cat("Step 2: Making predictions on test set...\n")

# Initialize results storage
test_predictions <- data.frame()
test_performance <- data.frame()

# Get list of trained models
model_files <- list.files("models/trained/", pattern = "*.rds", full.names = TRUE)

cat(sprintf("  Found %d trained models\n", length(model_files)))

# Predict for each drug
for(i in seq_len(nrow(training_log))) {
  drug <- training_log$drug[i]
  
  if(i %% 20 == 0) {
    cat(sprintf("  Processing drug %d/%d...\n", i, nrow(training_log)))
  }
  
  # Get response
  y_test <- drug_test[, drug]
  valid_idx <- !is.na(y_test)
  
  if(sum(valid_idx) < 10) {
    next  # Skip if too few test samples
  }
  
  y_test_valid <- y_test[valid_idx]
  
  # Get features
  if(use_features && drug %in% names(selected_features)) {
    features <- selected_features[[drug]]
    if(is.list(features)) features <- unlist(features)
    features <- as.character(features)
    features <- intersect(features, colnames(expr_test))
  } else {
    features <- colnames(expr_test)
  }
  
  if(length(features) < 10) {
    next
  }
  
  X_test_valid <- expr_test[valid_idx, features]
  
  # Load best model for this drug
  best_algo <- training_log$best_algorithm[i]
  model_file <- sprintf("models/trained/%s_%s.rds", 
                        gsub("[^[:alnum:]]", "_", drug), 
                        best_algo)
  
  if(!file.exists(model_file)) {
    next
  }
  
  model <- readRDS(model_file)
  
  # Make predictions
  tryCatch({
    predictions <- predict(model, newdata = X_test_valid)
    
    # Calculate test metrics
    test_r2 <- cor(predictions, y_test_valid)^2
    test_rmse <- sqrt(mean((predictions - y_test_valid)^2))
    test_mae <- mean(abs(predictions - y_test_valid))
    
    # Store predictions
    pred_df <- data.frame(
      drug = drug,
      sample = rownames(X_test_valid),
      observed = y_test_valid,
      predicted = predictions
    )
    test_predictions <- rbind(test_predictions, pred_df)
    
    # Store performance
    perf_df <- data.frame(
      drug = drug,
      algorithm = best_algo,
      n_test_samples = sum(valid_idx),
      n_features = length(features),
      train_r2 = training_log$best_rsquared[i],
      train_rmse = training_log$best_rmse[i],
      train_mae = training_log$best_mae[i],
      test_r2 = test_r2,
      test_rmse = test_rmse,
      test_mae = test_mae
    )
    test_performance <- rbind(test_performance, perf_df)
    
  }, error = function(e) {
    # Skip if prediction fails
  })
}

cat(sprintf("\n  ✓ Predictions complete for %d drugs\n", 
            nrow(test_performance)))
cat(sprintf("  ✓ Total predictions: %d\n\n", 
            nrow(test_predictions)))


# =====================================
# Step 3: Calculate Performance Metrics
# =====================================

cat("Step 3: Calculating performance metrics...\n")

# Overall statistics
train_r2_mean <- mean(test_performance$train_r2, na.rm = TRUE)
train_r2_sd <- sd(test_performance$train_r2, na.rm = TRUE)
test_r2_mean <- mean(test_performance$test_r2, na.rm = TRUE)
test_r2_sd <- sd(test_performance$test_r2, na.rm = TRUE)

cat("\nPerformance Summary:\n")
cat(sprintf("  Training R²: %.3f ± %.3f\n", train_r2_mean, train_r2_sd))
cat(sprintf("  Test R²:     %.3f ± %.3f\n", test_r2_mean, test_r2_sd))
cat(sprintf("  Generalization gap: %.3f\n", train_r2_mean - test_r2_mean))

# Check for overfitting
test_performance$overfit <- test_performance$train_r2 - test_performance$test_r2
overfit_drugs <- sum(test_performance$overfit > 0.1)

cat(sprintf("\n  Overfitting check:\n"))
cat(sprintf("    Drugs with gap >0.1: %d/%d (%.1f%%)\n", 
            overfit_drugs, nrow(test_performance),
            100*overfit_drugs/nrow(test_performance)))

# Per-algorithm performance
algo_perf <- test_performance %>%
  group_by(algorithm) %>%
  summarise(
    count = n(),
    train_r2 = mean(train_r2, na.rm = TRUE),
    test_r2 = mean(test_r2, na.rm = TRUE),
    gap = mean(train_r2 - test_r2, na.rm = TRUE)
  )

cat("\nPer-Algorithm Performance:\n")
print(algo_perf)
cat("\n")


# =================================
# Step 4: Save Results
# =================================

cat("Step 4: Saving results...\n")

# Save predictions
saveRDS(test_predictions, "results/predictions/test_predictions.rds")
write.csv(test_predictions, "results/predictions/test_predictions.csv", 
          row.names = FALSE)

# Save performance
saveRDS(test_performance, "results/predictions/test_performance.rds")
write.csv(test_performance, "results/predictions/test_performance.csv", 
          row.names = FALSE)

cat("  ✓ test_predictions.rds/csv saved\n")
cat("  ✓ test_performance.rds/csv saved\n\n")

# =================================
# Step 5: Create Visualizations
# =================================

cat("Step 5: Creating visualizations...\n\n")

# 5.1 Train vs Test R² Comparison
cat("  Creating train vs test comparison...\n")

png("figures/evaluation/01_train_vs_test_r2.png",
    width = 10, height = 8, units = "in", res = 300)

ggplot(test_performance, aes(x = train_r2, y = test_r2)) +
  geom_point(aes(color = algorithm), size = 3, alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", 
              color = "black", size = 1) +
  geom_smooth(method = "lm", se = TRUE, color = "blue", 
              linetype = "dotted") +
  scale_color_brewer(palette = "Set1", name = "Algorithm") +
  labs(title = "Training vs Test Set Performance",
       subtitle = sprintf("Mean train R² = %.3f, Mean test R² = %.3f", 
                          train_r2_mean, test_r2_mean),
       x = "Training R² (5-fold CV)",
       y = "Test R² (Hold-out Set)") +
  xlim(0, max(test_performance$train_r2) * 1.1) +
  ylim(0, max(test_performance$test_r2) * 1.1) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

dev.off()
cat("    ✓ Saved: 01_train_vs_test_r2.png\n")

# 5.2 R² Distribution (Train vs Test)
cat("  Creating R² distributions...\n")

png("figures/evaluation/02_r2_distributions.png",
    width = 12, height = 6, units = "in", res = 300)

r2_long <- test_performance %>%
  select(drug, train_r2, test_r2) %>%
  pivot_longer(cols = c(train_r2, test_r2), 
               names_to = "dataset", 
               values_to = "r2")

ggplot(r2_long, aes(x = r2, fill = dataset)) +
  geom_histogram(alpha = 0.6, bins = 30, position = "identity") +
  geom_vline(data = r2_long %>% group_by(dataset) %>% 
               summarise(mean_r2 = mean(r2, na.rm = TRUE)),
             aes(xintercept = mean_r2, color = dataset),
             linetype = "dashed", size = 1) +
  scale_fill_manual(values = c("train_r2" = "steelblue", 
                               "test_r2" = "coral"),
                    labels = c("Training (CV)", "Test (Hold-out)"),
                    name = "Dataset") +
  scale_color_manual(values = c("train_r2" = "steelblue", 
                                "test_r2" = "coral"),
                     labels = c("Training (CV)", "Test (Hold-out)"),
                     name = "Mean") +
  labs(title = "Distribution of R² Scores",
       subtitle = "Comparison between training and test performance",
       x = "R² Score", y = "Number of Drugs") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

dev.off()
cat("    ✓ Saved: 02_r2_distributions.png\n")

# 5.3 Predicted vs Observed (Best drugs)
cat("  Creating prediction plots...\n")

# Select top 6 drugs by test R²
top_drugs <- test_performance %>%
  arrange(desc(test_r2)) %>%
  head(6) %>%
  pull(drug)

png("figures/evaluation/03_predictions_top_drugs.png",
    width = 14, height = 10, units = "in", res = 300)

top_predictions <- test_predictions %>%
  filter(drug %in% top_drugs)

ggplot(top_predictions, aes(x = observed, y = predicted)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "red", size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  facet_wrap(~drug, scales = "free", ncol = 3) +
  labs(title = "Predicted vs Observed Drug Response",
       subtitle = "Top 6 drugs by test set R²",
       x = "Observed AAC",
       y = "Predicted AAC") +
  theme_bw(base_size = 10)

dev.off()
cat("    ✓ Saved: 03_predictions_top_drugs.png\n")

# 5.4 Algorithm Comparison
cat("  Creating algorithm comparison...\n")

png("figures/evaluation/04_algorithm_comparison.png",
    width = 12, height = 6, units = "in", res = 300)

test_perf_long <- test_performance %>%
  select(drug, algorithm, train_r2, test_r2) %>%
  pivot_longer(cols = c(train_r2, test_r2), 
               names_to = "dataset", 
               values_to = "r2")

ggplot(test_perf_long, aes(x = algorithm, y = r2, fill = dataset)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = c("train_r2" = "steelblue", 
                               "test_r2" = "coral"),
                    labels = c("Training", "Test"),
                    name = "Dataset") +
  labs(title = "Algorithm Performance Comparison",
       subtitle = "Training vs Test set performance by algorithm",
       x = "Algorithm", y = "R² Score") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()
cat("    ✓ Saved: 04_algorithm_comparison.png\n")

# 5.5 Overfitting Analysis
cat("  Creating overfitting analysis...\n")

png("figures/evaluation/05_overfitting_analysis.png",
    width = 12, height = 6, units = "in", res = 300)

test_performance$overfit_category <- cut(
  test_performance$overfit,
  breaks = c(-Inf, 0, 0.05, 0.1, 0.2, Inf),
  labels = c("Better on test", "Minimal (<0.05)", 
             "Moderate (0.05-0.1)", "High (0.1-0.2)", 
             "Severe (>0.2)")
)

ggplot(test_performance, aes(x = overfit_category, fill = overfit_category)) +
  geom_bar() +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  scale_fill_brewer(palette = "RdYlGn", direction = -1, 
                    name = "Overfitting") +
  labs(title = "Overfitting Analysis",
       subtitle = "Distribution of train-test performance gap",
       x = "Overfitting Category (Train R² - Test R²)",
       y = "Number of Drugs") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()
cat("    ✓ Saved: 05_overfitting_analysis.png\n")

# 5.6 Performance vs Features
cat("  Creating feature analysis...\n")

png("figures/evaluation/06_performance_vs_features.png",
    width = 10, height = 6, units = "in", res = 300)

ggplot(test_performance, aes(x = n_features, y = test_r2)) +
  geom_point(aes(color = algorithm), size = 3, alpha = 0.7) +
  geom_smooth(method = "loess", se = TRUE, color = "black", 
              linetype = "dashed") +
  scale_color_brewer(palette = "Set1", name = "Algorithm") +
  labs(title = "Test Performance vs Number of Features",
       subtitle = "Effect of feature count on predictive performance",
       x = "Number of Features Used",
       y = "Test R²") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

dev.off()
cat("    ✓ Saved: 06_performance_vs_features.png\n\n")


# =====================================
# Step 6: Generate Comprehensice Report
# =====================================

cat("Step 6: Generating evaluation report...\n")

report <- c(
  "=============================================================================",
  "GDSC Model Evaluation Report",
  "=============================================================================",
  "",
  paste("Date:", Sys.time()),
  "",
  "Data Summary:",
  sprintf("  - Training samples: 640"),
  sprintf("  - Test samples: 161"),
  sprintf("  - Drugs evaluated: %d", nrow(test_performance)),
  sprintf("  - Total predictions: %d", nrow(test_predictions)),
  "",
  "Performance Metrics:",
  "",
  "Training Set (5-fold CV):",
  sprintf("  - Mean R²:   %.3f ± %.3f", train_r2_mean, train_r2_sd),
  sprintf("  - Median R²: %.3f", median(test_performance$train_r2, na.rm = TRUE)),
  sprintf("  - Range:     [%.3f, %.3f]", 
          min(test_performance$train_r2, na.rm = TRUE),
          max(test_performance$train_r2, na.rm = TRUE)),
  "",
  "Test Set (Hold-out):",
  sprintf("  - Mean R²:   %.3f ± %.3f", test_r2_mean, test_r2_sd),
  sprintf("  - Median R²: %.3f", median(test_performance$test_r2, na.rm = TRUE)),
  sprintf("  - Range:     [%.3f, %.3f]", 
          min(test_performance$test_r2, na.rm = TRUE),
          max(test_performance$test_r2, na.rm = TRUE)),
  "",
  "Generalization:",
  sprintf("  - Mean train-test gap: %.3f", train_r2_mean - test_r2_mean),
  sprintf("  - Correlation (train vs test R²): %.3f",
          cor(test_performance$train_r2, test_performance$test_r2)),
  "",
  "Overfitting Analysis:",
  "",
  sprintf("  Drugs with excellent generalization (<0.05 gap): %d (%.1f%%)",
          sum(test_performance$overfit < 0.05),
          100*sum(test_performance$overfit < 0.05)/nrow(test_performance)),
  sprintf("  Drugs with moderate overfitting (0.05-0.1): %d (%.1f%%)",
          sum(test_performance$overfit >= 0.05 & test_performance$overfit < 0.1),
          100*sum(test_performance$overfit >= 0.05 & test_performance$overfit < 0.1)/nrow(test_performance)),
  sprintf("  Drugs with high overfitting (>0.1): %d (%.1f%%)",
          sum(test_performance$overfit >= 0.1),
          100*sum(test_performance$overfit >= 0.1)/nrow(test_performance)),
  "",
  "Algorithm Performance:",
  "",
  paste(capture.output(print(algo_perf)), collapse = "\n"),
  "",
  "Top 10 Best Performing Drugs (by test R²):",
  "",
  paste(sprintf("  %d. %s: Test R² = %.3f (Train R² = %.3f)",
                1:10,
                head(test_performance[order(-test_performance$test_r2), "drug"], 10),
                head(test_performance[order(-test_performance$test_r2), "test_r2"], 10),
                head(test_performance[order(-test_performance$test_r2), "train_r2"], 10)),
        collapse = "\n"),
  "",
  "Top 10 Worst Performing Drugs (by test R²):",
  "",
  paste(sprintf("  %d. %s: Test R² = %.3f (Train R² = %.3f)",
                1:10,
                head(test_performance[order(test_performance$test_r2), "drug"], 10),
                head(test_performance[order(test_performance$test_r2), "test_r2"], 10),
                head(test_performance[order(test_performance$test_r2), "train_r2"], 10)),
        collapse = "\n"),
  "",
  "Output Files:",
  "",
  "  Results:",
  "    - test_predictions.rds/csv",
  "    - test_performance.rds/csv",
  "",
  "  Figures:",
  "    - 01_train_vs_test_r2.png",
  "    - 02_r2_distributions.png",
  "    - 03_predictions_top_drugs.png",
  "    - 04_algorithm_comparison.png",
  "    - 05_overfitting_analysis.png",
  "    - 06_performance_vs_features.png",
  "",
  "=============================================================================",
  "Evaluation Completed Succesfully!",
  "=============================================================================",
  ""
)

writeLines(report, "reports/evaluation_report.txt")
cat("  ✓ Report saved to reports/evaluation_report.txt\n\n")

# =============================================================================
# SUMMARY
# =============================================================================

cat("=============================================================================\n")
cat("Evaluation Summary\n")
cat("=============================================================================\n\n")

cat(sprintf("Drugs evaluated: %d\n", nrow(test_performance)))
cat(sprintf("Test set R²: %.3f ± %.3f\n", test_r2_mean, test_r2_sd))
cat(sprintf("Generalization gap: %.3f\n", train_r2_mean - test_r2_mean))

if(train_r2_mean - test_r2_mean < 0.05) {
  cat("\n✓ EXCELLENT: Models generalize very well!\n")
} else if(train_r2_mean - test_r2_mean < 0.1) {
  cat("\n✓ GOOD: Models show acceptable generalization\n")
} else {
  cat("\n⚠ MODERATE: Some overfitting detected\n")
}

cat("\nGenerated files:\n")
cat("  Results:\n")
cat("    - test_predictions.rds/csv\n")
cat("    - test_performance.rds/csv\n")
cat("\n  Figures:\n")
cat("    - 01_train_vs_test_r2.png\n")
cat("    - 02_r2_distributions.png\n")
cat("    - 03_predictions_top_drugs.png\n")
cat("    - 04_algorithm_comparison.png\n")
cat("    - 05_overfitting_analysis.png\n")
cat("    - 06_performance_vs_features.png\n")
cat("\n  Report:\n")
cat("    - evaluation_report.txt\n\n")

cat("=============================================================================\n")
cat("Evaluation Completed Succesfully!\n")
cat("=============================================================================\n")

# =============================================================================
# End of Script
# =============================================================================