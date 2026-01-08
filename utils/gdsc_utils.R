# =============================================================================
# GDSC Utility Functions
# =============================================================================
# Purpose: Reusable helper functions for GDSC analysis
# Author: prbfarel
# Date: 2026-01-07
# =============================================================================

library(tidyverse)

# =============================================================================
# Data Loading Utilities
# =============================================================================

#' Load GDSC expression data
#' @param path Path to expression matrix
#' @return Matrix of gene expression values
load_expression <- function(path) {
  expr <- readRDS(path)
  message(sprintf("Loaded expression: %d samples × %d genes", 
                  nrow(expr), ncol(expr)))
  return(expr)
}

#' Load GDSC drug response data
#' @param path Path to drug response matrix
#' @return Matrix of AAC values
load_drug_response <- function(path) {
  drug <- readRDS(path)
  message(sprintf("Loaded drug response: %d samples × %d drugs", 
                  nrow(drug), ncol(drug)))
  return(drug)
}

# =============================================================================
# Data Quality Checks
# =============================================================================

#' Check for missing values
#' @param data Data matrix or data frame
#' @return Summary of missing values
check_missing <- function(data) {
  missing_count <- sum(is.na(data))
  missing_pct <- 100 * missing_count / length(data)
  
  message(sprintf("Missing values: %d (%.2f%%)", 
                  missing_count, missing_pct))
  
  # Per-column missing
  col_missing <- colSums(is.na(data))
  high_missing <- names(col_missing[col_missing > nrow(data) * 0.5])
  
  if(length(high_missing) > 0) {
    warning(sprintf("%d columns have >50%% missing values", 
                    length(high_missing)))
  }
  
  return(list(
    total_missing = missing_count,
    pct_missing = missing_pct,
    high_missing_cols = high_missing
  ))
}

#' Check for outliers using IQR method
#' @param x Numeric vector
#' @param threshold IQR multiplier (default: 3)
#' @return Logical vector indicating outliers
detect_outliers_iqr <- function(x, threshold = 3) {
  q1 <- quantile(x, 0.25, na.rm = TRUE)
  q3 <- quantile(x, 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  
  lower <- q1 - threshold * iqr
  upper <- q3 + threshold * iqr
  
  outliers <- x < lower | x > upper
  return(outliers)
}

# =============================================================================
# Feature Selection Utilities
# =============================================================================

#' Calculate correlation with p-values
#' @param x Predictor vector
#' @param y Response vector
#' @return List with correlation and p-value
cor_with_pvalue <- function(x, y) {
  test <- cor.test(x, y, method = "pearson")
  return(list(
    cor = test$estimate,
    pvalue = test$p.value
  ))
}

#' Select top N features by variance
#' @param expr Expression matrix
#' @param n Number of features to select
#' @return Vector of selected feature names
select_by_variance <- function(expr, n = 5000) {
  variances <- apply(expr, 2, var, na.rm = TRUE)
  top_features <- names(sort(variances, decreasing = TRUE)[1:n])
  return(top_features)
}

# =============================================================================
# Model Evaluation Utilities
# =============================================================================

#' Calculate regression metrics
#' @param observed Observed values
#' @param predicted Predicted values
#' @return List of metrics (R², RMSE, MAE)
calculate_metrics <- function(observed, predicted) {
  # Remove NA pairs
  valid <- !is.na(observed) & !is.na(predicted)
  obs <- observed[valid]
  pred <- predicted[valid]
  
  # Calculate metrics
  ss_res <- sum((obs - pred)^2)
  ss_tot <- sum((obs - mean(obs))^2)
  r2 <- 1 - (ss_res / ss_tot)
  rmse <- sqrt(mean((obs - pred)^2))
  mae <- mean(abs(obs - pred))
  
  return(list(
    r2 = r2,
    rmse = rmse,
    mae = mae,
    n = length(obs)
  ))
}

#' Calculate classification metrics (for binary response)
#' @param observed Observed binary labels
#' @param predicted Predicted probabilities
#' @param threshold Classification threshold
#' @return List of metrics
calculate_classification_metrics <- function(observed, predicted, 
                                             threshold = 0.5) {
  pred_class <- ifelse(predicted > threshold, 1, 0)
  
  tp <- sum(observed == 1 & pred_class == 1)
  tn <- sum(observed == 0 & pred_class == 0)
  fp <- sum(observed == 0 & pred_class == 1)
  fn <- sum(observed == 1 & pred_class == 0)
  
  accuracy <- (tp + tn) / (tp + tn + fp + fn)
  precision <- tp / (tp + fp)
  recall <- tp / (tp + fn)
  f1 <- 2 * (precision * recall) / (precision + recall)
  
  return(list(
    accuracy = accuracy,
    precision = precision,
    recall = recall,
    f1 = f1,
    confusion_matrix = matrix(c(tn, fp, fn, tp), nrow = 2)
  ))
}

# =============================================================================
# Visualization Utilities
# =============================================================================

#' Create publication-quality theme for ggplot2
#' @return ggplot2 theme object
theme_publication <- function() {
  theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.title = element_text(face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey90")
    )
}

#' Save plot with consistent settings
#' @param plot ggplot2 object
#' @param filename Output filename
#' @param width Width in inches
#' @param height Height in inches
save_plot <- function(plot, filename, width = 10, height = 6) {
  ggsave(filename, plot, width = width, height = height, 
         dpi = 300, units = "in")
  message(sprintf("Saved: %s", filename))
}

# =============================================================================
# File Management Utilities
# =============================================================================

#' Create directory if it doesn't exist
#' @param path Directory path
ensure_dir <- function(path) {
  if(!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    message(sprintf("Created directory: %s", path))
  }
}

#' Save object with message
#' @param object R object to save
#' @param path Output path
save_with_message <- function(object, path) {
  saveRDS(object, path)
  message(sprintf("Saved: %s", path))
}

# =============================================================================
# Logging Utilities
# =============================================================================

#' Print section header
#' @param text Header text
print_header <- function(text) {
  line <- paste(rep("=", 77), collapse = "")
  cat("\n", line, "\n", text, "\n", line, "\n\n", sep = "")
}

#' Print step message
#' @param step_num Step number
#' @param text Step description
print_step <- function(step_num, text) {
  cat(sprintf("Step %d: %s\n", step_num, text))
}

# =============================================================================
# Data Transformation Utilities
# =============================================================================

#' Z-score normalization
#' @param x Numeric vector or matrix
#' @return Normalized values
zscore <- function(x) {
  if(is.matrix(x)) {
    return(scale(x))
  } else {
    return((x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
  }
}

#' Min-max normalization
#' @param x Numeric vector
#' @param range Target range (default: 0-1)
#' @return Normalized values
minmax_normalize <- function(x, range = c(0, 1)) {
  min_val <- min(x, na.rm = TRUE)
  max_val <- max(x, na.rm = TRUE)
  
  normalized <- (x - min_val) / (max_val - min_val)
  normalized <- normalized * (range[2] - range[1]) + range[1]
  
  return(normalized)
}

# =============================================================================
# Statistical Utilities
# =============================================================================

#' Perform permutation test
#' @param x Group 1 values
#' @param y Group 2 values
#' @param n_permutations Number of permutations
#' @return P-value from permutation test
permutation_test <- function(x, y, n_permutations = 1000) {
  observed_diff <- mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE)
  
  combined <- c(x, y)
  n_x <- length(x)
  
  perm_diffs <- replicate(n_permutations, {
    shuffled <- sample(combined)
    mean(shuffled[1:n_x]) - mean(shuffled[(n_x+1):length(shuffled)])
  })
  
  p_value <- mean(abs(perm_diffs) >= abs(observed_diff))
  
  return(list(
    observed_diff = observed_diff,
    p_value = p_value,
    null_distribution = perm_diffs
  ))
}

# =============================================================================
# Export Message
# =============================================================================

message("✓ GDSC utility functions loaded successfully")
message("Available functions:")
message("  Data loading: load_expression(), load_drug_response()")
message("  Quality checks: check_missing(), detect_outliers_iqr()")
message("  Feature selection: select_by_variance(), cor_with_pvalue()")
message("  Evaluation: calculate_metrics(), calculate_classification_metrics()")
message("  Visualization: theme_publication(), save_plot()")
message("  File management: ensure_dir(), save_with_message()")
message("  Logging: print_header(), print_step()")
message("  Transformation: zscore(), minmax_normalize()")
message("  Statistics: permutation_test()")