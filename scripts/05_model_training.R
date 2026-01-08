# ============================================
# Script 05: Model Training Fixed
# Deskripsi: Script untuk training machine learning models untuk prediksi
#            drug response menggunakan gene expression data
# ============================================

# Input: 
#   - expr_train.rds, drug_train.rds
#   - selected_features_per_drug.rds (dari feature selection)
# Output: 
#   - Trained models per drug
#   - Performance metrics
#   - Training report
# Author: prbfarel
# Date: 2026-01-02
# ============================================

# Load required packages
library(tidyverse)
library(caret)
library(glmnet)        # Elastic Net
library(randomForest)  # Random Forest
library(kernlab)       # SVM
library(parallel)      # Parallel processing
library(doParallel)    # Parallel backend

# Create output directories
dir.create("models/trained", recursive = TRUE, showWarnings = FALSE)
dir.create("models/performance", recursive = TRUE, showWarnings = FALSE)
dir.create("reports", recursive = TRUE, showWarnings = FALSE)

# =================================
# Step 0: Configuration
# =================================

# Training configuration
config <- list(
  # Algorithms to train
  algorithms = c("elastic_net", "random_forest", "svm"),
  
  # Cross-validation settings
  cv_folds = 5,
  cv_repeats = 1,
  
  # Number of drugs to train (for testing, set to small number first)
  n_drugs_to_train = "all",  # Change to "all" for full training
  
  # Parallel processing
  use_parallel = TRUE,
  n_cores = detectCores() - 1,
  
  # Save models
  save_models = TRUE
)

cat("Training Configuration:\n")
cat(sprintf("  - Algorithms: %s\n", paste(config$algorithms, collapse = ", ")))
cat(sprintf("  - CV folds: %d\n", config$cv_folds))
cat(sprintf("  - Drugs to train: %s\n", 
            ifelse(config$n_drugs_to_train == "all", "ALL", config$n_drugs_to_train)))
cat(sprintf("  - Parallel cores: %d\n\n", config$n_cores))

# =================================
# Step 1: Load Data
# =================================

cat("Step 1: Loading data...\n")

expr_train <- readRDS("data/processed/expr_train.rds")
drug_train <- readRDS("data/processed/drug_train.rds")

# Load feature selection results (if available)
if(file.exists("data/processed/selected_features_per_drug.rds")) {
  selected_features <- readRDS("data/processed/selected_features_per_drug.rds")
  use_feature_selection <- TRUE
  cat("  ✓ Using selected features from feature selection\n")
} else {
  use_feature_selection <- FALSE
  cat("  ⚠ No feature selection file found, using all genes\n")
}

cat(sprintf("  - Expression: %d samples × %d genes\n", 
            nrow(expr_train), ncol(expr_train)))
cat(sprintf("  - Drug response: %d samples × %d drugs\n\n", 
            nrow(drug_train), ncol(drug_train)))

# =================================
# Step 2: Prepare Training
# =================================

cat("Step 2: Preparing training setup...\n")

# Select drugs to train
if(is.character(config$n_drugs_to_train) && config$n_drugs_to_train == "all") {
  drugs_to_train <- colnames(drug_train)
} else {
  n_drugs <- min(config$n_drugs_to_train, ncol(drug_train))
  drugs_to_train <- colnames(drug_train)[1:n_drugs]
}

cat(sprintf("  - Number of drugs to train: %d\n", length(drugs_to_train)))

# Setup parallel processing
if(config$use_parallel) {
  cl <- makeCluster(config$n_cores)
  registerDoParallel(cl)
  cat(sprintf("  ✓ Parallel processing enabled (%d cores)\n", config$n_cores))
} else {
  cat("  - Serial processing (no parallel)\n")
}

# Setup cross-validation
train_control <- trainControl(
  method = "cv",
  number = config$cv_folds,
  savePredictions = "final",
  returnResamp = "final",
  allowParallel = FALSE  # Disable parallel in trainControl
)

cat("  ✓ Cross-validation setup complete\n\n")

# =================================
# Step 3: Define Training Functions
# =================================

cat("Step 3: Preparing training setup...\n")

# Function to train Elastic Net
train_elastic_net <- function(X, y, train_control) {
  
  # Hyperparameter grid
  tune_grid <- expand.grid(
    alpha = seq(0, 1, 0.25),  # 0=Ridge, 1=LASSO, 0.5=Elastic Net
    lambda = 10^seq(-3, 1, length.out = 20)
  )
  
  # Train model
  model <- train(
    x = X,
    y = y,
    method = "glmnet",
    trControl = train_control,
    tuneGrid = tune_grid,
    metric = "RMSE"
  )
  
  return(model)
}

# Function to train Random Forest
train_random_forest <- function(X, y, train_control) {
  
  # Hyperparameter grid
  tune_grid <- expand.grid(
    mtry = c(sqrt(ncol(X)), ncol(X)/3, ncol(X)/2)
  )
  
  # Train model
  model <- train(
    x = X,
    y = y,
    method = "rf",
    trControl = train_control,
    tuneGrid = tune_grid,
    metric = "RMSE",
    ntree = 100,  # Number of trees
    importance = TRUE
  )
  
  return(model)
}

# Function to train Support Vector Machine (SVM)
train_svm <- function(X, y, train_control) {
  
  # Hyperparameter grid for SVM with RBF kernel
  tune_grid <- expand.grid(
    C = c(0.1, 1, 10),           # Regularization parameter
    sigma = c(0.001, 0.01, 0.1)  # Kernel bandwidth
  )
  
  # Train model
  model <- train(
    x = X,
    y = y,
    method = "svmRadial",  # SVM with Radial Basis Function kernel
    trControl = train_control,
    tuneGrid = tune_grid,
    metric = "RMSE"
  )
  
  return(model)
}


# ==================================
# Step 4: Train Models for Each Drug
# ==================================

cat("Step 4: Training models...\n\n")

# Initialize results storage
all_results <- list()
training_log <- data.frame()

# Start training timer
start_time <- Sys.time()

for(i in seq_along(drugs_to_train)) {
  drug <- drugs_to_train[i]
  
  cat(sprintf("=== Drug %d/%d: %s ===\n", i, length(drugs_to_train), drug))
  
  # Get response variable
  y <- drug_train[, drug]
  
  # Remove samples with NA (if any)
  valid_idx <- !is.na(y)
  
  if(sum(valid_idx) < 50) {
    cat("  ⚠ Skipping: Too few samples (<50)\n\n")
    next
  }
  
  y_valid <- y[valid_idx]
  
  # Get features for this drug
  if(use_feature_selection && drug %in% names(selected_features)) {
    features <- selected_features[[drug]]
    
    # Check if features is a list (from feature selection result)
    if(is.list(features)) {
      features <- unlist(features)
    }
    
    # Ensure features are character vector
    features <- as.character(features)
    
    # Ensure features exist in expr_train
    features <- intersect(features, colnames(expr_train))
  } else {
    features <- colnames(expr_train)
  }
  
  if(length(features) < 10) {
    cat("  ⚠ Skipping: Too few features (<10)\n\n")
    next
  }
  
  X_valid <- expr_train[valid_idx, features]
  
  cat(sprintf("  - Samples: %d\n", nrow(X_valid)))
  cat(sprintf("  - Features: %d\n", ncol(X_valid)))
  
  # Train each algorithm
  drug_results <- list()
  
  for(algo in config$algorithms) {
    cat(sprintf("  Training %s...\n", algo))
    
    algo_start <- Sys.time()
    
    tryCatch({
      if(algo == "elastic_net") {
        model <- train_elastic_net(X_valid, y_valid, train_control)
      } else if(algo == "random_forest") {
        model <- train_random_forest(X_valid, y_valid, train_control)
      } else if(algo == "svm") {
        model <- train_svm(X_valid, y_valid, train_control)
      }
      
      # Get best performance
      best_rmse <- min(model$results$RMSE, na.rm = TRUE)
      best_rsquared <- max(model$results$Rsquared, na.rm = TRUE)
      best_mae <- min(model$results$MAE, na.rm = TRUE)
      
      # Store results
      drug_results[[algo]] <- list(
        model = model,
        rmse = best_rmse,
        rsquared = best_rsquared,
        mae = best_mae,
        training_time = as.numeric(difftime(Sys.time(), algo_start, units = "secs"))
      )
      
      cat(sprintf("    ✓ R² = %.3f, RMSE = %.3f (%.1fs)\n", 
                  best_rsquared, best_rmse, drug_results[[algo]]$training_time))
      
    }, error = function(e) {
      cat(sprintf("    ✗ Training failed: %s\n", e$message))
      drug_results[[algo]] <- NULL
    })
  }
  
  # Select best model
  if(length(drug_results) > 0) {
    best_algo <- names(which.max(sapply(drug_results, function(x) x$rsquared)))
    best_model <- drug_results[[best_algo]]
    
    cat(sprintf("  ★ Best: %s (R² = %.3f)\n\n", best_algo, best_model$rsquared))
    
    # Store results
    all_results[[drug]] <- drug_results
    
    # Log performance
    training_log <- rbind(training_log, data.frame(
      drug = drug,
      best_algorithm = best_algo,
      n_samples = nrow(X_valid),
      n_features = ncol(X_valid),
      best_rmse = best_model$rmse,
      best_rsquared = best_model$rsquared,
      best_mae = best_model$mae,
      elastic_net_r2 = ifelse("elastic_net" %in% names(drug_results), 
                              drug_results$elastic_net$rsquared, NA),
      rf_r2 = ifelse("random_forest" %in% names(drug_results), 
                     drug_results$random_forest$rsquared, NA),
      svm_r2 = ifelse("svm" %in% names(drug_results), 
                      drug_results$svm$rsquared, NA)
    ))
    
    # Save best model
    if(config$save_models) {
      model_filename <- sprintf("models/trained/%s_%s.rds", 
                                gsub("[^[:alnum:]]", "_", drug), 
                                best_algo)
      saveRDS(best_model$model, model_filename)
    }
  } else {
    cat("  ✗ All algorithms failed for this drug\n\n")
  }
}

# Stop parallel cluster
if(config$use_parallel) {
  stopCluster(cl)
}

# Training duration
total_time <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))

cat("\n=============================================================================\n")
cat(sprintf("Training Completed in %.1f minutes\n", total_time))
cat("=============================================================================\n\n")


# ==================================
# Step 5: Analyze Results
# ==================================

cat("Step 5: Analyzing training results...\n")

# Summary statistics
cat("\nPerformance Summary:\n")
cat(sprintf("  - Drugs successfully trained: %d/%d\n", 
            nrow(training_log), length(drugs_to_train)))
cat(sprintf("  - Average R²: %.3f ± %.3f\n", 
            mean(training_log$best_rsquared, na.rm = TRUE),
            sd(training_log$best_rsquared, na.rm = TRUE)))
cat(sprintf("  - Average RMSE: %.3f ± %.3f\n", 
            mean(training_log$best_rmse, na.rm = TRUE),
            sd(training_log$best_rmse, na.rm = TRUE)))
cat(sprintf("  - Average MAE: %.3f ± %.3f\n", 
            mean(training_log$best_mae, na.rm = TRUE),
            sd(training_log$best_mae, na.rm = TRUE)))

# Algorithm comparison
cat("\nAlgorithm Performance:\n")
algo_comparison <- training_log %>%
  summarise(
    elastic_net_r2 = mean(elastic_net_r2, na.rm = TRUE),
    rf_r2 = mean(rf_r2, na.rm = TRUE),
    svm_r2 = mean(svm_r2, na.rm = TRUE)
  )
print(algo_comparison)

# Best algorithm frequency
cat("\nBest Algorithm Frequency:\n")
print(table(training_log$best_algorithm))

cat("\n")


# ==================================
# Step 6: Save Results
# ==================================

cat("Step 6: Saving results...\n")

# Save training log
saveRDS(training_log, "models/performance/training_log.rds")
write.csv(training_log, "models/performance/training_log.csv", row.names = FALSE)

# Save all results
saveRDS(all_results, "models/performance/all_models_results.rds")

cat("  ✓ training_log.rds saved\n")
cat("  ✓ training_log.csv saved\n")
cat("  ✓ all_models_results.rds saved\n")

if(config$save_models) {
  cat(sprintf("  ✓ %d best models saved\n", nrow(training_log)))
}

cat("\n")


# ==================================
# Step 7: Visualize Results
# ==================================

cat("Step 7: Creating visualizations...\n")

# 7.1 R² distribution
png("models/performance/01_rsquared_distribution.png",
    width = 10, height = 6, units = "in", res = 300)

ggplot(training_log, aes(x = best_rsquared)) +
  geom_histogram(bins = 20, fill = "steelblue", color = "white") +
  geom_vline(xintercept = mean(training_log$best_rsquared), 
             color = "red", linetype = "dashed", size = 1) +
  labs(title = "Distribution of Model Performance (R²)",
       subtitle = sprintf("Mean R² = %.3f", mean(training_log$best_rsquared)),
       x = "R² Score", y = "Number of Drugs") +
  theme_bw(base_size = 12)

dev.off()
cat("  ✓ Saved: 01_rsquared_distribution.png\n")

# 7.2 Algorithm comparison
png("models/performance/02_algorithm_comparison.png",
    width = 10, height = 6, units = "in", res = 300)

algo_long <- training_log %>%
  select(drug, elastic_net_r2, rf_r2, svm_r2) %>%
  pivot_longer(cols = -drug, names_to = "algorithm", values_to = "r2") %>%
  filter(!is.na(r2))

ggplot(algo_long, aes(x = algorithm, y = r2, fill = algorithm)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Algorithm Performance Comparison",
       x = "Algorithm", y = "R² Score") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

dev.off()
cat("  ✓ Saved: 02_algorithm_comparison.png\n")

# 7.3 Performance vs Features
png("models/performance/03_performance_vs_features.png",
    width = 10, height = 6, units = "in", res = 300)

ggplot(training_log, aes(x = n_features, y = best_rsquared)) +
  geom_point(aes(color = best_algorithm), size = 3, alpha = 0.7) +
  geom_smooth(method = "loess", se = TRUE, color = "black", linetype = "dashed") +
  scale_color_brewer(palette = "Set1", name = "Best Algorithm") +
  labs(title = "Model Performance vs Number of Features",
       x = "Number of Features", y = "R² Score") +
  theme_bw(base_size = 12)

dev.off()
cat("  ✓ Saved: 03_performance_vs_features.png\n\n")


# ==================================
# Step 8: Generate Report
# ==================================

cat("Step 8: Generating training report...\n")

report <- c(
  "=============================================================================",
  "GDSC Model Training REport",
  "=============================================================================",
  "",
  paste("Date:", Sys.time()),
  paste("Training duration:", sprintf("%.1f minutes", total_time)),
  "",
  "Configuration:",
  sprintf("  - Algorithms: %s", paste(config$algorithms, collapse = ", ")),
  sprintf("  - CV folds: %d", config$cv_folds),
  sprintf("  - Parallel cores: %d", config$n_cores),
  "",
  "Data:",
  sprintf("  - Training samples: %d", nrow(expr_train)),
  sprintf("  - Total genes: %d", ncol(expr_train)),
  sprintf("  - Total drugs: %d", ncol(drug_train)),
  sprintf("  - Drugs trained: %d", length(drugs_to_train)),
  "",
  "Results:",
  sprintf("  - Successful training: %d/%d drugs", 
          nrow(training_log), length(drugs_to_train)),
  "",
  "Performance Metrics:",
  "",
  "  R² Score:",
  sprintf("    - Mean: %.3f", mean(training_log$best_rsquared)),
  sprintf("    - Median: %.3f", median(training_log$best_rsquared)),
  sprintf("    - SD: %.3f", sd(training_log$best_rsquared)),
  sprintf("    - Min: %.3f", min(training_log$best_rsquared)),
  sprintf("    - Max: %.3f", max(training_log$best_rsquared)),
  "",
  "  RMSE:",
  sprintf("    - Mean: %.3f", mean(training_log$best_rmse)),
  sprintf("    - Median: %.3f", median(training_log$best_rmse)),
  "",
  "  MAE:",
  sprintf("    - Mean: %.3f", mean(training_log$best_mae)),
  sprintf("    - Median: %.3f", median(training_log$best_mae)),
  "",
  "Algorithm Comparison:",
  "",
  sprintf("  Elastic Net - Average R²: %.3f", 
          mean(training_log$elastic_net_r2, na.rm = TRUE)),
  sprintf("  Random Forest - Average R²: %.3f", 
          mean(training_log$rf_r2, na.rm = TRUE)),
  sprintf("  SVM - Average R²: %.3f", 
          mean(training_log$svm_r2, na.rm = TRUE)),
  "",
  "Best Algorithm Frequency:",
  "",
  paste(sprintf("  %s: %d drugs", 
                names(table(training_log$best_algorithm)),
                table(training_log$best_algorithm)), 
        collapse = "\n"),
  "",
  "Top 5 Best Performing Drugs:",
  "",
  paste(sprintf("  %d. %s: R² = %.3f (%s)", 
                1:5,
                head(training_log[order(-training_log$best_rsquared), "drug"], 5),
                head(training_log[order(-training_log$best_rsquared), "best_rsquared"], 5),
                head(training_log[order(-training_log$best_rsquared), "best_algorithm"], 5)),
        collapse = "\n"),
  "",
  "Output Files:",
  "",
  "  Performance data:",
  "    - training_log.rds / .csv",
  "    - all_models_results.rds",
  "",
  "  Visualizations:",
  "    - 01_rsquared_distribution.png",
  "    - 02_algorithm_comparison.png",
  "    - 03_performance_vs_features.png",
  "",
  sprintf("  Trained models: %d files in models/trained/", nrow(training_log)),
  "",
  "=============================================================================",
  "NEXT STEPS:",
  "  1. Review performance visualizations",
  "  2. Analyze feature importance for top models",
  "  3. Run 06_predict.R to evaluate on test set",
  "  4. Generate final report and visualizations",
  "============================================================================="
)

writeLines(report, "reports/training_report.txt")
cat("  ✓ Report saved to reports/training_report.txt\n\n")

# =============================================================================
# SUMMARY
# =============================================================================

cat("=============================================================================\n")
cat("TRAINING SUMMARY\n")
cat("=============================================================================\n\n")

cat(sprintf("Successfully trained: %d/%d drugs\n", 
            nrow(training_log), length(drugs_to_train)))
cat(sprintf("Average R²: %.3f\n", mean(training_log$best_rsquared)))
cat(sprintf("Best performing drug: %s (R² = %.3f)\n",
            training_log$drug[which.max(training_log$best_rsquared)],
            max(training_log$best_rsquared)))

cat("\nGenerated files:\n")
cat("  Performance:\n")
cat("    - training_log.rds / .csv\n")
cat("    - all_models_results.rds\n")
cat("\n  Visualizations:\n")
cat("    - 01_rsquared_distribution.png\n")
cat("    - 02_algorithm_comparison.png\n")
cat("    - 03_performance_vs_features.png\n")
cat("\n  Report:\n")
cat("    - training_report.txt\n")
cat(sprintf("\n  Models: %d files in models/trained/\n", nrow(training_log)))

cat("\nNext steps:\n")
cat("  - Review performance metrics\n")
cat("  - Run 06_predict.R for test set evaluation\n")
cat("  - Run 07_evaluate.R for detailed analysis\n\n")

cat("=============================================================================\n")
cat("Training Completed Succesfully!\n")
cat("=============================================================================\n")

# =============================================================================
# End of Script
# =============================================================================