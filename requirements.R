#Drug Response Prediction - Package Requirements

#Function to install and load packages
install_and_load <- function(packages) {
  for(pkg in packages) {
    if(!require(pkg, character.only = TRUE)) {
      install.packages(pkg, dedpendencies = TRUE)
      library(pkg, character.only = TRUE)
    } else {
      library(pkg, character.only = TRUE)
    }
  }
}

# CRAN packages
cran_packages <- c(
  "tidyverse",           # Data manipulation
  "caret",               # Machine learning framework
  "xgboost",             # Gradient boosting
  "randomForest",        # Random forest
  "glmnet",              # Elastic net/LASSO
  "ggplot2",             # Visualization
  "pheatmap",            # Heatmaps
  "corrplot",            # Correlation plots
  "pROC",                # ROC curves
  "gridExtra",           # Multiple plots
  "RColorBrewer",        # Color palletes
  "scales",              # Scaling functions
  "reshape2"             # Data reshaping
)

cat("Installing CRAN packages...\n")
install_and_load(cran_packages)

#Bioconductor packages
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_packages <- c(
  "PharmacoGx",           # Pharmacogenomics data
  "Biobase"              # Bioconductor base
)

cat("Installing Bioconductor packages...\n")
for(pkg in bioc_packages) {
  if(!require(pkg, character.only = TRUE)) {
    BiocManager::install(pkg)
    library(pkg, character.only = TRUE)
  } else {
    library(pkg, character.only = TRUE)
  }
}

cat("\n✅ packages installed successfully!\n")
cat("Loaded packages:\n")
print(sessionInfo())