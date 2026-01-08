#  RxPression

![RxPression](RxPression/results/figures/RxPression.png)

![Status](https://img.shields.io/badge/status-complete-success)
[![R](https://img.shields.io/badge/R-4.0+-blue.svg)](https://www.r-project.org/)
![Drugs](https://img.shields.io/badge/drugs-173-orange)
![Test R^2](https://img.shields.io/badge/test%20R%C2%B2-0.172-brightgreen)

##  Overview

Machine learning pipeline untuk memprediksi respons sel kanker terhadap obat kemoterapi berdasarkan profil ekspresi gen. Project ini menggunakan data dari **GDSC (Genomics of Drug Sensitivity in Cancer)** untuk membangun model prediktif yang dapat membantu personalized cancer therapy.

##  Tujuan Project
- Memprediksi IC50 (drug sensitivity) dari profil ekspresi gen
- Identifikasi biomarker genetik untuk drug response
- Perbandingan berbagai algoritma machine learning
- Interpretasi model untuk clinical insight

##  Dataset
**Source:** GDSC (Genomics of Drug Sensitivity in Cancer)
- **Gene Expression:** RNA-seq dari ~1000 cell lines
- **Drug Response:** IC50 values untuk 100+ anticancer drugs
- **Features:** ~17,000 genes
- **Samples:** ~700 cell lines dengan complete data


##  Technologies

- **Language:** R (4.0+)
- **Key Libraries:**
  - `PharmacoGx` - Pharmacogenomics data handling
  - `caret` - Machine learning framework
  - `xgboost` - Gradient boosting
  - `randomForest` - Random forest models
  - `glmnet` - Elastic net regression
  - `ggplot2` - Visualization
  - `tidyverse` - Data manipulation

##  Installation

```r
# Install required packages
install.packages(c("tidyverse", "caret", "xgboost", "randomForest", 
                   "glmnet", "ggplot2", "pheatmap", "corrplot"))

# Install Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("PharmacoGx", "Biobase"))
```

## 🚀 Quick Start
```r
# 1. Download data
source("scripts/01_data_download.R")

# 2. Preprocess data
source("scripts/02_data_preprocessing.R")

# 3. Feature selection
source("scripts/03_feature_selection.R")

# 4. Train models
source("scripts/04_model_training.R")

# 5. Evaluate models
source("scripts/05_model_evaluation.R")

# 6. Make predictions
source("scripts/06_prediction.R")
```

##  Workflow
```
Raw Data → Preprocessing → Feature Selection → Model Training → Evaluation → Prediction
   ↓            ↓               ↓                   ↓              ↓           ↓
 GDSC      Normalization   Variance Filter    Random Forest    R² Score    New Samples
           Missing Data    Correlation         XGBoost          RMSE
           Log Transform   PCA/LASSO           Elastic Net      MAE
```

### Model Performance (Expected)

| Model | R² | RMSE | MAE |
|-------|-----|------|-----|
| Random Forest | 0.65-0.75 | TBD | TBD |
| XGBoost | 0.70-0.80 | TBD | TBD |
| Elastic Net | 0.55-0.65 | TBD | TBD |

### Top Predictive Genes
(Will be updated after analysis)

## 📁 Project Structure
```
drug-response-prediction/
├── data/              # Raw and processed data
├── scripts/           # Analysis scripts
├── models/            # Trained models
├── results/           # Figures and tables
├── notebooks/         # Exploratory analysis
└── utils/             # Helper functions
```

##  Key Features

1. **Automated Pipeline:** End-to-end dari download hingga prediction
2. **Multiple Algorithms:** RF, XGBoost, Elastic Net, ensemble
3. **Feature Selection:** Variance filtering, correlation analysis, LASSO
4. **Cross-Validation:** K-fold CV untuk robust evaluation
5. **Visualization:** Heatmaps, feature importance, prediction plots
6. **Interpretability:** SHAP values, partial dependence plots

##  Usage Example
```r
# Load the trained model
model <- readRDS("models/saved_models/xgboost_model.rds")

# Prepare new gene expression data
new_data <- read.csv("data/new_samples.csv")

# Make predictions
predictions <- predict(model, newdata = new_data)

# IC50 values (lower = more sensitive)
print(predictions)
```

##  Contributing

Contributions are welcome! Please feel free to submit a Pull Request.
=======

##  References

1. Yang, W. et al. (2013). "Genomics of Drug Sensitivity in Cancer (GDSC)". *Nucleic Acids Research*.
2. Garnett, M.J. et al. (2012). "Systematic identification of genomic markers of drug sensitivity in cancer cells". *Nature*.

## 👤 Author
[Farel Immanuel]
- GitHub: [@prbfarel](https://github.com/prbfarel)
- Email: farelpurba09@gmail.com

##  Acknowledgments
- GDSC consortium for providing the data
- Bioconductor community for PharmacoGx package
---
