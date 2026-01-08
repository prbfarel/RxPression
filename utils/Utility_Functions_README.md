# 🔧 Utility Functions

This directory contains reusable helper functions for the GDSC analysis pipeline.

##  Files

### `gdsc_utils.R`
**Purpose:** Comprehensive utility functions for GDSC analysis

**Categories:**
1. **Data Loading** - Load expression and drug response data
2. **Quality Checks** - Missing values, outliers
3. **Feature Selection** - Variance-based, correlation
4. **Model Evaluation** - Metrics calculation
5. **Visualization** - Publication themes, plotting
6. **File Management** - Directory creation, saving
7. **Logging** - Formatted messages
8. **Transformation** - Normalization, scaling
9. **Statistics** - Permutation tests, significance

---

##  Usage

### Load utilities in your script:
```r
source("utils/gdsc_utils.R")
```

### Example usage:
```r
# Load data with automatic messages
expr <- load_expression("data/processed/expr_train.rds")

# Check data quality
missing_summary <- check_missing(expr)

# Calculate metrics
metrics <- calculate_metrics(observed, predicted)
print(paste("R² =", metrics$r2))

# Save plot
p <- ggplot(data) + geom_point() + theme_publication()
save_plot(p, "figures/my_plot.png")
```

---

##  Function Reference

### Data Loading
- `load_expression(path)` - Load expression matrix
- `load_drug_response(path)` - Load drug response data

### Quality Checks
- `check_missing(data)` - Count missing values
- `detect_outliers_iqr(x, threshold)` - Identify outliers

### Feature Selection
- `select_by_variance(expr, n)` - Select top N variable genes
- `cor_with_pvalue(x, y)` - Correlation with significance

### Evaluation
- `calculate_metrics(observed, predicted)` - Regression metrics
- `calculate_classification_metrics(obs, pred)` - Classification metrics

### Visualization
- `theme_publication()` - Clean publication theme
- `save_plot(plot, filename)` - Save with consistent settings

### File Management
- `ensure_dir(path)` - Create directory if needed
- `save_with_message(object, path)` - Save with confirmation

### Logging
- `print_header(text)` - Formatted section headers
- `print_step(num, text)` - Step messages

### Transformation
- `zscore(x)` - Z-score normalization
- `minmax_normalize(x, range)` - Min-max scaling

### Statistics
- `permutation_test(x, y, n)` - Non-parametric testing

---

##  Best Practices

**DO:**
✓ Source utilities at the start of scripts
✓ Use consistent function names
✓ Add error handling in scripts
✓ Document custom modifications

**DON'T:**
✗ Modify utility functions directly (make copies)
✗ Hard-code paths (use relative paths)
✗ Skip input validation

---

## Extending Utilities

To add new functions:

1. Add to appropriate section in `gdsc_utils.R`
2. Follow existing naming conventions
3. Include roxygen comments
4. Update this README
5. Test thoroughly

Example:
```r
#' Calculate custom metric
#' @param x Input vector
#' @return Metric value
my_metric <- function(x) {
  # Your code here
  return(result)
}
```

---

## 📝 Notes

- All functions include input validation
- Messages use `message()` for easy suppression
- Errors use `stop()` for clear debugging
- Warnings use `warning()` for non-critical issues

---