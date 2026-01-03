#Drug Response Prediction - Script 01: Download GDSC Data

#Purpose: Download gene expression and drug response data from GDSC
#Author: prbfarel
#Date: 2025-12-20

# ============================================

cat("\n Starting GDSC Data Download \n")

# Load required libraries
suppressPackageStartupMessages({
  library(PharmacoGx)
  library(SummarizedExperiment)
  library(tidyverse)
  library(reshape2)
  library(ggplot2)
  library(pheatmap)
})

# Check versions
cat(" Package versions:\n")
cat("   PharmacoGx:", as.character(packageVersion("PharmacoGx")), "\n")
cat("   SummarizedExperiment:", as.character(packageVersion("SummarizedExperiment")), "\n\n")

# Set output directory
output_dir <- "data/raw"
if(!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ============================================

# 1. Download GDSC Dataset
# =================================


cat(" Step 1: Downloading GDSC dataset...\n")
cat(" This may take 5-10 minutes depending on internet connection...")

# Define path to GDSC file
gdsc_file <- file.path(output_dir, "GDSC_2020(v2-8.2).rds")

# Check if already downloaded
if(file.exists(gdsc_file)) {
  cat("📂 Found existing GDSC file, loading from disk...\n")
  gdsc <- readRDS(gdsc_file)
  cat("✅ GDSC loaded successfully!\n\n")
} else {
  # Try download
  tryCatch({
    # List available PSets first
    cat(" Checking available datasets...\n")
    available <- availablePSets()
    gdsc_versions <- available[grep("GDSC", available$PSet.Name, ignore.case = TRUE), ]
    
    cat("Available GDSC versions:\n")
    print(gdsc_versions[, c("PSet.Name", "Date.Updated")])
    
    # Download latest version
    pset_name <- gdsc_versions$PSet.Name[1]
    cat("\nDownloading:", pset_name, "\n")
    
    gdsc <- downloadPSet(pset_name, saveDir = output_dir)
    
    # Save immediately
    saveRDS(gdsc, gdsc_file)
    cat("✅ GDSC downloaded and saved!\n\n")
    
  }, error = function(e) {
    cat("\n Auto-download failed:", e$message, "\n\n")
    cat("=" , rep("=", 60), "\n", sep = "")
    cat("️  MANUAL DOWNLOAD INSTRUCTIONS\n")
    cat("=" , rep("=", 60), "\n", sep = "")
    cat("\n1. Visit: https://orcestra.ca/pset/GDSC\n")
    cat("2. Download: GDSC_2020(v2-8.2).rds (~500 MB)\n")
    cat("3. Save to:", normalizePath(output_dir, mustWork = FALSE), "\n")
    cat("4. Rename to: gdsc_full.rds\n")
    cat("5. Re-run this script\n\n")
    
    stop("Please download GDSC manually and re-run.")
  })
}

# Verify GDSC object
cat(" GDSC Object Info:\n")
cat("   Class:", class(gdsc), "\n")
cat("   Cell lines:", length(cellNames(gdsc)), "\n")
cat("   Drugs:", length(drugNames(gdsc)), "\n\n")

# =================================
# 2. Extract Gene Expression Data
# =================================

cat("🧬 Step 2: Extracting gene expression data...\n")

# Get available molecular data types
mol_types <- mDataNames(gdsc)
cat("   Available molecular types:", mol_types, "\n")

# Extract RNA data - this returns a matrix directly
mol_profiles <- molecularProfiles(gdsc, mDataType = "rna")
cat("   Molecular profile class:", class(mol_profiles)[1], "\n")

# Since mol_profiles is already a matrix, just transpose it
# Current: genes (rows) × samples (cols)
# We want: samples (rows) × genes (cols)
expr_matrix <- as.data.frame(t(mol_profiles))

cat("   ✅ Expression matrix extracted\n")
cat("   - Samples:", nrow(expr_matrix), "\n")
cat("   - Genes:", ncol(expr_matrix), "\n")
cat("   - Missing values:", sum(is.na(expr_matrix)), 
    "(", round(mean(is.na(expr_matrix))*100, 2), "%)\n\n")


# =================================
# 3. Extract Drug Response Data
# =================================


cat("💊 Step 3: Extracting drug response data...\n")

# Get sensitivity info
sens_info <- sensitivityInfo(gdsc)
cat("   Total sensitivity measurements:", nrow(sens_info), "\n")

# Get sensitivity profiles
sens_profiles <- sensitivityProfiles(gdsc)
cat("   Sensitivity measures available:", 
    paste(colnames(sens_profiles), collapse = ", "), "\n")

# Choose measure (aac_recomputed is recommended)
if("aac_recomputed" %in% colnames(sens_profiles)) {
  measure <- "aac_recomputed"
} else if("ic50_recomputed" %in% colnames(sens_profiles)) {
  measure <- "ic50_recomputed"
} else {
  measure <- colnames(sens_profiles)[1]
}

cat("   Using measure:", measure, "\n")

# Create drug response dataframe using correct column names
drug_response_long <- data.frame(
  cellid = sens_info$sampleid,        # Use 'sampleid' column
  drugid = sens_info$treatmentid,     # Use 'treatmentid' column
  value = sens_profiles[, measure],
  stringsAsFactors = FALSE
)

cat("   - Total observations:", nrow(drug_response_long), "\n")

# Remove NA values
drug_response_long <- drug_response_long[!is.na(drug_response_long$value), ]
cat("   - After removing NA:", nrow(drug_response_long), "\n")

# Reshape to wide format (cell lines × drugs)
ic50_matrix <- reshape2::dcast(drug_response_long, 
                               cellid ~ drugid, 
                               value.var = "value",
                               fun.aggregate = mean)

rownames(ic50_matrix) <- ic50_matrix$cellid
ic50_matrix$cellid <- NULL

cat("   - Cell lines:", nrow(ic50_matrix), "\n")
cat("   - Drugs:", ncol(ic50_matrix), "\n")
cat("   - Missing values:", sum(is.na(ic50_matrix)), 
    "(", round(mean(is.na(ic50_matrix))*100, 2), "%)\n\n")


# =================================
# 4. Get Drug Information
# =================================


cat(" Step 4: Extracting drug information...\n")

drug_info <- drugInfo(gdsc)
cat("   - Total drugs in database:", nrow(drug_info), "\n")
cat("   - Columns:", ncol(drug_info), "\n\n")



# =================================
# 5. Get Cell Line Information
# =================================


cat(" Step 5: Extracting cell line information...\n")

cell_info <- cellInfo(gdsc)
cat("   - Total cell lines:", nrow(cell_info), "\n")
cat("   - Tissue types:", length(unique(cell_info$tissueid)), "\n\n")



# =================================
# 6. Processed Data (Align Data)
# =================================


cat(" Step 6: Aligning gene expression and drug response data...\n")

# Find common cell lines
common_cells <- intersect(rownames(expr_matrix), rownames(ic50_matrix))
cat("   - Cell lines in expression data:", nrow(expr_matrix), "\n")
cat("   - Cell lines in drug response data:", nrow(ic50_matrix), "\n")
cat("   - Common cell lines:", length(common_cells), "\n")

# Subset to common cell lines
expr_matrix_aligned <- expr_matrix[common_cells, ]
ic50_matrix_aligned <- ic50_matrix[common_cells, ]

cat("   ✅ Data aligned!\n\n")


# =================================
# 7. Save Processed Data
# =================================

cat(" Step 7: Saving data to disk...\n")

# Save full objects
saveRDS(gdsc, file.path(output_dir, "gdsc_full.rds"))
saveRDS(expr_matrix, file.path(output_dir, "gene_expression_raw.rds"))
saveRDS(ic50_matrix, file.path(output_dir, "drug_response_raw.rds"))

# Save aligned data
saveRDS(expr_matrix_aligned, file.path(output_dir, "gene_expression.rds"))
saveRDS(ic50_matrix_aligned, file.path(output_dir, "drug_response.rds"))

# Save metadata
saveRDS(drug_info, file.path(output_dir, "drug_info.rds"))
saveRDS(cell_info, file.path(output_dir, "cell_info.rds"))

# Save as CSV
write.csv(drug_info, file.path(output_dir, "drug_info.csv"), row.names = TRUE)
write.csv(cell_info, file.path(output_dir, "cell_info.csv"), row.names = TRUE)

# Save sample of expression data (first 100 genes for quick viewing)
write.csv(expr_matrix_aligned[1:min(10, nrow(expr_matrix_aligned)), 
                              1:min(100, ncol(expr_matrix_aligned))], 
          file.path(output_dir, "gene_expression_sample.csv"))

cat("✅ All data saved!\n\n")

# =================================
# 8. Data Summary
# =================================

cat("\n")
cat("=" , rep("=", 60), "\n", sep = "")
cat("Data Summary\n")
cat("=" , rep("=", 60), "\n", sep = "")

cat("\n Gene Expression Matrix:\n")
cat("   Dimensions:", nrow(expr_matrix_aligned), "samples ×", 
    ncol(expr_matrix_aligned), "genes\n")
cat("   Value range:", round(range(expr_matrix_aligned, na.rm = TRUE), 2), "\n")

cat("\n Drug Response Matrix:\n")
cat("   Dimensions:", nrow(ic50_matrix_aligned), "cell lines ×", 
    ncol(ic50_matrix_aligned), "drugs\n")
cat("   Value range:", round(range(ic50_matrix_aligned, na.rm = TRUE), 2), "\n")
cat("   Measure:", measure, "\n")

cat("\n Cell Lines by Tissue (Top 10):\n")
tissue_counts <- sort(table(cell_info$tissueid), decreasing = TRUE)
print(head(tissue_counts, 10))

cat("\n Sample Drugs:\n")
sample_drugs <- head(drug_info[, c("DRUG_ID", "DRUG_NAME", "TARGET_PATHWAY")], 10)
print(sample_drugs)

# =================================
# 9. Quick Visualization
# =================================


cat(" Creating quick visualizations...\n")

library(ggplot2)

# Create results/figures if not exists
if(!dir.exists("results/figures")) {
  dir.create("results/figures", recursive = TRUE)
}

# Plot 1: Distribution of drug responses
first_drug <- colnames(ic50_matrix_aligned)[1]
p1 <- ggplot(data.frame(Response = ic50_matrix_aligned[, first_drug]), 
             aes(x = Response)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "black", alpha = 0.7) +
  labs(title = paste("Drug Response Distribution:", first_drug),
       x = paste(measure, "value"), 
       y = "Frequency") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("results/figures/drug_response_distribution.png", 
       plot = p1, width = 8, height = 6, dpi = 300)

# Plot 2: Tissue distribution
tissue_df <- data.frame(Tissue = names(tissue_counts), 
                        Count = as.numeric(tissue_counts))
tissue_df <- tissue_df[1:min(15, nrow(tissue_df)), ]  # Top 15

p2 <- ggplot(tissue_df, aes(x = reorder(Tissue, Count), y = Count)) +
  geom_bar(stat = "identity", fill = "coral", alpha = 0.7) +
  coord_flip() +
  labs(title = "Cell Lines by Tissue Type",
       x = "Tissue Type", y = "Number of Cell Lines") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("results/figures/tissue_distribution.png", 
       plot = p2, width = 10, height = 6, dpi = 300)

# Plot 3: Gene expression heatmap (sample)
library(pheatmap)
set.seed(123)

# Sample genes and cells
sample_genes <- sample(1:ncol(expr_matrix_aligned), min(50, ncol(expr_matrix_aligned)))
sample_cells <- sample(1:nrow(expr_matrix_aligned), min(50, nrow(expr_matrix_aligned)))

# Extract subset
expr_subset <- expr_matrix_aligned[sample_cells, sample_genes]

# Remove genes with NA or Inf values
valid_genes <- apply(expr_subset, 2, function(x) all(is.finite(x)))
expr_subset <- expr_subset[, valid_genes]

# Remove cells with NA or Inf values
valid_cells <- apply(expr_subset, 1, function(x) all(is.finite(x)))
expr_subset <- expr_subset[valid_cells, ]

# Check if we have enough data
if(nrow(expr_subset) > 2 && ncol(expr_subset) > 2) {
  png("results/figures/expression_heatmap_sample.png", 
      width = 10, height = 8, units = "in", res = 300)
  
  pheatmap(t(expr_subset),
           show_rownames = FALSE,
           show_colnames = FALSE,
           main = "Gene Expression Heatmap (Sample)",
           color = colorRampPalette(c("blue", "white", "red"))(100),
           na_col = "grey")
  
  dev.off()
  cat("   ✅ Heatmap created\n")
} else {
  cat("   ⚠️  Not enough valid data for heatmap, skipping...\n")
}

cat("   ✅ Visualizations saved to results/figures/\n\n")

# =================================
# 10. Create Data Dictionary
# =================================

data_dict <- data.frame(
  File = c(
    "gdsc_full.rds",
    "gene_expression.rds",
    "drug_response.rds",
    "drug_info.rds",
    "cell_info.rds"
  ),
  Description = c(
    "Complete GDSC PharmacoSet object",
    paste("Aligned gene expression matrix -", nrow(expr_matrix_aligned), 
          "samples ×", ncol(expr_matrix_aligned), "genes"),
    paste("Aligned drug response matrix -", nrow(ic50_matrix_aligned), 
          "cell lines ×", ncol(ic50_matrix_aligned), "drugs"),
    "Drug annotations and information",
    "Cell line annotations and tissue types"
  ),
  Rows = c(
    "PharmacoSet",
    nrow(expr_matrix_aligned),
    nrow(ic50_matrix_aligned),
    nrow(drug_info),
    nrow(cell_info)
  ),
  Columns = c(
    "Multiple slots",
    ncol(expr_matrix_aligned),
    ncol(ic50_matrix_aligned),
    ncol(drug_info),
    ncol(cell_info)
  )
)

write.csv(data_dict, "data/raw/DATA_DICTIONARY.csv", row.names = FALSE)

cat("\n")
cat("=" , rep("=", 60), "\n", sep = "")
cat(" Download And Extraction Complete!\n")
cat("=" , rep("=", 60), "\n", sep = "")

cat("\n Files created in data/raw/:\n")
files_created <- list.files(output_dir, pattern = "\\.(rds|csv)$")
for(f in files_created) {
  cat("   -", f, "\n")
}

cat("\n Files created in results/figures/:\n")
figs_created <- list.files("results/figures", pattern = "\\.png$")
for(f in figs_created) {
  cat("   -", f, "\n")
}

cat("\n Ready for Step 4: Exploratory Data Analysis!\n\n")