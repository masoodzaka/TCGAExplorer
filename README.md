# TCGA Transcriptomics Pipeline

A comprehensive R-based bioinformatics pipeline for analyzing TCGA transcript-level data with GTEx normal tissue controls, performing differential expression analysis, cell type composition analysis, and building patient outcome prediction models.

## 📋 Pipeline Overview

This pipeline implements a complete workflow for TCGA transcriptomics analysis:

```
Data Acquisition
    ↓
Preprocessing (ComBat + DESeq2)
    ↓
Differential Expression Analysis (tidybulk)
    ↓
Cell Type Composition (CIBERSORT via tidybulk)
    ↓
Clinical Data Integration
    ↓
Survival Analysis & Prediction Models
    ├── Cox Proportional Hazards
    ├── Random Forest
    └── Logistic Regression
```

## 🎯 Key Features

### 1. **Multi-Indication Analysis**
- Analyzes ALL available TCGA cancer indications
- Compares cancer samples vs GTEx normal tissue
- Separate analysis for each indication

### 2. **Advanced Preprocessing**
- **Batch Effect Correction**: ComBat algorithm
- **Normalization**: DESeq2 with variance stabilization transformation (VST)
- **Filtering**: Removes low-count transcripts
- **Quality Control**: Comprehensive validation at each step

### 3. **Differential Expression Analysis**
- Uses `tidybulk` package for streamlined analysis
- Identifies significantly dysregulated transcripts
- Generates volcano plots and summary statistics
- Extracts top 100 significant transcripts per indication

### 4. **Cell Type Composition Analysis**
- CIBERSORT deconvolution via `tidybulk`
- Separate analysis for tumor and normal samples
- Generates heatmaps and boxplots
- Quantifies immune cell infiltration

### 5. **Clinical Integration**
- Merges transcriptomic data with clinical metadata
- Extracts survival information
- Prepares feature matrices for prediction
- 70/30 train/test split (stratified)

### 6. **Survival Prediction Models**
- **Cox Proportional Hazards**: C-index evaluation
- **Random Forest**: AUC and accuracy metrics
- **Logistic Regression**: Binary outcome prediction
- ROC curves and performance comparison

## 📁 Project Structure

```
tcga_pipeline/
├── 00_setup.R                      # Environment setup & utilities
├── 01_data_acquisition.R           # Download TCGA & GTEx data
├── 02_preprocessing.R              # ComBat + DESeq2 normalization
├── 03_differential_expression.R    # DE analysis with tidybulk
├── 04_cell_composition.R           # CIBERSORT deconvolution
├── 05_clinical_integration.R       # Clinical data merging
├── 06_survival_prediction.R        # Model building & evaluation
├── master_pipeline.R               # Orchestrates all stages
├── README.md                       # This file
│
├── data/
│   ├── raw/                        # Downloaded raw data
│   ├── processed/                  # Preprocessed data
│   └── clinical/                   # Clinical metadata
│
├── results/
│   ├── de_analysis/                # DE results & volcano plots
│   ├── cell_composition/           # Cell composition heatmaps
│   ├── survival_models/            # Model predictions & ROC curves
│   └── predictions/                # Final predictions
│
└── logs/
    └── pipeline_log.txt            # Detailed execution log
```

## 🚀 Quick Start

### Prerequisites

```r
# Install required packages (run once)
source("00_setup.R")
```

### Run Complete Pipeline

```r
# Execute the entire pipeline
source("master_pipeline.R")
main()
```

### Run Individual Stages

```r
# Stage 1: Data Acquisition
source("01_data_acquisition.R")
data_results <- main()

# Stage 2: Preprocessing
source("02_preprocessing.R")
preprocessed_data <- main()

# Stage 3: Differential Expression
source("03_differential_expression.R")
de_results <- main()

# Stage 4: Cell Composition
source("04_cell_composition.R")
cell_comp_results <- main()

# Stage 5: Clinical Integration
source("05_clinical_integration.R")
clinical_results <- main()

# Stage 6: Survival Prediction
source("06_survival_prediction.R")
survival_results <- main()
```

## 📊 Input Data

### TCGA Data
- **Source**: UCSC Xena Hub (tcgaHub)
- **Data Type**: Transcript-level expression (RSEM)
- **Indications**: All available cancer types
- **Format**: Automatically downloaded via UCSCXenaTools

### GTEx Data
- **Source**: UCSC Xena Hub (gtexHub)
- **Data Type**: Transcript-level expression
- **Tissues**: All available tissues
- **Format**: Automatically downloaded via UCSCXenaTools

### Clinical Data
- **Source**: TCGA clinical metadata
- **Variables**: Survival time, event status, demographics, pathology
- **Format**: TSV/CSV files

## 📈 Output Files

### Preprocessing Results
- `{indication}_raw_counts.rds` - Raw count matrix
- `{indication}_normalized_counts.rds` - VST-normalized counts
- `{indication}_metadata.rds` - Sample metadata
- `{indication}_dds.rds` - DESeq2 object
- `preprocessing_summary.csv` - Summary statistics

### Differential Expression Results
- `{indication}_de_results.rds` - Full DE results
- `{indication}_de_results.csv` - DE results (CSV)
- `{indication}_volcano.pdf` - Volcano plot
- `{indication}_top_transcripts.txt` - Top 100 significant transcripts
- `de_summary.csv` - DE summary table

### Cell Composition Results
- `{indication}_tumor_composition.rds` - Tumor cell fractions
- `{indication}_normal_composition.rds` - Normal cell fractions
- `{indication}_tumor_composition.pdf` - Tumor heatmap
- `{indication}_normal_composition.pdf` - Normal heatmap
- `{indication}_composition_boxplot.pdf` - Boxplot comparison

### Clinical Integration Results
- `{indication}_feature_matrix.rds` - Feature matrix
- `{indication}_train_test_splits.rds` - Train/test splits
- `clinical_integration_summary.csv` - Integration summary

### Survival Model Results
- `cox_models.rds` - Cox model objects
- `rf_models.rds` - Random Forest models
- `logistic_models.rds` - Logistic regression models
- `cox_predictions.rds` - Cox predictions
- `rf_predictions.rds` - RF predictions
- `logistic_predictions.rds` - Logistic predictions
- `{indication}_roc_curves.pdf` - ROC curves
- `model_performance_summary.csv` - Performance metrics

## 🔧 Configuration

Edit `master_pipeline.R` to customize pipeline execution:

```r
PIPELINE_CONFIG <- list(
  download_data = TRUE,              # Download fresh data
  apply_combat = TRUE,               # Apply ComBat correction
  deseq2_normalize = TRUE,           # DESeq2 normalization
  perform_de_analysis = TRUE,        # DE analysis
  perform_cell_composition = TRUE,   # Cell composition
  build_survival_models = TRUE,      # Survival models
  output_dir = "results",
  log_file = "logs/pipeline_log.txt"
)
```

## 📊 Analysis Parameters

### Preprocessing
- **Batch Effect Correction**: ComBat (parametric prior)
- **Normalization**: DESeq2 VST
- **Filtering**: CPM > 1 in ≥10% of samples
- **Transformation**: Variance stabilization

### Differential Expression
- **Method**: edgeR (via tidybulk)
- **Comparison**: Cancer vs Normal tissue
- **Significance Threshold**: FDR < 0.05, |log2FC| > 1
- **Top Transcripts**: 100 per indication

### Cell Type Composition
- **Method**: CIBERSORT (quantiseq alternative)
- **Signature**: Immune cell reference
- **Analysis**: Separate tumor and normal

### Survival Models
- **Train/Test Split**: 70/30 (stratified)
- **Cox Model**: C-index evaluation
- **Random Forest**: 500 trees, AUC metric
- **Logistic Regression**: Binary outcome (median survival)

## 📝 Logging

All pipeline execution is logged to `logs/pipeline_log.txt`:

```
[2024-01-15 10:30:45] INFO: TCGA Transcriptomics Pipeline - Master Execution
[2024-01-15 10:30:45] INFO: Timestamp: 2024-01-15 10:30:45
[2024-01-15 10:30:46] INFO: Stage: Data Acquisition
...
```

## 🐛 Troubleshooting

### Issue: Package Installation Fails
```r
# Install packages manually
install.packages("package_name")
BiocManager::install("bioconductor_package")
```

### Issue: Data Download Fails
- Check internet connection
- Verify Xena Hub is accessible
- Check available disk space

### Issue: Memory Issues
- Process indications one at a time
- Reduce number of features
- Use smaller train/test split

### Issue: Missing Survival Data
- Check clinical data format
- Verify column names
- Ensure survival columns exist

## 📚 References

### Key Packages
- **tidybulk**: https://github.com/tidyomics/tidybulk
- **DESeq2**: https://bioconductor.org/packages/DESeq2/
- **UCSCXenaTools**: https://github.com/ropensci/UCSCXenaTools
- **sva**: https://bioconductor.org/packages/sva/

### Methods
- ComBat: Johnson et al. (2007) Biostatistics
- DESeq2: Love et al. (2014) Genome Biology
- CIBERSORT: Newman et al. (2015) Nature Methods
- Random Forest: Breiman (2001) Machine Learning

## 📄 License

This pipeline is provided as-is for research purposes.

## 👤 Author

Created for comprehensive TCGA transcriptomics analysis with clinical integration and outcome prediction.

## 🤝 Contributing

For improvements or bug reports, please document:
1. Error message and traceback
2. Which stage failed
3. System information (R version, OS)
4. Reproducible example

## ⚠️ Important Notes

1. **Data Size**: Full TCGA + GTEx data can be >100GB. Ensure sufficient disk space.
2. **Computation Time**: Complete pipeline may take 24-48 hours depending on hardware.
3. **Internet Connection**: Required for initial data download.
4. **Memory Requirements**: Minimum 16GB RAM recommended.
5. **Clinical Data**: Ensure clinical files are in correct format (TSV/CSV).

## 🎓 Citation

If you use this pipeline, please cite:
- tidybulk: Mangiola et al. (2021) F1000Research
- DESeq2: Love et al. (2014) Genome Biology
- UCSCXenaTools: Liu et al. (2019) Bioinformatics

---

**Last Updated**: January 2024
**Version**: 1.0
**Status**: Production Ready
