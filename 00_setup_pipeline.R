# ============================================================================
# TCGA Transcriptomics Pipeline - Setup & Environment Configuration
# ============================================================================
# This script sets up the project environment, installs required packages,
# and creates utility functions for the entire pipeline
# ============================================================================

# Set working directory
setwd("~/workspace/rstudio/tcga_explorer")

# Create project directory structure
create_project_structure <- function() {
  dirs <- c(
    "data/raw",
    "data/processed",
    "data/clinical",
    "results/de_analysis",
    "results/cell_composition",
    "results/survival_models",
    "results/predictions",
    "logs",
    "scripts"
  )
  
  for (dir in dirs) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE, showWarnings = FALSE)
      cat(sprintf("✓ Created directory: %s\n", dir))
    }
  }
}

# ============================================================================
# Install and Load Required Packages
# ============================================================================

install_and_load_packages <- function() {
  
  # CRAN packages
  cran_packages <- c(
    "tidyverse",
    "data.table",
    "parallel",
    "doParallel",
    "foreach",
    "glue",
    "logger",
    "survival",
    "survminer",
    "randomForest",
    "caret",
    "pROC",
    "ggplot2",
    "ComplexHeatmap",
    "circlize",
    "igraph"
  )
  
  # Bioconductor packages
  bioc_packages <- c(
    "UCSCXenaTools",
    "DESeq2",
    "tidybulk",
    "sva",
    "limma",
    "edgeR",
    "immunedeconv",
    "SingleCellExperiment",
    "SummarizedExperiment"
  )
  
  # Install CRAN packages
  cat("\n=== Installing CRAN Packages ===\n")
  for (pkg in cran_packages) {
    if (!require(pkg, character.only = TRUE)) {
      cat(sprintf("Installing %s...\n", pkg))
      install.packages(pkg, dependencies = TRUE, quiet = TRUE)
      library(pkg, character.only = TRUE)
    } else {
      cat(sprintf("✓ %s already installed\n", pkg))
    }
  }
  
  # Install Bioconductor packages
  cat("\n=== Installing Bioconductor Packages ===\n")
  
  # Check if BiocManager is installed
  if (!require("BiocManager", character.only = TRUE)) {
    install.packages("BiocManager", quiet = TRUE)
    library("BiocManager")
  }
  
  for (pkg in bioc_packages) {
    if (!require(pkg, character.only = TRUE)) {
      cat(sprintf("Installing %s from Bioconductor...\n", pkg))
      BiocManager::install(pkg, ask = FALSE, update = FALSE, quiet = TRUE)
      library(pkg, character.only = TRUE)
    } else {
      cat(sprintf("✓ %s already installed\n", pkg))
    }
  }
  
  cat("\n✓ All packages loaded successfully!\n")
}

# ============================================================================
# Setup Logging Framework
# ============================================================================

setup_logging <- function() {
  library(logger)
  
  # Create log file
  log_file <- "logs/pipeline_log.txt"
  
  # Configure logger
  log_appender(appender_tee(log_file))
  log_threshold(DEBUG)
  
  log_info("=== TCGA Transcriptomics Pipeline Started ===")
  log_info("Timestamp: {Sys.time()}")
  log_info("Working Directory: {getwd()}")
  
  return(log_file)
}

# ============================================================================
# Utility Functions
# ============================================================================

# Function to check internet connection
check_internet <- function() {
  tryCatch({
    readLines("http://www.google.com", n = 1)
    return(TRUE)
  }, error = function(e) {
    return(FALSE)
  })
}

# Function to safely download data with retry logic
safe_download <- function(url, destfile, max_retries = 3) {
  library(logger)
  
  for (attempt in 1:max_retries) {
    tryCatch({
      log_info("Downloading from: {url} (Attempt {attempt}/{max_retries})")
      download.file(url, destfile, mode = "wb", quiet = TRUE)
      log_info("✓ Successfully downloaded to: {destfile}")
      return(TRUE)
    }, error = function(e) {
      log_warn("Download attempt {attempt} failed: {e$message}")
      if (attempt < max_retries) {
        Sys.sleep(5)  # Wait 5 seconds before retry
      }
    })
  }
  
  log_error("Failed to download after {max_retries} attempts: {url}")
  return(FALSE)
}

# Function to validate data integrity
validate_data <- function(data, name, min_rows = 100, min_cols = 10) {
  library(logger)
  
  if (is.null(data) || nrow(data) == 0) {
    log_error("Data validation failed for {name}: Empty or NULL data")
    return(FALSE)
  }
  
  if (nrow(data) < min_rows) {
    log_warn("Data {name} has fewer rows than expected: {nrow(data)} < {min_rows}")
  }
  
  if (ncol(data) < min_cols) {
    log_warn("Data {name} has fewer columns than expected: {ncol(data)} < {min_cols}")
  }
  
  log_info("✓ Data validation passed for {name}: {nrow(data)} rows × {ncol(data)} columns")
  return(TRUE)
}

# Function to save data with metadata
save_data <- function(data, filename, description = "") {
  library(logger)
  
  filepath <- file.path("data/processed", filename)
  
  tryCatch({
    if (grepl("\\.rds$", filename)) {
      saveRDS(data, filepath)
    } else if (grepl("\\.csv$", filename)) {
      write.csv(data, filepath, row.names = TRUE)
    } else if (grepl("\\.tsv$", filename)) {
      write.table(data, filepath, sep = "\t", row.names = TRUE)
    }
    
    log_info("✓ Saved {filename} ({description}): {nrow(data)} rows × {ncol(data)} columns")
    return(filepath)
  }, error = function(e) {
    log_error("Failed to save {filename}: {e$message}")
    return(NULL)
  })
}

# Function to load data with validation
load_data <- function(filename) {
  library(logger)
  
  filepath <- file.path("data/processed", filename)
  
  if (!file.exists(filepath)) {
    log_error("File not found: {filepath}")
    return(NULL)
  }
  
  tryCatch({
    if (grepl("\\.rds$", filename)) {
      data <- readRDS(filepath)
    } else if (grepl("\\.csv$", filename)) {
      data <- read.csv(filepath, row.names = 1)
    } else if (grepl("\\.tsv$", filename)) {
      data <- read.table(filepath, sep = "\t", row.names = 1)
    }
    
    log_info("✓ Loaded {filename}: {nrow(data)} rows × {ncol(data)} columns")
    return(data)
  }, error = function(e) {
    log_error("Failed to load {filename}: {e$message}")
    return(NULL)
  })
}

# ============================================================================
# Main Execution
# ============================================================================

main <- function() {
  cat("\n╔════════════════════════════════════════════════════════════════╗\n")
  cat("║  TCGA Transcriptomics Pipeline - Environment Setup             ║\n")
  cat("╚════════════════════════════════════════════════════════════════╝\n\n")
  
  # Create project structure
  cat("Creating project directory structure...\n")
  create_project_structure()
  
  # Setup logging
  cat("\nSetting up logging framework...\n")
  log_file <- setup_logging()
  
  # Install and load packages
  cat("\nInstalling and loading required packages...\n")
  cat("(This may take several minutes on first run)\n")
  install_and_load_packages()
  
  # Check internet connection
  cat("\nChecking internet connection...\n")
  if (check_internet()) {
    cat("✓ Internet connection available\n")
  } else {
    cat("⚠ Warning: No internet connection detected\n")
  }
  
  cat("\n╔════════════════════════════════════════════════════════════════╗\n")
  cat("║  ✓ Environment Setup Complete!                                ║\n")
  cat("║  Log file: logs/pipeline_log.txt                              ║\n")
  cat("╚════════════════════════════════════════════════════════════════╝\n\n")
  
  # Export utility functions to global environment
  list(
    create_project_structure = create_project_structure,
    install_and_load_packages = install_and_load_packages,
    setup_logging = setup_logging,
    check_internet = check_internet,
    safe_download = safe_download,
    validate_data = validate_data,
    save_data = save_data,
    load_data = load_data
  )
}

# Run setup
if (!interactive()) {
  main()
} else {
  utils <- main()
  cat("\nUtility functions available in 'utils' object\n")
}
