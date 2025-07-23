# COVID-19 scRNA-seq pseudobulk analysis in R
# Author: Jean-Marie Ravel, Institut Curie
# Focus: Construction of a Seurat object for downstream analysis

# -------------------------------
# STEP 0: Setup and Environment
# -------------------------------
rm(list=ls(all.names=TRUE)); gc(); if (dev.cur() > 1) dev.off(); cat("\014")

# Install required libraries if missing
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (pkg in c("HDF5Array", "DelayedArray", "DropletUtils", "GEOquery")) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg)
}

# Load libraries
suppressMessages({
  library(Seurat)
  library(Matrix)
  library(tidyverse)
  library(DropletUtils)
  library(data.table)
  library(GEOquery)
})

# -------------------------------
# STEP 1: Load and Prepare Data
# -------------------------------
dataset <- "GSE149689"
base_path <- "/Users/jean-marie/Documents/Privé/Curie/covid19"
data_dir <- file.path(base_path, "data", dataset)
setwd(data_dir)

# Rename files if needed
# file.rename("GSE149689_matrix.mtx.gz", "matrix.mtx.gz")
# file.rename("GSE149689_features.tsv.gz", "features.tsv.gz")
# file.rename("GSE149689_barcodes.tsv.gz", "barcodes.tsv.gz")

# Load matrix
data_10x <- Read10X(data.dir = getwd())

# Create Seurat object
seurat_obj <- CreateSeuratObject(counts = data_10x, project = dataset, min.cells = 3, min.features = 200)
seurat_obj <- UpdateSeuratObject(seurat_obj)
seurat_obj <- JoinLayers(seurat_obj)

# -------------------------------
# STEP 2: Add Sample Information
# -------------------------------
barcode_suffix <- sub(".*-(\\d+)$", "\\1", colnames(seurat_obj))

sample_mapping <- data.frame(
  suffix = as.character(1:20),
  sample_id = paste0("Sample", 1:20),
  condition = c(
    "nCoV", "nCoV",   # Sample1-2
    "Flu", "Flu",     # Sample3-4
    "Normal",         # Sample5
    "Flu", "Flu", "Flu",  # Sample6-8
    "nCoV", "nCoV", "nCoV", "nCoV",  # Sample9-12
    "Normal", "Normal",             # Sample13-14
    "nCoV", "nCoV", "nCoV", "nCoV",  # Sample15-18
    "Normal", "nCoV"                # Sample19-20
  ),
  stringsAsFactors = FALSE
)

# Merge cell barcodes with sample info
cell_metadata <- data.frame(
  cell = colnames(seurat_obj),
  suffix = barcode_suffix
) %>%
  left_join(sample_mapping, by = "suffix")

# Add to Seurat metadata
seurat_obj$sample_id  <- cell_metadata$sample_id
seurat_obj$condition  <- cell_metadata$condition

# -------------------------------
# STEP 3: Merge with GEO Metadata
# -------------------------------
gset <- getGEO(dataset, GSEMatrix = TRUE)
pheno <- pData(gset[[1]])

metadata <- pheno %>%
  dplyr::select(title, geo_accession, 
         `subject group:ch1`, 
         `subject status:ch1`, 
         `cell type:ch1`, 
         `gender:ch1`) %>%
  dplyr::mutate(
    sample_id = str_extract(title, "Sample\\s*\\d+") %>% 
      str_replace_all("\\s+", "")  # e.g., "Sample 1" → "Sample1"
  )

colnames(metadata)[3:6] <- c("group", "status", "cell_type", "gender")


# Merge with sample mapping
sample_metadata <- left_join(sample_mapping, metadata, by = "sample_id") %>%
  mutate(
    status = ifelse(is.na(status), group, status),
    status = str_replace_all(status, "COVID-19 patient", ""),
    status = str_replace_all(status, "patient", ""),
    status = str_replace_all(status, "healthy ", "")
  )

# -------------------------------
# STEP 4: Add All Metadata to Seurat Object
# -------------------------------
# Join back to each cell
cell_metadata <- left_join(cell_metadata, sample_metadata, by = c("suffix", "sample_id", "condition"))

# Ensure correct order
cell_metadata <- cell_metadata[match(colnames(seurat_obj), cell_metadata$cell), ]

# Assign metadata columns
seurat_obj$geo_accession <- cell_metadata$geo_accession
seurat_obj$group         <- cell_metadata$group
seurat_obj$status        <- cell_metadata$status
seurat_obj$sample_type   <- cell_metadata$cell_type
seurat_obj$gender        <- cell_metadata$gender
seurat_obj@meta.data$patient = seurat_obj@meta.data$geo_accession

# -------------------------------
# STEP 5: Save Seurat Object
# -------------------------------
output_path <- file.path(base_path, "processed", dataset)
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

saveRDS(seurat_obj, file = file.path(output_path, "seurat_obj_with_metadata.rds"))
