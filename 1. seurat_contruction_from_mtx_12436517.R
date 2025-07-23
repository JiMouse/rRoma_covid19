# COVID-19 scRNA-seq pseudobulk analysis in R
# Author: Jean-Marie Ravel, Institut Curie
# Focus: Import 12436517 dataset as Seurat raw object
# raw data and covid_nbt_main.rds were downloaded from 
#   https://figshare.com/articles/dataset/COVID-19_severity_correlates_with_airway_epithelium-immune_cell_interactions_identified_by_single-cell_analysis/12436517
# original paper : https://www.nature.com/articles/s41587-020-0602-4

#clean environment
rm(list=ls(all.names=TRUE)); gc(); if (dev.cur() > 1) dev.off(); cat("\014")

#Install libraries
# R 4.4.3
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("Seurat", quietly = TRUE)) install.packages("Seurat")
if (!requireNamespace("Seurat", quietly = TRUE)) install.packages("tidyverse")

# Load libraries
suppressMessages({
  library(Seurat)
  library(Matrix)
  library(tidyverse)
})

# STEP 1A: Load raw matrix and create Seurat object
path = "/Users/jean-marie/Documents/Privé/Curie/covid19"
setwd(path)
data_dir <- str_c(path, "data/raw_matrices", sep="/")

# Load multiple patient matrices and merge
sample_folders <- list.dirs(data_dir, full.names = TRUE, recursive = FALSE)
seurat_list <- list()
for (folder in sample_folders) {
  sample_name <- basename(folder)
  sample_data <- Read10X(data.dir = folder)
  seurat_obj_tmp <- CreateSeuratObject(counts = sample_data, project = sample_name, min.cells = 3, min.features = 200)
  seurat_obj_tmp$patient_id <- sample_name
  seurat_list[[sample_name]] <- seurat_obj_tmp
}
seurat_obj_raw <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = names(seurat_list), project = "12436517")

# Join Layers
seurat_obj_raw <- JoinLayers(seurat_obj_raw)

# Save RDS containing raw counts
saveRDS(object = seurat_obj_raw, file = "processed/seurat_obj_raw_matrices.rds")

#########################################################################
# Part 2 : add phenotypes
# Load raw counts
seurat_obj_raw <- readRDS("processed/seurat_obj_raw_matrices.rds")

# load annotated seurat (covid_nbt_main)
seurat_obj <- readRDS("data/covid_nbt_main.rds")
seurat_obj <- UpdateSeuratObject(seurat_obj)

# Step 1 : rename barcode to match
# Rename barcodes to fit
colnames(seurat_obj_raw) <- gsub("(-\\d+)+$", "", colnames(seurat_obj_raw))  # strip head "1-"

#Remove duplicate and then rename barcodes for raw object
seurat_obj_names<- str_c(seurat_obj$sample, gsub("(_\\d+)+$", "", colnames(seurat_obj)), sep = "_")  # prepend sample ID
seurat_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[!duplicated(seurat_obj_names)])
colnames(seurat_obj) <- str_c(seurat_obj$sample, gsub("(_\\d+)+$", "", colnames(seurat_obj)), sep = "_")  # prepend sample ID

# Step 2: Find matching barcodes
common_barcodes <- intersect(colnames(seurat_obj_raw), colnames(seurat_obj))
common_genes <- intersect(rownames(seurat_obj_raw[["RNA"]]), rownames(seurat_obj[["RNA"]]))

# Step 3: Subset both objects to matching cells
seurat_obj_raw <- subset(seurat_obj_raw, cells = common_barcodes)
seurat_obj <- subset(seurat_obj, cells = common_barcodes)

seurat_obj_raw[["RNA"]] <- subset(seurat_obj_raw[["RNA"]], features = common_genes)
seurat_obj[["RNA"]] <- subset(seurat_obj[["RNA"]], features = common_genes)

# STEP 4 Extract and replace raw counts
counts <- LayerData(seurat_obj_raw[["RNA"]], layer = "counts")

seurat_obj[["RNA"]] <- SetAssayData(
  object = seurat_obj[["RNA"]],
  layer = "counts",
  new.data = counts
)

# Save RDS
saveRDS(object = seurat_obj, file="processed/seurat_obj_with_metadata.rds")
