# COVID-19 scRNA-seq pseudobulk analysis in R
# Author: Jean-Marie Ravel Institut Curie
# Focus: Construction of a Seurat object for following pipeline

#clean environment
rm(list=ls(all.names=TRUE)); gc(); if (dev.cur() > 1) dev.off(); cat("\014")

#Install libraries
# R 4.4.3
# Install if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) BiocManager::install("BiocManager")
if (!requireNamespace("HDF5Array", quietly = TRUE)) BiocManager::install("HDF5Array")
if (!requireNamespace("DelayedArray", quietly = TRUE)) BiocManager::install("DelayedArray")
if (!requireNamespace("DropletUtils", quietly = TRUE)) BiocManager::install("DropletUtils")
if (!requireNamespace("GEOquery", quietly = TRUE))  BiocManager::install("GEOquery")

# Load libraries
suppressMessages({
  library(Seurat)
  library(Matrix)
  library(tidyverse)
  library(DropletUtils)
  library(data.table)
  library(GEOquery)
})


# STEP 0: Define path
path = "/Users/jean-marie/Documents/Privé/Curie/covid19"
setwd(path)
dataset <- "GSE216020"

# ------------------------------
# 2. Download and Extract Data
# ------------------------------

# Download the data
setwd(str_c(path, "data", sep="/"))
#getGEOSuppFiles("GSE216020", makeDirectory = F)

# Unzip the files
#untar("GSE216020/GSE216020_RAW.tar", exdir = "GSE216020/raw")

# ------------------------------
# 3. Load Series Matrix (Metadata)
# ------------------------------
gset <- getGEO("GSE216020", GSEMatrix = TRUE)
pheno <- pData(gset[[1]])
str(pheno)
#write_tsv(x = pheno, file = "GSE216020/pheno.tsv")

# Clean and simplify metadata
metadata <- pheno %>%
  select(title, geo_accession, description, "genotype:ch1", "time:ch1")

head(metadata)

# ------------------------------
# 4. Read 10X Raw Data
# ------------------------------
library(Matrix)
library(Seurat)
library(tidyverse)

raw_dir <- "GSE216020/raw"
files <- list.files(raw_dir, full.names = TRUE)

# Extract unique GSM IDs
gsm_ids <- unique(str_extract(basename(files), "GSM\\d+"))

seurat_list <- list()

for (gsm in gsm_ids) {
  message("Processing sample: ", gsm)
  
  # Identify the corresponding files
  barcode_file <- list.files(raw_dir, pattern = paste0(gsm, ".*barcodes.tsv.gz$"), full.names = TRUE)
  feature_file <- list.files(raw_dir, pattern = paste0(gsm, ".*features.tsv.gz$"), full.names = TRUE)
  matrix_file  <- list.files(raw_dir, pattern = paste0(gsm, ".*matrix.mtx.gz$"), full.names = TRUE)
  
  # Read files
  mat <- readMM(matrix_file)
  features <- read.delim(gzfile(feature_file), header = FALSE)
  barcodes <- read.delim(gzfile(barcode_file), header = FALSE)
  
  # Assign row/column names to matrix
  rownames(mat) <- make.unique(features$V2)  # V2 = gene name
  colnames(mat) <- barcodes$V1
  
  # Create Seurat object
  seu <- CreateSeuratObject(counts = mat, project = "GSE216020", min.cells = 3, min.features = 200)
  
  # Add sample ID
  seu$sample_id <- gsm
  
  # Add metadata if available
  if (gsm %in% metadata$geo_accession) {
    sample_meta <- metadata %>% filter(geo_accession == gsm)
    seu$patient <- sample_meta$description
    seu$status <- sample_meta[["genotype:ch1"]]
    seu$series <- sample_meta[["time:ch1"]]
    seu$title <- sample_meta$title
  }
  
  # Store
  seurat_list[[gsm]] <- seu
}
seurat_obj <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = names(seurat_list))
seurat_obj<-UpdateSeuratObject(seurat_obj)

# Join Layers
seurat_obj <- JoinLayers(seurat_obj)
seurat_obj

#Save
setwd(path)
output_path <- str_c("processed/", dataset)
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

saveRDS(object = seurat_obj, file = str_c("processed", dataset, "seurat_obj_with_metadata.rds", sep="/"))
