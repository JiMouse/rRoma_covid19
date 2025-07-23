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
if (!requireNamespace("GEOquery", quietly = TRUE))  BiocManager::install("GEOquery")

# Load libraries
suppressMessages({
  library(Seurat)
  library(Matrix)
  library(tidyverse)
  library(data.table)
  library(GEOquery)
})


# STEP 0: Define path
path = "C:/Users/Administrateur/Documents/covid19/"#"/Users/jean-marie/Documents/Privé/Curie/covid19"
setwd(path)
dataset <- "GSE161918"

# ------------------------------
# 1. Load RDS file directly
# ------------------------------
seurat_obj <- readRDS(file = "data/GSE161918/GSE161918_AllBatches_SeuratObj.rds")
seurat_obj@meta.data %>% names()
seurat_obj@meta.data$mergedcelltype %>% table()
seurat_obj@meta.data$coursecelltype %>% table()
seurat_obj@meta.data$adjustedcelltype %>% table()

seurat_obj$patient <- seurat_obj$Subject
seurat_obj$status <- seurat_obj$severity
levels(seurat_obj$status) <- c(levels(seurat_obj$status), "HC")
seurat_obj$status[is.na(seurat_obj$status)] <- rep("HC", length(seurat_obj$status[is.na(seurat_obj$status)]))
seurat_obj$status %>% table()
seurat_obj@meta.data$celltype <- seurat_obj@meta.data$mergedcelltype                                                   

#Upadte to seurat v5
seurat_obj <- UpdateSeuratObject(seurat_obj)

# Join Layers
seurat_obj <- JoinLayers(seurat_obj)

#Save
setwd(path)
output_path <- str_c("processed/", dataset)
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

saveRDS(object = seurat_obj, file = str_c("processed", dataset, "seurat_obj_with_metadata.rds", sep="/"))


################# OR #############################

# ------------------------------
# 2. Download and Extract Data
# ------------------------------

# Download the data #NOT WORKING
setwd(str_c(path, "data", sep="/"))
#getGEOSuppFiles(dataset, makeDirectory = F, fetch_files=T, filter_regex="_RAW.tar")

# Manually download the file 
# url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE161nnn/GSE161918/suppl/GSE161918_RAW.tar"
# download.file(url, destfile = "GSE161918_RAW.tar", mode = "wb")

# Unzip the files
raw_dir <- str_c(dataset, "raw", sep= "/")
#untar(str_c(dataset, str_c(dataset, "_RAW.tar"), sep="/"), exdir = raw_dir)

# ------------------------------
# 3. Load Series Matrix (Metadata)
# ------------------------------
# gset <- getGEO(dataset, GSEMatrix = TRUE)
# pheno <- pData(gset[[1]])
# str(pheno)
# write_tsv(x = pheno, file = str_c(dataset, "pheno.tsv", sep="/"))

pheno <- read_tsv(str_c(dataset, "pheno.tsv", sep="/"), show_col_types = F)

# Clean and simplify metadata
metadata <- pheno %>%
  dplyr::select(title, geo_accession)

head(metadata)

# ------------------------------
# 4. Read 10X Raw Data
# ------------------------------
files <- list.files(raw_dir, full.names = TRUE, pattern = "*.h5")

# Load all into a list
seurat_list <- lapply(files, function(f) {
  cat(f)
  data <- Read10X_h5(f)
  if(type(data)=="list"){
    CreateSeuratObject(counts = data$Gene, project = f)
  }
  else CreateSeuratObject(counts = data, project = f)
})

# Optionally add sample names
for (i in seq_along(seurat_list)) {
  seurat_list[[i]]$sample <- gsub("_filtered_feature_bc_matrix.h5", "", files[i])
}

# Merge all Seurat objects
seurat_obj <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = basename(tools::file_path_sans_ext(files)), project = "GSE161918")

# ------------------------------
# 4. Update, Join and Export Seurat object
# ------------------------------


seurat_obj<-UpdateSeuratObject(seurat_obj)

# Join Layers
seurat_obj <- JoinLayers(seurat_obj)
seurat_obj

#Save
setwd(path)
output_path <- str_c("processed/", dataset)
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

saveRDS(object = seurat_obj, file = str_c("processed", dataset, "seurat_obj_with_metadata.rds", sep="/"))
