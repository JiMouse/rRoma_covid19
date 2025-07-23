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

# Load libraries
suppressMessages({
  library(Seurat)
  library(Matrix)
  library(tidyverse)
  library(DropletUtils)
  library(data.table)
})

# STEP 0: Define path
path = "/Users/jean-marie/Documents/Privé/Curie/covid19"
setwd(path)
dataset <- "GSE176269"

# STEP 1 : load expression matrix
#raw_counts <- read_tsv("data/GSE176269/GSE176269_CovidStudy_rawCounts_061721.txt", col_names = TRUE, show_col_types = FALSE)
#raw_counts <- column_to_rownames(raw_counts, var = names(raw_counts)[1])
#saveRDS(object = raw_counts, file = "data/GSE176269/GSE176269_CovidStudy_rawCounts_061721.rds")
raw_counts <- readRDS("data/GSE176269/GSE176269_CovidStudy_rawCounts_061721.rds")
dim(raw_counts)

#### CONVERT AS MATRIX IN CHUNK CAUSE OF MEMORY ISSUES.
# Function to convert block to sparse
make_sparse_block <- function(df_block, start_row = 1) {
  # Find non-zero entries
  nonzero_idx <- which(df_block != 0, arr.ind = TRUE)
  
  # Construct sparse matrix
  sparse_mat <- sparseMatrix(
    i = nonzero_idx[, 1] + (start_row - 1),
    j = nonzero_idx[, 2],
    x = df_block[nonzero_idx],
    dims = c(nrow(raw_counts), ncol(raw_counts)),
    dimnames = list(rownames(raw_counts), colnames(raw_counts))
    #dimnames = list(rownames(df_block), colnames(raw_counts))
    
  )
  
  return(sparse_mat)
}


chunk_size <- 10
n_chunks <- ceiling(nrow(raw_counts) / chunk_size)


# Prepare to collect i, j, x for sparseMatrix()
i_list <- list()
j_list <- list()
x_list <- list()

for (chunk_id in seq_len(n_chunks)) {
  row_start <- (chunk_id - 1) * chunk_size + 1
  row_end <- min(chunk_id * chunk_size, nrow(raw_counts))
  
  cat("Processing rows:", row_start, "to", row_end, "\n")
  
  chunk <- raw_counts[row_start:row_end, , drop = FALSE]
  nonzero_idx <- which(chunk != 0, arr.ind = TRUE)
  
  if (nrow(nonzero_idx) > 0) {
    # Adjust row indices to global matrix
    i_list[[chunk_id]] <- row_start - 1 + nonzero_idx[, 1]
    j_list[[chunk_id]] <- nonzero_idx[, 2]
    x_list[[chunk_id]] <- chunk[nonzero_idx]
  }
}

# Combine all parts
i_all <- unlist(i_list)
j_all <- unlist(j_list)
x_all <- unlist(x_list)

# Construct sparse matrix once
sparse_all <- sparseMatrix(
  i = i_all,
  j = j_all,
  x = x_all,
  dims = c(nrow(raw_counts), ncol(raw_counts)),
  dimnames = list(rownames(raw_counts), colnames(raw_counts))
)

# Save as HDF5 using DropletUtils
DropletUtils::write10xCounts("data/GSE176269/tenx_out", sparse_all)

#Create seurat object
seurat_obj <- Read10X("data/GSE176269/tenx_out") %>% CreateSeuratObject()

# check
rownames(seurat_obj@assays$RNA) %>% tail(100)
dim(seurat_obj@assays$RNA)

# Load metadata
meta_data <- read_tsv("data/GSE176269/GSE176269_CovidStudy_phenotype_061721.txt", show_col_types = FALSE)
meta_data <- column_to_rownames(meta_data, var = names(meta_data)[1])

# Construct SeuratObject
seurat_obj <- AddMetaData(object = seurat_obj, metadata = meta_data)
seurat_obj@meta.data %>% str

#Update metadata feature to match following scripts
seurat_obj$celltype = seurat_obj$cellType
#seurat_obj$patient = seurat_obj$sampID
seurat_obj$patient = str_c("id_", seurat_obj@meta.data$donorID)
seurat_obj$series <- case_when(
  grepl("[a]$", seurat_obj$sampID) ~ 2,
  grepl("[b]$", seurat_obj$sampID) ~ 3,
  grepl("[c]$", seurat_obj$sampID) ~ 4,
  TRUE ~ 1
)

seurat_obj<-UpdateSeuratObject(seurat_obj)

saveRDS(object = seurat_obj, file = str_c("processed", dataset, "seurat_obj_with_metadata.rds", sep="/"))
