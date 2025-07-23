# COVID-19 scRNA-seq pseudobulk analysis in R
# Author: Jean-Marie Ravel Institut Curie
# Focus: Annotation + pseudobulk for macrophages + normalization check + export

#clean environment
rm(list=ls(all.names=TRUE)); gc(); if (dev.cur() > 1) dev.off(); cat("\014")

#Install libraries
# R 4.4.3

# Load libraries
suppressMessages({
  library(Seurat)
  library(Matrix)
  library(tidyverse)
  library(sctransform)
})

# STEP 0: Define path
path = "/Users/jean-marie/Documents/Privé/Curie/covid19" #"C:/Users/Administrateur/Documents/covid19/"#"/Users/jean-marie/Documents/Privé/Curie/covid19"
setwd(path)
dataset <- "GSE149689" #"12436517" "GSE176269" #"GSE171524" #"GSE216020" #"GSE161918" "GSE149689"

# STEP 1 load Seurat object
seurat_obj = readRDS(str_c("processed",dataset,"seurat_obj_with_metadata.rds", sep="/"))
seurat_obj

# STEP 2: Basic filtering and normalization
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# OPTIONAL WITH SCTtransform : https://satijalab.org/seurat/articles/sctransform_vignette.html
# but useful for 
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 20)

# STEP 3 - Normalization
# Option 1 : Apply SCTransform
# Option 1a
#seurat_obj <- SCTransform(seurat_obj, verbose = FALSE)

#option 1b if memory failed
# seurat_list <- SplitObject(seurat_obj, split.by = "patient")
# seurat_list <- lapply(seurat_list, function(obj) {
#   Seurat::SCTransform(obj, verbose = FALSE, vars.to.regress = "percent.mt")
# })

# Save RDS
#saveRDS(object = seurat_list, file = str_c("processed",dataset,"seurat_obj_tmp.rds", sep="/"))
#dataset <- "GSE216020" ; seurat_list = readRDS(str_c("processed",dataset,"seurat_obj_tmp.rds", sep="/"))

# Merge Option 1
#seurat_obj <- merge(x = seurat_list[[1]], y = seurat_list[-1])
#seurat_obj[["RNA"]] <- JoinLayers(seurat_obj[["RNA"]])
#seurat_obj <- JoinLayers(seurat_obj)
#DefaultAssay(seurat_obj) <- "SCT"

# Option 2 : Normalize (in stead of SCT)
seurat_obj <- NormalizeData(object = seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)

# STEP 3 : PCA and UMAP
seurat_obj <- FindVariableFeatures(seurat_obj)#, assay = "SCT")
seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30) #20

# STEP 4: Clustering
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30) #20
seurat_obj <- FindClusters(seurat_obj, resolution = 0.8)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.3)
DimPlot(seurat_obj, label = TRUE)

# Save RDS
saveRDS(object = seurat_obj, file = str_c("processed", dataset, "seurat_obj_with_clusters_features.rds", sep = "/"))
