# COVID-19 scRNA-seq pseudobulk analysis in R
# Author: Jean-Marie Ravel Institut Curie
# Focus: SingleR Annotation for Seurat object with PCA, UMAP etc.
# https://pmc.ncbi.nlm.nih.gov/articles/PMC9709423/?utm_source=chatgpt.com#s2

#clean environment
rm(list=ls(all.names=TRUE)); gc(); if (dev.cur() > 1) dev.off(); cat("\014")

# General information
#https://moodle.france-bioinformatique.fr/pluginfile.php/907/mod_resource/content/4/EBAII_Cell-type_annotation.html#automatic-annotation-using-a-reference-scrnaseq
#SingleR works better with the normalized data

#Install libraries
# R 4.4.3
if (!requireNamespace("SingleR", quietly = TRUE)) BiocManager::install("SingleR")
if (!requireNamespace("edgeR", quietly = TRUE)) BiocManager::install("edgeR")
if (!requireNamespace("celldex", quietly = TRUE)) BiocManager::install("celldex")

# Load libraries
suppressMessages({
  library(Seurat)
  library(Matrix)
  library(SingleR)
  library(celldex)
  library(edgeR)
  library(tidyverse)
})

# STEP 0: Define path
path = "/Users/jean-marie/Documents/Privé/Curie/covid19" #"C:/Users/Administrateur/Documents/covid19/" #"/Users/jean-marie/Documents/Privé/Curie/covid19"
setwd(path)
dataset_list <- c("GSE149689") #"12436517", "GSE176269", "GSE216020")#, "GSE161918")#, "GSE171524")
#dataset <- "GSE176269" #"12436517" "GSE176269" #"GSE171524" #"GSE216020" #"GSE161918"

for (dataset in dataset_list) {

dir_path <- str_c("processed", dataset, sep="/")

# STEP 1: load and update to current version of Seurat
seurat_obj = readRDS(str_c(dir_path,"seurat_obj_with_clusters_features.rds", sep="/"))
#DefaultAssay(seurat_obj) <- "SCT"

cluster_label = str_c(DefaultAssay(seurat_obj),"_snn_res.0.8") #0.8

# STEP 2: LOAD REFERENCE
# ref_ImmGen <- celldex::ImmGenData(ensembl = FALSE)
hpca.ref <- celldex::HumanPrimaryCellAtlasData(ensembl = FALSE)
#hpca.ref <- readRDS("data/ref_HumanPrimaryCellAtlasData.rds")
#ref_blueprint <- celldex::BlueprintEncodeData()

# STEP 3: Extraction of the normanised data
clust_labels_col = "singler_clust_labels"
#equivalent to seurat_obj[["SCT"]]$scale.data
norm.data = Seurat::GetAssayData(
  object = seurat_obj,
  assay = "RNA", #RNA if using normalize, else SCT
  layer = "data" #normalised count matrix #data if using normalize
)

# STEP 4: prediction using clustering information 
# OPTION 1 : based on cluster
# singleR.results  =
#   SingleR::SingleR(
#     test = norm.data,
#     clusters = seurat_obj[[cluster_label]][, 1], #seurat_obj$RNA_snn_res.0.3
#     ref = hpca.ref,
#     labels = hpca.ref$label.main, #label.main, #label.fine,
#     #assay.type.test = "logcounts",
#     #assay.type.ref = "logcounts",
#     #BPPARAM = BiocParallel::SerialParam()
#     )
# # if cluster used -> return one cell type name for each cluster
# 
# table(singleR.results$labels) %>% sort(., decreasing = TRUE)
# summary(is.na(singleR.results$pruned.labels))
# 
# # IF cluster is used
# # Save the name of future annotation
# 
# # Create a column with this name in the metadata and fill it with the cluster levels of each cell
# seurat_obj@meta.data[[clust_labels_col]] = seurat_obj@meta.data[[cluster_label]]
# # Fill associate each cluster with its annotation 
# levels(seurat_obj@meta.data[[clust_labels_col]]) = singleR.results$labels


# OPTION 2 : not cluster-based
singleR.results  =
  SingleR::SingleR(
    test = norm.data,
    ref = hpca.ref,
    labels = hpca.ref$label.main, 
  )
#seurat_obj@meta.data[[clust_labels_col]] <- seurat_obj@meta.data[[cluster_label]] #create column
seurat_obj@meta.data[[clust_labels_col]] <- singleR.results$labels #fill it

#check cluster
table(seurat_obj@meta.data[[clust_labels_col]]) %>% sort()


# STEP: Check corelation between celltype and single R for mrMa
# annotation_df <- tibble(
#   celltype = seurat_obj$celltype,
#   singleR = seurat_obj$SingleR_label
# ) %>% 
#   filter(celltype == "M1-mac")
# 
# celltype_table <- table(annotation_df$celltype, annotation_df$singleR) %>% 
#   as.data.frame(.) %>%
#   arrange(desc(Freq))
# celltype_table


# set color palette to use
seeable_palette = setNames(
  c(RColorBrewer::brewer.pal(name = "Dark2", n = 8),
    c(1:(length(unique(singleR.results$labels)) - 8))),
  nm = names(sort(table(singleR.results$labels), decreasing = TRUE)))

singleR_cluster_plot = Seurat::DimPlot(
  object = seurat_obj, 
  reduction = "umap", 
  group.by = clust_labels_col,
  pt.size = 2,
  label = FALSE, 
  cols = seeable_palette
) #+ ggplot2::theme(legend.position = "none")

ann_cell_plot = Seurat::DimPlot(
  object = seurat_obj, 
  reduction = "umap", 
  group.by = cluster_label,
  pt.size = 2,
  label = TRUE,
  repel = TRUE,
) + ggplot2::theme(legend.position = "none")

#print(singleR_cluster_plot + ann_cell_plot)

# Save RDS
saveRDS(object = seurat_obj, file = str_c(dir_path, "seurat_obj_SingleR_annotated.rds", sep="/"))

}

# STEP 4: Normalization quality check
# pdf("processed/seurat_obj_annotated_Normalization_QC_Plots.pdf")
# VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# ElbowPlot(seurat_obj)
# DimPlot(seurat_obj, reduction = "pca", group.by = "SingleR_label", label = TRUE)
# dev.off()