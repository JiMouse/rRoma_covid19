### Set up
rm(list=ls()) ; cat("\014") ; gc() ; graphics.off() #dev.off()

library(tidyverse)
library(Seurat)
library(edgeR)

# based on edger User Guide
# 4.10 Single cell RNA-seq differential expression with pseudobulking
# https://www.nature.com/articles/s41597-022-01236-2
# https://hbctraining.github.io/scRNA-seq/lessons/04_SC_quality_control.html

source(file = "C:\\Users\\Administrateur\\Documents\\PhD\\covid19\\_NaviCell\\Scripts\\scRNA\\AveragePseudoBulk.R")

#### select dataset ######################################
# dataset = "12436517"
dataset = "GSE176269"

#### import parameters ######################################
if(dataset == "12436517"){
  filepath = "12436517_processed_final.rds"
  cluster_id = "predicted.ann_level_3" #celltype
  sample_id = "patient"
  group = "severity"
}

if(dataset == "GSE176269"){
  filepath = "CovidStudy_rawCounts_061721.rds"
  cluster_id = "cellType"
  sample_id = "sampID"
  group = "status"
}

##
#Input is a rds file containing a Seurat object.
##

## Import dataset ############################################
path = paste0("C:/Users/Administrateur/Documents/PhD/covid19/_Navicell/Datasets/", dataset, "/")
setwd(dir = paste0(path, "output/"))
seurat = readRDS(paste0(path, filepath))
celltype_builder = FALSE

# dataset specific paramaters #################################
if (dataset == "12436517") {
  seurat = subset(x=seurat, subset = sampleType == "NS")
}

#### Export groups ############################################
groups_df <- seurat@meta.data %>% rownames_to_column("barcodes")
patient_groups_df <- groups_df[c(group, sample_id)] %>%
  na.omit() %>%
  distinct() %>%
  # mutate(severity=str_replace_all(severity, "critical", "xcritical")) %>%
  arrange(get(group)) #%>%
  # mutate(severity=str_replace_all(severity, "xcritical", "critical"))

write_tsv(x = patient_groups_df, file = "patient_group.tsv")

#### QC control scRNA ###########################
if (!("percent.mt" %in% colnames(seurat@meta.data))) {
  seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
}

# Visualize QC metrics as a violin plot
# VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0)

# Cell-level filtering
filtered_seurat <- subset(seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 25)

# Gene-level filtering
# Extract counts
counts <- GetAssayData(object = filtered_seurat, slot = "counts")

# Output a logical vector for every gene on whether the more than zero counts per cell
nonzero <- counts > 0

# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10

# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]

# Reassign to filtered Seurat object
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)

##### Grouping and filtering
seurat <- filtered_seurat
# Create new groups
if (celltype_builder) {
  seurat@meta.data[["celltype2"]] = seurat@meta.data[[cluster_id]] %>%
    str_replace_all(celltype_builder)
} else {
  seurat@meta.data[["celltype2"]] = seurat@meta.data[[cluster_id]]
}

# filter cells with poor celltype annotation
if(cluster_id == "predicted.ann_level_3") {
  seurat <- subset(seurat, subset = predicted.ann_level_3.score > 0.95)
}

#update cluster_id
cluster_id = "celltype2"

#remove cells with NA in cluster
seurat <- subset(seurat, subset = celltype2 != "NA") #subset NA

# set cluster_id as default ident
seurat <- SetIdent(seurat, value = cluster_id)

# update metadata to replace all "- and _" by space
seurat@meta.data[[cluster_id]] <- seurat@meta.data[[cluster_id]] %>% str_replace_all(., "-|_", " ") %>% factor

#check cell composition
# seurat@meta.data[[cluster_id]] %>% table(useNA = "ifany")

# Create pseudo-bulk samples
y<-Seurat2PB(seurat, sample=sample_id, cluster=cluster_id)

# extract raw count matrix
lcpm=y$counts

# Not every cluster is present in all samples; create a vector that represents how to split samples
pb = t(lcpm)
splitf <- sapply(stringr::str_split(rownames(pb), 
                                    pattern = "_",  
                                    n = 2), #Maximum number of pieces to return
                 `[`, 2)

pb <- split.data.frame(pb, factor(splitf))%>%
  lapply(function(u)
    magrittr::set_colnames(t(u),
                           stringr::str_remove(rownames(u), "cluster")))

names(pb) <- names(pb) %>% stringr::str_remove("cluster")
# names(pb)

# export by cluster
for(celltype in names(pb)) {
  ## extract counts
  df = pb[[celltype]]
  colnames(df) <- colnames(df) %>% str_split_i(., "_", 1)
  
  ## create DGE object
  dge <- DGEList(df)
  # dge$samples
  # summary(dge$samples$lib.size)
  # barplot(dge$samples$lib.size*1e-6, names=1:(dge$samples$lib.size %>% length))
  
  ## filter out sample with low lib.size
  m=dge$samples$lib.size %>% median(.)/10 #threshold = median/10
  keep.samples<-dge$samples$lib.size>m
  dge<-dge[,keep.samples]
  
  #filter out lowly expressed genes.
  keep.genes<-filterByExpr(dge, group=dge$samples$cluster)#, min.count=10,min.total.count=20)
  dge<-dge[keep.genes,,keep=FALSE]
  
  if(sum(dge$counts) == 0 | dim(dge)[1] < 10) next #skip if no counts
  
  dge <- calcNormFactors(dge, method = "TMM")
  logCPM <- edgeR::cpm(dge, log=TRUE)
  
  # export TSV and boxplot
  png(filename = paste0(dataset, "_", celltype, "_boxplot.png"), width = 1200, height = 800, units = "px", pointsize = 14)
  boxplot(logCPM)
  dev.off()
  
  logCPM <- logCPM %>% as.matrix() %>% as_tibble(rownames="gene")
  write_tsv(x = logCPM, file = paste0(dataset, "_", celltype, "_pseudobulk", ".tsv"))
  
  #export average for Navicell visualisation
  exportAveragePseudoBulk(df, patient_groups_df, celltype, sample_id, group, subdir, dataset)
  
  print(paste0(celltype, " export is done !"))
}
print('All of the subtypes have been successfully exported!')
