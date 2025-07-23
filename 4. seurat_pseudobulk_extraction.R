# 4th script : export pseudobulk count from seurat object
#last update : 21/07/2025
#https://www.bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
#https://hbctraining.github.io/scRNA-seq/lessons/pseudobulk_DESeq2_scrnaseq.html

#clean environment
rm(list=ls(all.names=TRUE)); gc(); if (dev.cur() > 1) dev.off(); cat("\014")

#Install Matrix.utils library
#install.packages('devtools')
#devtools::install_url("https://cran.r-project.org/src/contrib/Archive/Matrix.utils/Matrix.utils_0.9.7.tar.gz")

# Load libraries
suppressMessages({
  library(Seurat)
  library(Matrix)
  library(tidyverse)
  library(edgeR)
  library(SingleCellExperiment)
  library(Matrix.utils)
})

# STEP 0: Define path
path = "/Users/jean-marie/Documents/Privé/Curie/covid19" #"C:/Users/Administrateur/Documents/covid19/" #"/Users/jean-marie/Documents/Privé/Curie/covid19"
setwd(path)

dataset_list <- c("GSE149689")# not working, too heavy, "GSE161918")#, "GSE171524")
dataset <- "GSE149689" #"12436517" "GSE176269" #"GSE171524" #"GSE216020" #"GSE161918"

for(dataset in dataset_list) {
  dir_path <- str_c("processed", dataset, sep="/")
  
  # OPTIONS DEPENDING OF DATASET
  series_subset = 1
  if(series_subset!=1) dataset=str_c(dataset,series_subset,sep="_")
  
  # DEFINE OUTPUT PATH
  output_path <- str_c("pseudobulk/", dataset)
  dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
  
  #STEP 1 load RDS object
  seurat_obj <- readRDS(file = str_c(dir_path, "seurat_obj_SingleR_annotated.rds", sep="/"))
  seurat_obj@meta.data$singleR.labels <- seurat_obj@meta.data$singler_clust_labels
  seurat_obj$newcell <- ifelse(seurat_obj$singleR.labels =='Macrophage' | seurat_obj$singleR.labels =='Monocyte', 'Mono_Mac', seurat_obj$singleR.labels)
  
  message("SingleR annotation is used by default!")
  seurat_obj$celltype <- seurat_obj$newcell
  seurat_obj@meta.data %>% str()
  
  # #celltype or clust_labels_col
  # if (is.null(seurat_obj@meta.data$celltype)) {
  #   message("No celltype. SingleR annotation is used!")
  #   seurat_obj@meta.data$celltype <- seurat_obj@meta.data$singler_clust_labels
  # }
  
  # subset to keep only NS samples of 1 time
  if(dataset=="12436517"){
    seurat_obj@meta.data$status = seurat_obj@meta.data$severity
    seurat_obj = subset(x=seurat_obj, subset = location == "NS")
    seurat_obj = subset(x=seurat_obj, subset = series == series_subset)
  }
  if(dataset=="GSE176269_2" && series_subset == 2){
    #seurat_obj = subset(x=seurat_obj, subset = series == series_subset)
    seurat_obj$celltype <- case_when(
      grepl("mac", seurat_obj$celltype) ~ "macrophage",
      grepl("G5", seurat_obj$celltype) ~ "neutrophile",
      grepl("neut", seurat_obj$celltype) ~ "neutrophile",
      grepl("epi|ciliated|basal|hillock|goblet", seurat_obj$celltype) ~ "epithelial",
      grepl("B", seurat_obj$celltype) ~ "B cell",
      grepl("T", seurat_obj$celltype) ~ "T cell",
      .default = as.character(seurat_obj$celltype)
    )
  }
  
  seurat_obj@meta.data$celltype %>% table() %>% sort()
  seurat_obj@meta.data$patient  %>% table()
  
  # Extract raw counts and metadata to create SingleCellExperiment object
  counts <- LayerData(seurat_obj[["RNA"]], layer = "counts") #LayerData(seurat_obj[["RNA"]], layer = "counts") #seurat_obj[["RNA"]]@counts
  metadata <- seurat_obj@meta.data
  
  # Set up metadata as desired for aggregation and DE analysis
  metadata$cluster_id <- factor(gsub("_|\\+", replacement = "-", seurat_obj@meta.data[["celltype"]])) #
  metadata$sampleInfo <- factor(seurat_obj@meta.data[["patient"]])
  
  # Create single cell experiment object
  sce <- SingleCellExperiment(assays = list(counts = counts), 
                              colData = metadata)
  
  # Identify groups for aggregation of counts
  groups <- colData(sce)[, c("cluster_id", "sampleInfo")]
  
  # Acquiring necessary metrics for aggregation across cells in a sample
  # Named vector of cluster names
  kids <- purrr::set_names(levels(sce$cluster_id))
  
  # Total number of clusters
  nk <- length(kids)
  
  # Named vector of sample names
  sids <- purrr::set_names(levels(sce$sampleInfo))
  
  # Total number of samples 
  ns <- length(sids)
  ns
  # Generate sample level metadata
  ## Determine the number of cells per sample
  table(sce$sampleInfo)
  
  ## Turn named vector into a numeric vector of number of cells per sample
  n_cells <- as.numeric(table(sce$sampleInfo))
  
  ## Determine how to reoder the samples (rows) of the metadata to match the order of sample names in sids vector
  m <- match(sids, sce$sampleInfo)
  
  ## Create the sample level metadata by combining the reordered metadata with the number of cells corresponding to each sample.
  ei <- data.frame(colData(sce)[m, ], 
                   n_cells, row.names = NULL) %>% 
    dplyr::select(-"cluster_id")
  #ei
  
  # Subset metadata to only include the cluster and sample IDs to aggregate across
  groups <- colData(sce)[, c("cluster_id", "sampleInfo")]
  
  # Aggregate across cluster-sample groups
  pb <- Matrix.utils::aggregate.Matrix(t(counts(sce)), 
                         groupings = groups, fun = "sum") 
  
  # Not every cluster is present in all samples; create a vector that represents how to split samples
  splitf <- sapply(stringr::str_split(rownames(pb), 
                                      pattern = "_",  
                                      n = 2), 
                   `[`, 1)
  # Turn into a list and split the list into components for each cluster and transform, 
  #so rows are genes and columns are samples and make rownames as the sample IDs
  pb <- split.data.frame(pb, 
                         factor(splitf)) %>%
    lapply(function(u)
      magrittr::set_colnames(t(u),
                             stringr::str_extract(rownames(u), "(.*)"))) #.* (?<=_)[:alnum:]+
  
  #Check number of cells
  cell_counts <- table(sce$cluster_id, sce$sampleInfo)
  cell_counts
  sce$sampleInfo %>% table() %>% sum()
  
  # Export to TSV
  output_file <- str_c(output_path, "/", "_cell_counts.tsv")
  cell_counts_df <- as.data.frame.matrix(cell_counts)
  cell_counts_df <- tibble::rownames_to_column(cell_counts_df, var = "cluster_id")
  write.table(cell_counts_df, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Extract sample with minimum 10 cells and export TSV for each celltype
  # Minimum required cells
  min_cells <- 5 # 5 10
  min_total_counts <- 5e4
  names(pb)
  
  cluster = "macrophage" #"M1-mac" #"unk-neut" #"M2-mac" #"nrMa"
  
  # STEP: Process each cluster
  for (cluster in names(pb)) {
    mat <- pb[[cluster]]  # genes × samples matrix
    #min_total_counts <- 5e4
    
    # Extract sample names from column names
    sample_ids <- gsub(paste0("^", cluster, "_"), "", colnames(mat))
    
    # Step 2: Look up cell counts
    counts_for_cluster <- cell_counts[cluster, sample_ids] #BUG
    keep <- as.numeric(counts_for_cluster) >= min_cells
    
    # Step 3: Filtered samples with less than 10 cells
    mat_filtered <- mat[, keep, drop = FALSE]
    colnames(mat_filtered) <- sample_ids[keep]  # remove prefix in column names
    
    # Step 4: Filter columns (samples) with total counts ≥ 50 000
    sample_totals <- colSums(mat_filtered)
    keep_samples <- sample_totals >= min_total_counts
    mat_filtered <- mat_filtered[, keep_samples, drop = FALSE]
    
    # Filter genes with low expression across samples
    keep_genes <- rowSums(mat_filtered >= 5) >= 2  # at least 5 counts in 2 samples
    mat_filtered <- mat_filtered[keep_genes, ]
    
    #checks
    #sample_totals %>% sort() %>% plot()
    summary(as.vector(mat_filtered))
    nrow(mat_filtered)
    summary(sample_totals)
    sample_totals %>% sort()
    dim(pb[[cluster]])
    
    # Step 6: Global normalization using edgeR (CPM)
    if (length(mat_filtered) > 0 && ncol(mat_filtered) > 2) { #if at least 3 samples
  
      
      # library_sizes <- colSums(mat_filtered)
      # barplot(sort(library_sizes), las = 2, main = paste0(cluster, ": library sizes"), ylab = "Total counts")
      # abline(h = 50000, col = "red", lty = 2)
      # 
      
      dge <- edgeR::DGEList(counts = mat_filtered)
      dge <- edgeR::calcNormFactors(dge)
      log_cpm_matrix <- edgeR::cpm(dge, log = TRUE) #linear normalisation
      
      # Step 6: Prepare data for ggplot
      cpm_df <- as.data.frame(log_cpm_matrix)
      cpm_df$gene <- rownames(cpm_df)
      cpm_long <- pivot_longer(cpm_df, cols = -gene, names_to = "patient_id", values_to = "log2_cpm")
      
      # Step 7: QA expression Plot
      png(str_c(output_path, "/", cluster, "_Pseudobulk_GlobalNormalization_Boxplot.png"))
      print(
        ggplot(cpm_long, aes(x = patient_id, y = log2_cpm, fill = patient_id)) +
          geom_boxplot(outlier.shape = NA) +
          theme_minimal() +
          labs(title = str_c( cluster, "- Pseudobulk Global Normalization"),
               y = "log2(CPM)",
               x = "Patient ID") +
          theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
      )
      dev.off()
      
    # Step 8: Export
      output_file <- str_c(output_path, "/", cluster, "_pseudobulk.tsv")
      write.table(log_cpm_matrix, file = output_file, sep = "\t", quote = FALSE, col.names = NA)
      
      # export raw counts
      output_file <- str_c(output_path, "/", cluster, "_raw_pseudobulk.tsv")
      write.table(mat_filtered, file = output_file, sep = "\t", quote = FALSE, col.names = NA)
      
      message("Exported: ", output_file)
    } else {
      message("Skipped (not enough cells in any sample): ", cluster)
    }
  }
  
  #Export group file if not exist
  patient_group.exist <- list.files(path = str_c("data/", dataset, sep="/"), pattern = "patient_group.tsv") %>% length() > 0
  
  #if(!patient_group.exist){
    patient_group <- metadata %>% as_tibble() %>%
      dplyr::select(patient, status) %>%
      unique()
    output_path <- str_c("data/", dataset, sep="/")
    dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
    patient_group <- write_tsv(patient_group, file = str_c(output_path, "patient_group.tsv", sep="/"))
    }
#}