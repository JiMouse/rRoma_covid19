# rROMA analysis script
# Prerequisites: rROMA package, MSigDB hallmark gene sets, normalized expression matrix (logCPM), and patient metadata

#clean environment
rm(list=ls(all.names=TRUE)); gc(); if (dev.cur() > 1) dev.off(); cat("\014")

#Install libraries
# R 4.4.3
# Load libraries
suppressMessages({
  #if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  if (!require("biomaRt", quietly = TRUE)) BiocManager::install("biomaRt")
  if (!require("scater", quietly = TRUE)) BiocManager::install("scater")
  if(!require("tidyverse")) install.packages("tidyverse")
  if(!requireNamespace("devtools")){install.packages("devtools")}
  if(!require("rRoma")){devtools::install_github("sysbio-curie/rRoma")}
  if(!require("RColorBrewer")) install.packages("RColorBrewer")
})

# Load required packages
suppressMessages({
  library(biomaRt)
  library(scater)
  library(rRoma)
  # library(msigdbr)
  library(tidyverse)
  library(RColorBrewer)
})

# STEP 0: Define path
path = "/Users/jean-marie/Documents/Privé/Curie/covid19" #"C:/Users/Administrateur/Documents/covid19/" #"/Users/jean-marie/Documents/Privé/Curie/covid19"
setwd(path)
dataset <- "GSE149689" #"12436517" "GSE176269" #"GSE171524" #"GSE216020" #"GSE161918" "aggregated"

input_path = str_c("pseudobulk", dataset, sep = "/")

#Define output path
output_path = str_c("rRoma", dataset, "pbmc_AllHall_full", sep = "/")
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

#List of cells to analyse
celltype_list <- list.files(input_path, pattern = "_pseudobulk.tsv") %>% gsub("_pseudobulk.tsv","",.)

filter_valid_modules <- function(ModuleList, ExpressionMatrix, min_genes = 10) {
  Filter(function(m) {
    valid_genes <- intersect(m$Genes, rownames(ExpressionMatrix))
    #remove if lenth < min_genes
    if (length(valid_genes) < min_genes) return(FALSE)
    mat <- ExpressionMatrix[valid_genes, ]
    
    # if all genes of the module have no variances
    if(all(apply(mat, 1, sd) == 0)) return(FALSE)
    
    # Try PCA safely
    tryCatch({
      prcomp(t(mat))
      TRUE
    }, error = function(e) FALSE)
  }, ModuleList)
}

# Step 4: Prepare hallmark modules
AllHall <- rRoma::SelectFromMSIGdb("HALLMARK")
AllHall <- lapply(AllHall, function(x){
  x$Name <- sub("HALLMARK_", "", x$Name)
  print(x$Name)
  x
})
target_names <- c("ANGIOGENESIS","EPITHELIAL_MESENCHYMAL_TRANSITION",
                  "EPITHELIAL_MESENCHYMAL_TRANSITION","TGF_BETA_SIGNALING",
                  "TNFA_SIGNALING_VIA_NFKB","INTERFERON_ALPHA_RESPONSE", 
                  "INTERFERON_GAMMA_RESPONSE")
filtered_Hall <- AllHall[sapply(AllHall, function(x) x$Name %in% target_names)]
#AllHall <- filtered_Hall

# Step 3: Load patient metadata (must include group labels)
patient_group_file <- str_c("data", dataset, "patient_group.tsv", sep = "/")
patient_group <- read_table(patient_group_file, show_col_types = F)
patient_group <- column_to_rownames(patient_group, "patient")

#test
celltype = celltype_list[1]

# function
rRoma_analysis <- function(dataset, celltype){

  #Define output path
  output_path = str_c(output_path, celltype, sep = "/")
  
  # Step 1: Load expression matrix (from previous pseudobulk step)
  expr_file <- str_c(input_path, str_c(celltype, "_pseudobulk.tsv"), sep = "/")
  expr_matrix <- read.table(expr_file, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
  expr_matrix <- as.matrix(expr_matrix)
  
  # Select common patients
  common_patients <- intersect(colnames(expr_matrix), rownames(patient_group))
  expr_matrix <- expr_matrix[, common_patients]
  patient_group <- patient_group[common_patients, , drop = FALSE]
  group <- patient_group[["status"]] %>% sort()
  names(group) <- rownames(patient_group)
  
  # Filter only valid modules
  #safe_AllHall <- filter_valid_modules(AllHall, expr_matrix)
  
  # if(length(safe_AllHall) == 0) {
  #   print("no module to test")
  #   return()
  #   }
  
  # Step 5 : Run rROMA
  try({
    
      rRoma.output <- rRoma::rRoma.R(
        ExpressionMatrix = expr_matrix, 
        ModuleList = AllHall, #bug 33 #safe_AllHall
        #PCADims = 2, #default
        #MinGenes = 10, #default
        #SampleFilter = FALSE,
        SuppressWarning = TRUE
        )
    
    # output matrix table
    rRomaMatrix <- rRoma.output$ModuleMatrix %>% as_tibble(rownames = NA)
    rRomaMatrix <- rownames_to_column(rRomaMatrix, var="Module") %>% as_tibble()
    write_tsv(x = rRomaMatrix, file = str_c(output_path, "_rRomaMatrix.tsv"))
    
    # 3.4.1 Shifted Modules
    shifted.modules <- which(rRoma.output$ModuleMatrix[, "ppv Median Exp"] <= 0.05)
    if(!is_empty(shifted.modules)) {
      png(filename = str_c(output_path,"_shifted.png"), width = 1600, height = 800, units = "px", pointsize = 14)
      rRoma::Plot.Genesets.Samples(rRoma.output, Selected = shifted.modules,
                          GroupInfo = group, cluster_cols = T) %>% print()
      dev.off()
    
    # rRoma::Plot.Genesets.vs.Sampled(rRoma.output, Selected = shifted.modules, Plot = "PC1Median")
    }
    
    # 3.4.2 Overdispersed Modules
    overdispersed.modules <- which(rRoma.output$ModuleMatrix[, "ppv L1"] <= 0.05 & rRoma.output$ModuleMatrix[, "ppv Median Exp"] > 0.05)
    
    if(!is_empty(overdispersed.modules)) {
      png(filename = str_c(output_path, "_overdispersed_non_shifted", ".png"), width = 1600, height = 800, units = "px", pointsize = 14)
      rRoma::Plot.Genesets.Samples(rRoma.output, Selected = overdispersed.modules, 
                            GroupInfo = group, cluster_cols = TRUE) %>% print()
      dev.off()
    }
    
    # Get top genes for overdispersed and shifted modules
    # GeneMat <- rRoma::GetTopContrib(rRoma.output,
    #                          Selected = overdispersed.modules,
    #                          nGenes = 0.1, OrderType = "Abs", Mode = "Wei", Plot = FALSE)
    # topGenes = GeneMat$Table
    # topGenes %>% head(.,10)
    
    # GlobalCompareAcrossSamples(RomaData = rRoma.output, Selected = c(overdispersed.modules, shifted.modules),
    #                            Plot = "both", Groups = group)
    
    
  })
}

for (celltype in celltype_list) {
  rRoma_analysis(dataset, celltype)
}
