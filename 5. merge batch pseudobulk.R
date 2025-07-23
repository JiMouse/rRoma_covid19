library(edgeR)
library(limma)
library(tidyverse)

# Define paths
path = "/Users/jean-marie/Documents/Privé/Curie/covid19" #"C:/Users/Administrateur/Documents/covid19/" #"/Users/jean-marie/Documents/Privé/Curie/covid19"
setwd(dir = path)
base_dir <- "/Users/jean-marie/Documents/Privé/Curie/covid19/pseudobulk"
data_dir <- "/Users/jean-marie/Documents/Privé/Curie/covid19/data"

#############
#12436517 : nasopharyngeal protected specimen brush (bronchial lavage available)
#GSE161918 : PBMC
#GSE176269 : nasal wash (fluA)
#GSE216020 : PBMC
#GSE149689 : PBMC (flu A)
#############

# List all datasets
datasets <- list.dirs(base_dir, full.names = FALSE, recursive = FALSE)#[-c("Archived_pseudobulk")]
datasets <- datasets[!(datasets %in% c("Archived_pseudobulk", "aggregated"))]
#datasets <- c("12436517", "GSE176269")
datasets <- c("GSE161918", "GSE216020", "GSE149689") 

# Initialize a list to store cell type data
celltype_list <- list()

for (dataset in datasets) {
  dataset_path <- file.path(base_dir, dataset)
  files <- list.files(dataset_path, pattern = "_raw_pseudobulk.tsv", full.names = TRUE)
  
  # Load patient group info
  patient_info_path <- file.path(data_dir, dataset, "patient_group.tsv")
  patient_group <- read_tsv(patient_info_path, show_col_types = FALSE)
  patient_group <- patient_group %>% 
    mutate(patient=patient %>% str_replace_all("-","_"))
  patient_group$status[patient_group$status %in% c("control", "healthy", "Healthy", "HC")] <- "control"
  patient_group$status[patient_group$status %in% c("critical", "Critical", "Severe")] <- "severe"
  patient_group$status[patient_group$status %in% c("mild","Mild", "moderate", "Moderate")] <- "moderate"
  
  for (file in files) {
    # Extract cell type from filename
    cell_type <- sub("_raw_pseudobulk.tsv", "", basename(file))
    
    # Read the counts
    counts <- read.table(file, header = TRUE, row.names = 1, sep = "\t")

    # # Harmonize cell type names (this is where you’ll manually customize)
    # cell_type <- str_replace_all(cell_type, pattern = " ", "-")
    # if (cell_type %in% c("Macrophage1", "Macrophage2")) {
    #   cell_type <- "Ma"
    # }
    
    # Add dataset ID to column names to make unique
    colnames(counts) <- paste0(dataset, "_", colnames(counts)) %>% str_replace_all("-|\\.","_")
    
    # Store counts per cell type
    if (!cell_type %in% names(celltype_list)) {
      celltype_list[[cell_type]] <- list()
    }
    
    celltype_list[[cell_type]][[dataset]] <- list(
      counts = counts,
      info = patient_group
    )
  }
}

# names(celltype_list[["B-cell"]][["12436517"]])
# celltype_list[["B-cell"]][["12436517"]]$counts %>% names()

# Final output list
normalized_corrected_counts <- list()

# Process each cell type
cell_type = names(celltype_list)[1] #test

for (cell_type in names(celltype_list)) {
  message("Processing cell type: ", cell_type)
  
  # Merge all counts by gene name
  merged_counts <- Reduce(function(df1, df2) {
    # Ensure gene name is not already a column
    if (!"gene" %in% colnames(df1)) df1 <- df1 %>% rownames_to_column("gene")
    if (!"gene" %in% colnames(df2)) df2 <- df2 %>% rownames_to_column("gene")
    
    full_join(df1, df2, by = "gene")
  }, lapply(celltype_list[[cell_type]], function(x) x$counts))
  
  # Restore rownames and fill NAs
  rownames(merged_counts) <- merged_counts$gene
  counts_all <- merged_counts[, -1]
  counts_all[is.na(counts_all)] <- 0
  
  # Combine patient info
  sample_names <- colnames(counts_all)
  patient_infos <- do.call(rbind, lapply(names(celltype_list[[cell_type]]), function(ds) {
    info <- celltype_list[[cell_type]][[ds]]$info
    info$sample <- paste0(ds, "_", info$patient)
    info$dataset <- ds
    info
  }))
  
  # Match and filter
  patient_infos <- patient_infos[match(sample_names, patient_infos$sample), ]
  stopifnot(all(patient_infos$sample == colnames(counts_all)))
  
  # edgeR normalization
  dge <- DGEList(counts = counts_all, group = patient_infos$status)
  dge <- calcNormFactors(dge, method = "TMM")
  
  if (length(unique(patient_infos$status)) > 1 && nrow(patient_infos) > 2) {
    # Voom with design matrix
    design <- model.matrix(~ status, data = patient_infos)
    v <- voom(dge, design)
    
    # Batch correction only if more than one dataset
    # Add the 'design = design' argument to protect the biological signal
    if (length(unique(patient_infos$dataset)) > 1) {
      v_corrected <- removeBatchEffect(v$E, batch = patient_infos$dataset, design = design)
    } else {
      message("Only one dataset for cell type ", cell_type, " — skipping batch correction.")
      v_corrected <- v$E
    }
    
  } else {
    message("Skipping voom for cell type ", cell_type, 
            ": too few samples or only one condition level.")
    
    # Fallback normalization: log2(CPM + 1)
    v_corrected <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
  }
  
  # Store normalized matrix
  normalized_corrected_counts[[cell_type]] <- v_corrected
  
}

#Save output
lapply(names(normalized_corrected_counts), function(ct) {
  write.table(normalized_corrected_counts[[ct]], file = paste0("pseudobulk/aggregated/", ct, "_pseudobulk.tsv"), sep = "\t", quote = FALSE)
})

############
# Collect all patient group information in one dataframe
patient_group_full <- do.call(rbind, lapply(names(celltype_list), function(cell_type) {
  lapply(names(celltype_list[[cell_type]]), function(ds) {
    df <- celltype_list[[cell_type]][[ds]]$info
    df$sample <- df$patient
    df$patient <- paste0(ds, "_", df$sample)
    df$dataset <- ds
    df
  }) %>% bind_rows()
})) %>% distinct(sample, .keep_all = TRUE)
patient_group_full$status %>% table()

# Export full patient group in data
write_tsv(x = patient_group_full, file = str_c(data_dir, "aggregated/patient_group.tsv", sep="/"))

