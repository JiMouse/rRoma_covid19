rm(list=ls()) ; cat("\014") ; gc() ; graphics.off()

setwd(dir = "C:/Users/Administrateur/Documents/PhD/covid19/")
wd = getwd()
set.seed(2016)

library(tidyverse)
library(rRoma)
# https://sysbio-curie.github.io/rROMA/index.html


dataset_list = list( #GSE210694
  c(dataset = "12436517" , barcodes_name = "patient" , group_name = "severity")
  ,c(dataset = "GSE150728" , barcodes_name = "Donor.orig" , group_name = "Status")
  ,c(dataset = "GSE161918" , barcodes_name = "Subject" , group_name = "cond_group")
  ,c(dataset = "GSE171524" , barcodes_name = "donor_id" , group_name = "group")
  ,c(dataset = "GSE174072" , barcodes_name = "Donor" , group_name = "Status")
  ,c(dataset = "GSE176269", barcodes_name = "sampID", group_name = "status")
  ,c(dataset = "GSE176269_without_fluA", barcodes_name = "sampID", group_name = "status")
  ,c(dataset = "GSE149689", barcodes_name = "patient", group_name = "status")
  ,c(dataset = "GSE157789", barcodes_name = "patient", group_name = "status")
  ,c(dataset = "GSE210694", barcodes_name = "patient", group_name = "status")
)
length(dataset_list)
# rRoma one dataset
d=1 #1 6 7 8
dataset=dataset_list[[d]][["dataset"]]
barcodes_name=dataset_list[[d]][["barcodes_name"]]
group_name=dataset_list[[d]][["group_name"]]

module_name="AllHall" #AllHall #REACTOME_SARS_COV_2_GMT.gmt
# module_name="ACSN"

subdir = paste0("_Navicell/Datasets/",dataset,"/output/")
plotScores=TRUE ; exportTSV=TRUE ; includShifted = FALSE ; module=module_name
pseudobulk = TRUE ; p = "_pseudobulk.tsv"
celltype_list = list.files(path=subdir, pattern = p) %>% 
  str_remove_all(str_c(dataset, "_")) %>%
  str_remove_all(p)
print(celltype_list)

# celltype_list = "nrMa"
# celltype = celltype_list

for(celltype in celltype_list[1:length(celltype_list)]) {

  ########################## Analysis #############
  graphics.off()
  datafile = paste0(dataset, "_", celltype, "_pseudobulk")
  export_name = str_c(dataset, "_", celltype, "_", module, "_pseudobulk")
  subdir_output = str_c("/_Navicell/Datasets/", dataset, "/rRoma/")
  
  # Loading dataset
  counts <- read_tsv(file = paste0(subdir, datafile, ".tsv"), show_col_types = FALSE)
  counts <- column_to_rownames(counts, "gene")
  
  counts_average <- read_tsv(file = paste0(subdir, datafile, "_average.tsv"), show_col_types = FALSE)
  # load as character
  
  # loading groups
  patient_group = "patient_group.tsv"
  sample_labels = read_tsv(file = paste0(subdir, patient_group), show_col_types = FALSE) #%>% arrange(group)
  group <- sample_labels[[group_name]]
  names(group) <- sample_labels[[barcodes_name]] #paste0(sample_labels[[barcodes_name]], "_", celltype)
  # group <- group %>% sort()
  group <-group[names(group) %in% names(counts)]

  # 3.2 Selecting the module list : Load a GMT file from a file
  subdir_gmt = "_Navicell/GMT/ACSN/"
  if(module == "ACSN") {
    # gmt_file_select = c("NK" = "natural_killer_cell_master.gmt",
    #                     "Ma" = "macrophages_mdsc_cells_master.gmt",
    #                     "DC" = "dendritic_cell_master.gmt",
    #                     "LT" = "adaptive_immune_master.gmt")
    # gmt_file = gmt_file_select[celltype]
    gmt_file = "rcd_master.gmt"
    modules <- ReadGMTFile(FileLocation = paste0(subdir_gmt, gmt_file))
    modules_name = gmt_file
  } else if(module == "AllHall"){
    AllHall <- SelectFromMSIGdb("HALLMARK")
    AllHall <- lapply(AllHall, function(x){
      x$Name <- sub("HALLMARK_", "", x$Name)
      x
    })
    modules = AllHall
    modules_name = "AllHall"
  } else {
    dir_gmt = str_c("C:/Users/Administrateur/Documents/PhD/covid19/_Navicell/GMT/", module)
    modules <- ReadGMTFile(FileLocation = dir_gmt)
    modules_name = module
  }
  
  ##### performing rROMA ##########
  ifrm <- function(obj, env = globalenv()) {
    obj <- deparse(substitute(obj))
    if(exists(obj, envir = env)) {
      rm(list = obj, envir = env)
    }
  }
  
  ifrm(rRoma.output) #clean object
  try(rRoma.output <- rRoma.R(ExpressionMatrix = counts, ModuleList = modules,
                          # GeneOutDetection = "L1OutExpOut",
                          # GeneOutThr = 5,
                          # GeneSelMode = "All",
                          # SampleFilter = F,
                          PCSignMode = 'none', #default
                          # PCSignMode = "UseMeanExpressionAllWeights" #from Loredana
                          # PCSignMode = "UseMeanExpressionKnownWeights"
                          ), silent=FALSE)
  
  
  if(!exists("rRoma.output")) next
  if(is_empty(rRoma.output)) next
  
  # create new folder
  subdir_output = str_c(wd,subdir_output, modules_name,"_", celltype, "/")
  dir.create(subdir_output, showWarnings = FALSE)
  
  ##### export graphs ############
  
  outputFilename = paste0(subdir_output, "rROMa.output_", export_name)
  rRomaMatrix <- rRoma.output$ModuleMatrix %>% as_tibble(rownames = NA)
  rRomaMatrix <- rownames_to_column(rRomaMatrix, var="Module") %>% as_tibble()
  if(exportTSV) {write_tsv(x = rRomaMatrix, file = paste0(outputFilename, ".tsv"))}
  
  ##grouped
  # png(filename = paste0(subdir_output, export_name,"_grouped.png"), width = 1200, height = 800, units = "px", pointsize = 14)
  # rRoma.output_grouped <- Plot.Genesets.Samples(RomaData = rRoma.output,
  #               Selected = SelectGeneSets(RomaData = rRoma.output, VarThr = 5e-2,
  #                                         VarMode = "Wil", VarType = "Over"),
  #               GenesetMargin = 20, SampleMargin = 14, cluster_cols = FALSE,
  #               GroupInfo = group, AggByGroupsFL = c("mean"),
  #               HMTite = "Overdispersed genesets (Fixed center)")
  # dev.off()
  
  ##all
  png(filename = paste0(subdir_output, export_name,"_all.png"), width = 1200, height = 800, units = "px", pointsize = 14)
  Plot.Genesets.Samples(rRoma.output, GroupInfo = group, cluster_cols = TRUE)
  dev.off()
  
  # rRomaMatrix_grouped <- rRoma.output_grouped$mean %>% as_tibble(rownames = NA)
  # rRomaMatrix_grouped <- rownames_to_column(rRomaMatrix_grouped, var="Module") %>% as_tibble()
  # write_tsv(x = rRomaMatrix_grouped, file = paste0(outputFilename, "_grouped.tsv"))
  
  # 3.4.1 Shifted Modules
  shifted.modules <- which(rRoma.output$ModuleMatrix[, "ppv Median Exp"] <= 0.05)
  if(!is_empty(shifted.modules)) {
    png(filename = paste0(subdir_output, export_name,"_shifted.png"), width = 1200, height = 800, units = "px", pointsize = 14)
    Plot.Genesets.Samples(rRoma.output, Selected = shifted.modules,
                          GroupInfo = group, cluster_cols = TRUE)
    dev.off()
    
    GeneMat_shifted <- GetTopContrib(rRoma.output,
                             Selected = shifted.modules, #c(overdispersed.modules, shifted.modules), #
                             nGenes = 0.1, OrderType = "Abs", Mode = "Wei", Plot = FALSE)
    # topGenes = rownames_to_column(GeneMat$Table, var = "rank")
    topGenes_shifted = GeneMat_shifted$Table
    
    #left_join with average expression
    topGenes_shifted <- left_join(topGenes_shifted, counts_average, by=c("Gene"="gene")) %>%
      mutate(across(where(is.numeric), round, 3))
    write_tsv(x = topGenes_shifted, file = paste0(subdir_output, "rROMa_GetTopContrib_", export_name, "shifted", ".tsv"))
  }
  # 3.4.2 Overdispersed Modules
  if(includShifted == FALSE) {
    overdispersed.modules <- which(rRoma.output$ModuleMatrix[, "ppv L1"] <= 0.1 & rRoma.output$ModuleMatrix[, "ppv Median Exp"] > 0.05)
    tag = "_overdispersed_non_shifted"
    } else {
    overdispersed.modules <- which(rRoma.output$ModuleMatrix[, "ppv L1"] <= 0.1)
    tag = "_overdispersed"
    }
  if(!is_empty(overdispersed.modules)) {
  
    if(plotScores) {
      png(filename = paste0(subdir_output, export_name, tag, ".png"), width = 1200, height = 800, units = "px", pointsize = 14)
      Plot.Genesets.Samples(rRoma.output, Selected = overdispersed.modules,
                            GroupInfo = group, cluster_cols = TRUE)
      dev.off()
    }
  
    # Get top genes for overdispersed and shifted modules
    ifrm(GeneMat)
    try(GeneMat <- GetTopContrib(rRoma.output,
                                 Selected = overdispersed.modules,
                                 nGenes = 0.1, OrderType = "Abs", Mode = "Wei", Plot = FALSE),
          silent=TRUE)
    if(exists("GeneMat")) {
      topGenes = GeneMat$Table
    
      #left_join with average expression
      topGenes <- left_join(topGenes, counts_average, by=c("Gene"="gene")) %>%
        mutate(across(where(is.numeric), round, 3))
        # mutate(across(\(x) where(is.numeric), round, 3))
      
      write_tsv(x = topGenes, file = paste0(subdir_output, "rROMa_GetTopContrib_", export_name, tag, ".tsv"))
    }
    print(paste0("celltype : ", celltype, " ; module = ", module))
  
    # 3.4.3.3 Looking at the orientation of the PCs
    # m = shifted.modules #c(overdispersed.modules, shifted.modules)
    # i=1
    # for(i in 1:length(m)) {
      # print(m[i])
    #   png(filename = paste0(subdir_output, export_name, "_", names(m[i]), "_GeneWeight.png"))#, width = 1600, height = 800, units = "px", pointsize = 14)
      # PlotGeneWeight(RomaData = rRoma.output, PlotGenes = 30,
      #                ExpressionMatrix = counts, LogExpression = FALSE,
      #                Selected = m[i])
      # PlotSampleProjections(RomaData = rRoma.output, PlotSamples = 30,
      #                       ExpressionMatrix = counts, LogExpression = FALSE,
      #                       Selected = m[i])
    #   dev.off()
    # }
    } else {
    print("no overdispersed modules")
    }
  # print(subdir_output)

}
