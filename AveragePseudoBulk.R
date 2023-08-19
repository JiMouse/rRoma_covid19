exportAveragePseudoBulk = function (df, patient_groups_df, celltype, sample_id, group, subdir, dataset) {
  df <- df %>% as.matrix() %>% as_tibble(rownames="gene")
  patient_groups_df <- patient_groups_df %>% select(all_of(sample_id), all_of(group))
  group_cond = patient_groups_df[[group]] %>% unique() #list of groups
  genes = df$gene
  
  # select groups
  g = list()
  selected_cond = c()
  
  for(x in 1:length(group_cond)) {
    g[[x]] = patient_groups_df[[sample_id]][patient_groups_df[[group]] == group_cond[x]]
    selected_col = g[[x]][g[[x]] %in% names(df)]
    
    if(!rlang::is_empty(selected_col)) {
      df[[group_cond[x]]] = rowMeans(select(df, all_of(selected_col)), na.rm = TRUE)
      selected_cond = c(selected_cond, group_cond[x])
    }
  }
  df <- df %>% select(all_of(selected_cond))

  # normalize
  dge <- DGEList(df)
  dge <- calcNormFactors(dge, method = "TMM")
  cpm <- edgeR::cpm(dge)
  cpm <- cpm %>% as.matrix() %>% as_tibble() %>%
    mutate(gene = genes) %>%
    select(gene, everything())
  
  outputname = paste0(dataset, "_", celltype, "_pseudobulk_average.tsv")
  write_tsv(x = cpm, file = outputname)
}
