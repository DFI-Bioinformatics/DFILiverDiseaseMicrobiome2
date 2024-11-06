function_compute_alpha_div = function(df_mp_species) {
  tmatm = df_mp_species %>%
    select(seq_id, species_mp, pctseqs) %>%
    pivot_wider(
      names_from = species_mp,
      values_from = pctseqs,
      values_fn = sum,
      values_fill = 0
    ) %>%
    column_to_rownames(var = "seq_id")
  
  tmatm.tran <- round(tmatm / rowSums(tmatm) * 1000000)        
  invm <- vegan::diversity(tmatm, index = "invsimpson")   
  
  invdfm <- tibble(seq_id = names(invm), invSimp =  invm)          
  
  if (length(invm) < 2) {
      invdfm <- invdfm %>%         
        mutate(seq_id = rownames(tmatm))
                        }        
  
  alpha.m <- estimateR(tmatm.tran) %>%     
    t() %>%     
    as.data.frame() %>%     
    rownames_to_column(var = "seq_id") %>%     
    left_join(invdfm)   
  
  return(alpha.m)
}