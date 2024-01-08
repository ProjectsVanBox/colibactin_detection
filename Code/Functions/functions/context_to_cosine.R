# calculate cosine similarites from a list of occurrences: 

context_to_cosine = function(list, profile){ 
  
  ctx_list = rbindlist(list, idcol = "sample")
  
  # get cosine for full profile: 
  counts_all = ctx_list %>% 
    group_by(sample, trinucleotide) %>% 
    dplyr::count() 
  matrix_all = counts_all %>% 
    pivot_wider(names_from = trinucleotide, values_from = n, values_fill = 0, names_expand = TRUE) %>%
    column_to_rownames("sample") %>% dplyr::select(all_of(TRIPLETS_48)) %>% 
    as.matrix()
  
    # get cosine for AA nucleotides 
  counts_AA = ctx_list %>%
    filter(select == "AA") %>% 
    group_by(sample, trinucleotide) %>% 
    dplyr::count() 
  matrix_AA = counts_AA %>% 
    pivot_wider(names_from = trinucleotide, values_from = n, values_fill = 0, names_expand = TRUE) %>%
    column_to_rownames("sample") %>% dplyr::select(all_of(TRIPLETS_48)) %>% 
    as.matrix()
  cs_all = apply(matrix_all, 1, MutationalPatterns::cos_sim, profile)
  cs_AA = apply(matrix_AA, 1, MutationalPatterns::cos_sim, profile)
  
  
  fraction_AA = rowSums(matrix_AA) / rowSums(matrix_all)
  
  df = data.table(id = names(cs_all), 
                  `all T>N mutations` = cs_all,
                  `T>N mutations -3-4AA` = cs_AA, 
                  fraction_AA = fraction_AA)
  return(df)
  
  
}
