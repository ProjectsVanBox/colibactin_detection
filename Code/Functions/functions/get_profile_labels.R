# function to make the script a bit shorter: 

get_profile_labels = function(ext_context, TN_signature) { 
  
  context_counts = ext_context %>% 
    group_by(name, select, trinucleotide) %>% 
    dplyr::count() %>% 
    pivot_wider(values_from = n, names_from = trinucleotide, values_fill = 0) %>% 
    ungroup() %>% as.data.table()
  
  # get the trinucleotide profile for the total samples: 
  type_counts_all = context_counts[, lapply(.SD, sum), by = c("name"),
                                   .SDcols = names(context_counts) %like% "[.>.]"]
  type_counts_all = type_counts_all %>% 
    column_to_rownames("name") %>% 
    t() %>% as.data.frame()
  type_counts_all = type_counts_all[TRIPLETS_48,]
  cosines_TN = sapply(type_counts_all, cos_sim, TN_signature)
  spearman_TN = sapply(type_counts_all, cor, TN_signature, method = "spearman") # get cosine similarities to SBS T>N fraction
  spearman_TN_pval = sapply(type_counts_all, \(x) cor.test(x, TN_signature, method = "spearman")$p.value) %>% 
    p.adjust(method = 'fdr') %>% 
    formatC(format = "e", digits = 2)# get spearman correlation
  
  
  TN_trinuc_counts = context_counts %>% filter(select == "AA") %>% 
    dplyr::select(-select) %>% column_to_rownames("name") %>% 
    t() %>% as.data.frame()
  TN_trinuc_counts = TN_trinuc_counts[TRIPLETS_48,]
  
  cosines_TN_AA = sapply(TN_trinuc_counts, cos_sim, TN_signature)
  spearman_TN_AA = sapply(TN_trinuc_counts, cor, TN_signature, method = "spearman") 
  spearman_TN_AA_pval = sapply(TN_trinuc_counts, \(x) cor.test(x, TN_signature, method = "spearman")$p.value) %>% 
    p.adjust(method = 'fdr') %>% 
    formatC(format = "e", digits = 2)# get cosine similarities to SBS T>N fraction
  
  
  label_df = data.frame(name = factor(names(cosines_TN), levels = names(cosines_TN)), 
                        cosine_TN = cosines_TN,
                        spearman_TN = spearman_TN,
                        cosine_TN_AA = cosines_TN_AA,
                        spearman_TN_AA = spearman_TN_AA, 
                        label_cosine = paste0("all T>N: ", round(cosines_TN, 2),
                                              "\nT>N -3-4AA: ", round(cosines_TN_AA, 2)),
                        label_spearman = paste0(round(spearman_TN, 2) ,
                                       "\n", round(spearman_TN_AA, 2)), 
                        label_pval = paste0(spearman_TN_pval, "\n", spearman_TN_AA_pval))
  
  return(label_df)
}

