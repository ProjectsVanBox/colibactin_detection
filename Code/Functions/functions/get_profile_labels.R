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
  cosines_TN = sapply(type_counts_all, cos_sim, TN_signature) # get cosine similarities to SBS T>N fraction
  
  TN_trinuc_counts = context_counts %>% filter(select == "AA") %>% 
    dplyr::select(-select) %>% column_to_rownames("name") %>% 
    t() %>% as.data.frame()
  TN_trinuc_counts = TN_trinuc_counts[TRIPLETS_48,]
  cosines_TN_AA = sapply(TN_trinuc_counts, cos_sim, TN_signature) %>% 
    as.data.frame() %>% rownames_to_column("name")
  colnames(cosines_TN_AA)[2] = "cs"
  
  label_df = data.frame(name = factor(cosines_TN_AA$name, levels = cosines_TN_AA$name), 
                        cosine = paste0("all T>N: ",
                                       round(cosines_TN, 2), 
                                       "\nT>N -3-4AA: ", 
                                       round(cosines_TN_AA$cs, 2)))
  
  return(label_df)
}

