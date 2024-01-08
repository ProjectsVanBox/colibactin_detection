# Nissle specific functions:
plot_dinuc_profile = function(context, dinuc = "AA", order) {
  context$select = ifelse(context$pos34 == dinuc, dinuc, "other")
  context$select = factor(context$select, levels = c("other", dinuc))
  
  counts = context %>%  dplyr::count(select,
                                         type,
                                         name,
                                         trinucleotide,
                                         .drop = FALSE,
                                         name = "n_sbs") %>% as_tibble()
  rel_counts =  pivot_wider(counts, names_from = name, values_from = n_sbs)
  
  
  if (length(unique(counts$name)) > 1) {
    rel_counts[-1:-3] =  prop.table(rel_counts[-1:-3] %>% as.matrix(), 2)
    rel_counts_long = pivot_longer(
      rel_counts,
      cols = -1:-3,
      names_to = "name",
      values_to = "n_sbs"
    ) 
    rel_counts_long = rel_counts_long %>% 
      mutate(name = factor(name, levels = order))
  }
  
  else if (length(unique(counts$name)) == 1) {
    rel_counts[4] =  prop.table(rel_counts[4] %>% as.matrix(), 2)
    rel_counts_long = pivot_longer(
      rel_counts,
      cols = 4,
      names_to = "name",
      values_to = "n_sbs"
    ) 
    rel_counts_long = rel_counts_long %>% 
      mutate(name = factor(name, levels = order))
  }
  
    profile_figure = ggplot(rel_counts_long) + 
      geom_bar(aes( x = trinucleotide,
                    y = n_sbs,
                    fill = type, alpha = select),
               stat = "identity") +
      
      facet_grid(name ~ . , scales = "free_y") + 
      theme_BM() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
      scale_alpha_manual(values = c(0.3, 1)) +
      scale_fill_manual(values = COLORS6[4:6]) + 
      ggtitle(dinuc) +
      theme(legend.position = "none", axis.text.x = element_blank()) +
      xlab("") + ylab("Relative Contribution")
  return(profile_figure)
}


get_dinuc_enrichment = function(dinuc_mat) {
  
  ftest = vector("numeric", nrow(dinuc_mat))
  names(ftest) = rownames(dinuc_mat)
  index_ftest = which(rowSums(dinuc_mat) != 0)
  
  ftest[rowSums(dinuc_mat) == 0] = 1  # samples with no exonic mutations. P-value = 1 
  
  
  for (i in index_ftest)  {
    print(i)
    test_counts = c(dinuc_mat[i, 1], sum(dinuc_mat[i, -1]))
    ctrl_counts = c(sum(dinuc_mat[-i, 1]), sum(dinuc_mat[-i, -1]))
    mat = rbind(test_counts, ctrl_counts)
    ftest[[i]] = fisher.test(mat, alternative = "greater")["p.value"][[1]]
  }
  
  # adjust p-value with fd correction
  ftest_adj_log =ftest %>% 
    p.adjust(method = "fdr") %>% 
    log10()*-1
  
  ftest_adj_log[is.infinite(ftest_adj_log)] = 300 # set al infinite values to 300 (max output)
  
  return(ftest_adj_log)
}



