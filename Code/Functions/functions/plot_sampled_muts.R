# Boostrap_muts

plot_sampled_muts = function(mut_list, name = "none") {
  set.seed(123456)
  
  n_EcN = sum(mut_list$name == "EcN") # get the number of EcN muts to sample in the specific mutation list
  
  list = lapply(1:1000, \(x)  {
    tmp = mut_list %>% 
      group_by(name) %>% 
      slice_sample(n = n_EcN, replace = TRUE)
    tmp$bin = x
    return(tmp)})
  
  total_list = rbindlist(list) %>% 
    group_by(select, name, bin) %>% 
    dplyr::count()
  
  TN_muts = ggplot(total_list %>% filter(select == "AA"), aes(y = n, x = name, fill = name)) +
    geom_boxplot(outlier.shape = NA)  + 
    geom_jitter(size = 0.5, width = 0.05, alpha = 0.6) + 
    xlab("") + ylab("Count T>N mutations with -3-4AA ") +
    ggpubr::stat_pwc() + theme_bw() + labs(fill = "Exposure\ncondition") + 
    theme(legend.position = "none") + scale_y_continuous(limits = c(0, NA)) 
  
  if (name != "none") {
    
    TN_muts = TN_muts + ggtitle(name) + theme(plot.title =  element_text(hjust =0.5) )
  
  }
  
  return(TN_muts)
}
