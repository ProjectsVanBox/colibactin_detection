# Figure 2 function
# function to plot the analyses for figure 3 two times: For PTA-exposed samples and for Clonal expansion

# USAGE: 
# 1. Clone the "https://github.com/ProjectsVanBox/colibactin_detection" repository or download as .zip file
# Set the working directory in the 'setwd()' function as the working directory of this script in the following line and you should be good to go:
setwd("C:/Users/Axel Rosendahl Huber/OneDrive/Nissle_manuscript/Nissle/")


library(ggseqlogo)
library(ggplot2)
library(gtools)
library(ggh4x)
library(patchwork)
library(ggridges)
library(dtplyr)
library(tidyverse)

# source function and data loading scripts: 
source("Code/Load_data.R")
source("Code/Functions/Utils.R")
source("Code/Functions/Nissle_functions.R")

contexts_TN = list()
for (type in names(context_list)) {
  ctx_table= context_list[[type]] %>%
    rbindlist() %>%
    distinct() %>%
    filter(grepl("^T", type))
  contexts_TN[[type]] = ctx_table
}

# select only the PTA-samples
cat_PTA = categories %>%
  filter(method == "PTA")

# select only the clonal expansion-generated samples
cat_CE = categories %>% 
  filter(method == "Clonal Expansion")
cat_CE$injection = factor(cat_CE$injection, levels = c("Control", "EcN", "EcC"))

contexts_TN_PTA = list()
for (type in unique(cat_PTA$injection)) {
  ctx_table= context_list[[type]] %>%
    rbindlist(idcol = "samplename") %>%
    filter(samplename %in% cat_PTA$name) %>% 
    dplyr::select(-samplename) %>% 
    distinct() %>%
    filter(grepl("^T", type))
  contexts_TN_PTA[[type]] = ctx_table
}

contexts_TN_CE = list()
for (type in unique(cat_CE$injection)) {
  ctx_table= context_list[[type]] %>%
    rbindlist(idcol = "samplename") %>%
    filter(samplename %in% cat_CE$name) %>% 
    dplyr::select(-samplename) %>% 
    distinct() %>%
    filter(grepl("^T", type))
  contexts_TN_CE[[type]] = ctx_table
}

# testing variables only - uncomment if you want to test the individual lines in the plotting function
# TN_contexts = contexts_TN_PTA
# cat = cat_PTA
# name = 'PTA'

TRIPLETS_48 = TRIPLETS_96[49:96]
SBS88_TN <- as.data.table(signatures) %>% dplyr::slice(49:96) %>% pull("SBS88")

plot_figures_2 = function(TN_contexts, cat, name) {
  
  ####### -3 -4 2bp upstream motif 
  ext_context = rbindlist(TN_contexts, idcol = "name") %>% 
    mutate(pos34 =  substr(context, 7,8)) %>% 
    mutate(trinucleotide = factor(trinucleotide, levels = TRIPLETS_48)) %>% 
    mutate(select =  ifelse(pos34 == "AA", "AA", "other") %>% 
             factor(levels = c("other", "AA"))) %>% 
    mutate(name = name %>% factor(levels = levels(cat$injection)))
  
  # 1. compare to the total level of trinucleotides
  label_df = get_profile_labels(ext_context, SBS88_TN)
  # Add the explaining label on the plot:
  label_df$label_cosine = paste0("SBS88:    cosine similarity\n", label_df$label_cosine)
  label_df$label_spearman = paste0("spearman\n", label_df$label_spearman)
  label_df$label_pval = paste0("pval\n", label_df$label_pval)
  
  plot_profile_absolute = function(mut_list) {
    
    mut_list %>% 
      mutate(select =  fct_recode(select, "-3-4AA" = "AA")) %>% 
      group_by(trinucleotide, select) %>% 
      ggplot(aes(x = trinucleotide, alpha = select, fill = type)) +
      geom_bar(width = 0.7 ) +
      facet_grid(name ~ . , scales  = "free_y") + 
      theme_BM()  + 
      scale_alpha_manual(values = c(0.3, 1)) +
      scale_fill_manual(values = COLORS6[4:6]) + 
      theme(legend.position = "top",
            legend.box = "horizontal", 
            axis.text.x = element_text(size= 6.5, angle = 90, vjust = 0.5, hjust=1),
            legend.box.background = element_rect(colour = "black", fill = NA), 
            legend.background =element_blank(),
            legend.text = element_text(size=7),
            legend.title = element_blank(),
            strip.background = element_blank(), 
            strip.text.y = element_blank(),
            legend.key.size = unit(7, "points"),
            plot.margin = margin(unit(c(3, 8, 8, 8), "points"))) +
      xlab("") + ylab("Mutation count")
    
  }
  
  F3e_AA_context_profile = plot_profile_absolute(ext_context) +
    labs(alpha = "nucs at pos -3-4", fill = NULL) + 
    ggpp::geom_text_npc(data = label_df, npcx = 0.98, npcy = 0.9,
                        aes(label = label_pval), size = 2.5, hjust = 1) +
    ggpp::geom_text_npc(data = label_df, npcx = 0.88, npcy = 0.9,
                        aes(label = label_spearman), size = 2.5, hjust = 1) +
    ggpp::geom_text_npc(data = label_df, npcx = 0.78, npcy = 0.9,
                       aes(label = label_cosine), size = 2.5, hjust = 1) + 
      ggpp::geom_text_npc(data = label_df, npcx = 0.02, npcy = 0.9,
                        aes(label = name), size = 3.5, hjust = 0)
    
  F3e_AA_context_profile
  
  # 1. all muts
  EcN_Ctrl_muts = ext_context %>% 
    filter(name %in% c("Control", "EcN"))
  TN_muts = plot_sampled_muts(EcN_Ctrl_muts)

  # 'Monte Carlo histogram plotting: 
  # see if it is possible to perform thousand of T-tests to empirically validate the mutagenic activity:
  names = levels(ext_context$name)[1:length(unique(ext_context$name))]
  sim_list = mapply(1:(200 * length(names)), names, FUN = \(i,n) {
    print(i)
    n_observed = sum(ext_context$name == n)
    tmp = ext_context %>% 
      filter(name == n) %>% 
      slice_sample(n = n_observed, replace = TRUE)
    tmp$bin = i
    tmp$injection = n
    return(list(tmp))
  })  
  
  sim_contexts = rbindlist(sim_list) %>% 
      group_by(select, name, bin) 
  
  sim_data = sim_contexts %>% 
      dplyr::count() %>% 
    pivot_wider(values_from = n, names_from = c(select)) %>% 
    ungroup() %>% 
    mutate(name = factor(name, levels  = c("Control", "EcN", "EcC", "19H2", "2F8")))
  
  true_mat = ext_context %>% group_by(select, name) %>% 
    dplyr::count() %>% 
    pivot_wider(names_from = select, values_from = n) %>% 
    as.data.frame() %>% 
    mutate(name = as.character(name))
  
  sim_data$pval = NA
  sim_mat = as.matrix(sim_data %>% dplyr::select(other, AA))
  ctrl_vals = true_mat[1, 2:3] %>% as.numeric()
  
  sim_data$pval = sapply(1:nrow(sim_data), FUN =  \(i) {
    m = matrix(c(ctrl_vals, sim_mat[i,]), nrow = 2, byrow = TRUE)
    fisher.test(m, alternative = "greater")$p.value
  })
  
  list_occurrences = split(ext_context, ext_context$name)
  occurrence_mat = sapply(list_occurrences, function(x) table(x$select))

  true_mat$name = factor(true_mat$name, levels = levels(sim_data$name))
  true_mat$pval = 1
  for (i in 2:nrow(true_mat)) {
    m = matrix(c(ctrl_vals, as.numeric(true_mat[i,2:3])), nrow = 2, byrow = TRUE)
    true_mat$pval[i] = fisher.test(m, alternative = "greater")$p.value
    
  }
  
  histogram_fisher = ggplot(sim_data %>% filter(name != "Control"),
                            aes(x = -log10(pval))) + 
    geom_histogram(aes(color = name, fill = name), binwidth = 0.1) + 
    geomtextpath::geom_textvline(data = true_mat %>% filter(name != "Control"),
                                 aes(xintercept = -log10(pval),
                                     label = paste0(name, " p-value: ", format(pval, digits = 3)), 
                                     linetype = "dashed"), size = 3, hjust = 1) + 
    lemon::facet_rep_grid(. ~ name, scale = "free_x") + 
    scale_color_manual(values = c("#00e8fc", "#f96e46", "#f9c846", "#ffe3e3", "#545863")) + 
    scale_fill_manual(values = c("#00e8fc", "#f96e46", "#f9c846", "#ffe3e3", "#545863")) + 
    geom_histogram(data = sim_data %>% filter(name == "Control") %>% dplyr::select(-name),
                   mapping = aes(x = -log10(pval)), binwidth = 0.1,  fill = "#545863", color = "#545863") + 
    xlab("-log10 p-value") +
    theme_BM() + theme(legend.position = "none")
  histogram_fisher
  
  # get cosine similarities to the SBS88 profile:
  cosine_list = lapply(1:length(names), \(x) {
         index = seq(0, length(sim_list)-1, length(names))
         context_to_cosine(sim_list[index + x], SBS88_TN)
  })
  names(cosine_list) = names
  cosine_df = rbindlist(cosine_list, idcol = "Exposure")
  
  cs_real_data = context_to_cosine(split(ext_context, ext_context$name), SBS88_TN) %>% 
    mutate(Exposure = id)
  
  cosines = rbindlist(list(simulations = cosine_df,
                           `real data` = cs_real_data), idcol = "type", use.names = TRUE) %>% 
    dplyr::select(-`T>N mutations -3-4AA`)
  
  cosines_long = pivot_longer(cosines, cols = -c(id,fraction_AA, type, Exposure)) %>% 
    arrange(desc(type)) %>% 
    mutate(Exposure = factor(Exposure, levels = c("Control", "EcN", "EcC", "19H2", "2F8")))
  
  simulation_plot = ggplot(cosines_long %>% filter(Exposure != "Control"), aes(
    x = value, y = fraction_AA, fill = Exposure, size = type, alpha = type)) + 
    geom_point(shape = 21, color = "black") + 
    lemon::facet_rep_grid(. ~ Exposure) + 
    geom_point(data = cosines_long %>% filter(Exposure == 'Control') %>% dplyr::select(-Exposure), 
               mapping = aes(x = value, y = fraction_AA), fill = '#545863', color = "black") + 
    scale_alpha_manual(values = c(1, 0.3)) + 
    scale_size_manual(values = c(5, 1)) + 
    scale_fill_manual(values = c("#00e8fc", "#f96e46", "#f9c846", "#ffe3e3")) + 
    theme_BM() + 
    ylab("fraction AA") + 
    xlab("cosine similarity to SBS88") + 
    ggtitle(name) + 
    theme(legend.position = "none", 
          plot.margin = margin(unit(c(5.5, 5.5, 5.5, 5.5), "points")))
  simulation_plot
  
  supplementary_figure_4 = simulation_plot

  # figure 2D - test the p-values for the differet dinucleotides
  # test p-value for enrichment
  list = split(ext_context, ext_context$name)
  occurrence_mat = sapply(list, function(x) table(x$select))
  AA_occurrence_mat_percentages = occurrence_mat[2,] / colSums(occurrence_mat)
  
  # get the enrichment scores for only A mutations at the -3 site: 
  list3 = list %>% rbindlist() %>% 
    mutate(pos3 = substr(pos34, 2,2))
  
  occurrence_mat_A = table(list3$name, list3$pos3) %>% as.data.frame.matrix()
  A_occurrence_mat_percentages = occurrence_mat_A[,1] / rowSums(occurrence_mat_A)
  AA_occurrence_mat_percentages /A_occurrence_mat_percentages
  
  # perform fisher test on all motif enrichments
  dinuc_counts = sapply(list, function(x) table(x$pos34))
  
  fisher_table = matrix(NA, nrow = 16, ncol = length(unique(cat$injection))) %>%
    `colnames<-`(unique(cat$injection)) %>% 
    as.data.frame() %>% 
    mutate(dinucs = dinucs) %>% 
    dplyr::select(-Control)
  
  for (injection in colnames(fisher_table %>% dplyr::select(-dinucs))) {
    
    fisher_table[, injection] = sapply(1:nrow(fisher_table), \(x) {
      select = dinuc_counts[x,c("Control", injection)]
      other = colSums(dinuc_counts[-x,c("Control", injection)])
      mat = rbind(other, select)
      fisher.test(mat, alternative = "greater")$p.value
    })
  }
  
  print(paste0("EcC enrichment ",name, " = ", fisher_table %>% filter(dinucs == "AA") %>% dplyr::select("EcC")))
  print(paste0("EcC enrichment ",name, " = ", fisher_table %>% filter(dinucs == "AA") %>% dplyr::select("EcN")))
  numcols = sapply(fisher_table, is.numeric)
  
  print(fisher_table %>% filter(dinucs == "AA"))
  fisher_table_m = pivot_longer(fisher_table, cols = -dinucs, values_to = "pval")
  
  # make the fisher table for the simulated data
  total_counts = table(ext_context$name) %>% 
    as.data.frame() %>% 
    `colnames<-`(c("injection", "total_muts"))
  
  simulated_counts = sim_list %>% 
    rbindlist() %>% 
    mutate(bin_index = ceiling(bin/length(names))) %>%
    group_by(bin_index, injection, pos34) %>% 
    dplyr::count() %>% 
    pivot_wider(names_from = pos34, values_from = n, values_fill = 0)
  
  sim_cnts =left_join(simulated_counts, total_counts, by = "injection")    %>% 
    ungroup() %>% 
    as.data.table()
  
  fisher_result = list()
  iterative_fisher = function(sim_counts) {
  
  strains_test = unique(sim_counts %>% filter(injection != "Control") %>% pull(injection))
  for (sel_strain in strains_test) {
      
      print(sel_strain)
      cnts_strain = sim_counts %>% 
        filter(injection %in% c(sel_strain, "Control")) %>% 
        mutate(injection = case_when(injection == sel_strain ~ "strain", 
                                     .default = injection)) %>% 
        arrange(bin_index, injection)
      
      for (sel_pos34 in unique(ext_context$pos34)) {
        
        print(sel_pos34)
        cnts_sel = cnts_strain %>% 
          dplyr::select(bin_index, injection, all_of(sel_pos34), total_muts)
        colnames(cnts_sel)[3] = "selpos"
        cnts_mat = cnts_sel %>% 
          mutate(noselect = total_muts - selpos) %>% 
          dplyr::select(-total_muts) %>% 
          pivot_wider(names_from = injection, values_from = c(selpos, noselect)) %>% 
          dplyr::select(selpos_strain, selpos_Control, noselect_strain, noselect_Control) %>% 
          as.matrix()
        
        dt_pval = data.table(bin_index = 1:200, dinucs = sel_pos34, name = sel_strain)
        dt_pval$pval = fisher_test = apply(cnts_mat, 1, \(x) 
              fisher.test(matrix(x, nrow = 2), alternative = "greater")$p.value)
        index = paste0(sel_pos34, sel_strain)
        fisher_result[[index]] = dt_pval
      }
      
    }
    return(rbindlist(fisher_result))
  }
  
  fisher_result_table = iterative_fisher(sim_cnts) %>% 
    mutate(type = "simulation")
  
  fisher_table_m = fisher_table_m %>% 
    mutate(type = "observed data")
  
  total_result = rbind(fisher_result_table %>% dplyr::select(-bin_index), fisher_table_m)
  
  F3d_dinc_enrichment = ggplot(total_result, aes(x = dinucs, y = -log10(pval), color = type)) + 
    geom_jitter(stat = "identity", width = 0.1) +
    geom_hline(yintercept = -log10(0.05), color = "grey") + 
    facet_grid(name ~ ., scales = "free_y") + 
    scale_color_manual(values = c("black", "grey")) +
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7), 
          legend.position = c(0.6, 0.95), legend.direction="horizontal", 
          legend.background = element_blank(), 
          legend.text = element_text(size = 7),
          legend.key.size = unit(2, units = "mm"),
          plot.margin = margin(unit(c(3, 3, 3, 3), "points"))) + 
    labs(x = "", y = "-log10 pvalue\nenrichment dinucleotide", color = "")
  F3d_dinc_enrichment
  
  # plot the enrichment of samples at specific motifs: 
  colnames(dinuc_counts) = paste0(colnames(dinuc_counts), "\n", colSums(dinuc_counts)," T>N SBS")
  dinuc_long = prop.table(dinuc_counts, 2) %>% 
    data.table(keep.rownames = "Dinucleotide") %>%
    mutate(color = ifelse(Dinucleotide == "AA", "AA", "noAA")) %>% 
    pivot_longer(cols = c(-Dinucleotide, -color), names_to = "Condition", values_to = "relative frequency") %>% 
    mutate(Condition = factor(Condition, levels = colnames(dinuc_counts)))
  
  F3c_dinuc_frequencies = ggplot(dinuc_long) + 
    geom_bar(aes(x = Dinucleotide, y = `relative frequency`, fill = color), stat = "identity") +
    facet_grid(. ~ Condition , scales = "free_y") + 
    theme_classic() +  panel_border()  + 
    scale_y_continuous(expand = c(0, 0), limits = c(0,0.75), breaks = seq(0, 0.75, 0.25)) + 
    scale_fill_manual(values = c("darkgreen", "gray30")) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,  size = 6 ), 
          strip.text.x = element_text(size = 7),
          legend.position = "none", panel.spacing.x = unit(1, units = "mm"), 
          plot.margin = margin(unit(c(3, 3, 3, 3), "points")))  +
    ylab("relative frequency") + xlab("")
  F3c_dinuc_frequencies
  
  # # plot sequence logo's
  TNctx = list()
  for (type in levels(cat$injection)) {
    type_context = TN_contexts[[type]]$context
    TNctx[[type]] = type_context
  }
  
  TNctx = TNctx[c("Control", "EcN", "EcC", "19H2", "2F8")]
  TNctx = TNctx[1:length(levels(cat$injection))]
  
  labels = data.frame(seq_group = factor(names(TNctx), names(TNctx)), 
                      label = paste0(names(TNctx), "\nT>N SBS:\n", lengths(TNctx)))
  
  F3a_seqlogo_plots = ggseqlogo(TNctx) + 
    annotate('rect', xmin = 10.5, xmax = 11.5, ymin = 0, ymax = 0.3, fill='grey') +  
    scale_x_continuous(breaks = c(1:21), labels= c(-10:10)) + 
    scale_y_continuous(limits = c(0,1.2), expand = c(0,0), breaks = c(0, 0.5, 1)) + 
    ggpp::geom_text_npc(data = labels, aes(label = label),
                        npcx = 0.02, npcy = 0.95, hjust = 0, size = 3) + 
    xlab("") +
    theme_BM() + 
    theme( strip.background = element_blank(),
                        strip.text.x = element_blank(),
           axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5, hjust=1),
           plot.margin = margin(20, 8, 8, 8, unit = "pt"))
  F3a_seqlogo_plots
  
  #################################################
  # Final test: perform pairwise enrichment of nucleotides in each category 
  #################################################
  # 5 different cat (-4A, -3A, -2 A, -1A, +1T)
  contexts = rbindlist(TN_contexts , idcol = "exposure")
  contexts_pattern = strsplit(contexts$context,"") %>% unlist() %>% matrix(ncol = 21, byrow = T)
  contexts = cbind(contexts_pattern, contexts)
  
  # generate all unique combinations of selected peaks from the EcC motif
  enrichments = c("7A", "8A", "9A", "10A", "12T")
  nucs_select_2 = combinations(n = 5,r = 2,  enrichments, repeats.allowed = F)  %>% t() %>% as.data.frame() %>%  as.list()
  nucs_select_3 = combinations(n = 5, r = 3,v =  enrichments, repeats.allowed = F) %>% t() %>% as.data.frame() %>%  as.list()
  nucs_select_4 = combinations(n = 5, r = 4,v =  enrichments, repeats.allowed = F) %>% t() %>% as.data.frame() %>%  as.list()
  
  list_total_enrichments = c(enrichments, nucs_select_2, nucs_select_3, nucs_select_4, list(enrichments))
  recode_names = lapply(list_total_enrichments, function(x) dplyr::recode(x, `7A` = "-4A", `8A` = "-3A", `9A` = "-2A", `10A` = "-1A", `12T` = "+1T"))
  names(list_total_enrichments) = sapply(recode_names, paste0, collapse = " ")
  
  # combinations of two
  test_contexts = contexts %>% as.data.frame() %>% filter(exposure %in% c("EcC", "Control"))
  colnames(test_contexts)[1:21] = 1:21
  
  test_list = list()
  
  set.seed(12356)
  test_table = data.frame(EcC = rep(NA, length(list_total_enrichments)))
  rownames(test_table) = names(list_total_enrichments)
  
  sim_contexts_pattern = strsplit(sim_contexts$context,"") %>% unlist() %>% matrix(ncol = 21, byrow = T)
  colnames(sim_contexts_pattern) = 1:21
  sim_contexts_bases = cbind(sim_contexts_pattern, sim_contexts)  
  sim_ctrl_contexts = sim_contexts_bases %>% dplyr::filter(name == "Control")
  
  pvalue_list = list()
  for (i in 1:length(list_total_enrichments)) {
    nucs = list_total_enrichments[[i]]
    name = names(list_total_enrichments)[[i]]
    pos = parse_number(nucs) %>% as.character()
    base = gsub(".*[0-9]", "", nucs)
    
    # 1. check for the observed data
    base_check = test_contexts[,pos] == base 
    
    if (length(nucs) > 1 ) { 
      base_check %>% as.matrix()
      idx  = apply(base_check, 1,all)
    } else { idx = base_check}
    
    motif_match = test_contexts$exposure[idx] %>% table()
    motif_nomatch = test_contexts$exposure[!idx] %>% table()
    mat = rbind(motif_nomatch, motif_match) %>% as.matrix()
    
    # fisher exact test for EcC vs control
    test_table[i,1] = fisher.test(mat[,c("Control", "EcC")], alternative = "greater")$p.value
    
    # 2. check for the simulated data
    base_check = sim_ctrl_contexts[,pos] == base 
    if (length(nucs) > 1 ) { 
      base_check %>% as.matrix()
      idx  = apply(base_check, 1,all)
    } else { idx = base_check}
    
    motif_match = sim_ctrl_contexts[idx, c("name","bin")] %>% group_by(bin) %>% dplyr::count()
    motif_nomatch = sim_ctrl_contexts[!idx, c("name","bin")] %>% group_by(bin) %>% dplyr::count()
    motif = rbindlist(list(motif_match = motif_match, motif_nomatch = motif_nomatch), idcol = "motif") %>% 
      pivot_wider(names_from = motif, values_from = n, values_fill = 0)
    
    pvals = vector("numeric", 200)
    # perform the fisher test for the control data:
    
    motif$control_match = mat[2,1]
    motif$control_nomatch = mat[1,1]
  
    # perform test for the total pvalue 
    pvals = apply(motif, 1, \(x) fisher.test(matrix(x[-1], nrow =2), alternative = "greater")$p.value)
    pvalue_list[[name]] = pvals 
  }
  
  test_tibble = test_table %>% 
    rownames_to_column("context") %>% 
    arrange(EcC)
  observed_data = test_tibble %>% 
    pivot_longer(EcC, names_to = "type", values_to = "pvalue")
  
  colnames(observed_data)
  simulation_data = as.data.frame(pvalue_list) %>%
    pivot_longer(cols = everything(), values_to = "pvalue", names_to = "context") %>%
    mutate(context = rep(names(list_total_enrichments), 200),
           type = "simulation")
  
  pval_contexts = rbind(observed_data, simulation_data) %>% 
    mutate(context = factor(context, levels = test_tibble$context), 
           `-log10 p-value` = -log10(pvalue))
  
  position_enrichment_plot = ggplot(pval_contexts, aes(x = context, y = `-log10 p-value`, color = type)) + geom_point() + 
    theme_classic() + 
    scale_y_continuous(expand = c(0, 3), limits = c(0, NA )) + 
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.title.y = element_text(size = 11), 
          plot.margin = margin(unit(c(5.5, 5.5, 0, 20), "points")), 
          legend.position = c(0.7,0.7), legend.box.background = element_rect(color = "black")) + 
    scale_color_manual(values = c("black", "grey")) + 
    xlab(NULL) + labs(color = NULL)
  
  legend_mat = matrix('N', nrow = 6, ncol = 31) 
  colnames(legend_mat) = list_total_enrichments
  
  for (i in 1:ncol(legend_mat)) { 
    name = names(list_total_enrichments)[i]
    index = gsub("[[:alpha:]]|\\+", "", name) %>% strsplit(split = " ") 
    index = as.numeric(index[[1]]) + 5
    # index = index*-1 + 7
    base = gsub("[^[:alpha:]]", "",name) %>% strsplit("")
    legend_mat[index,i] = base[[1]]
  }
  
  colnames(legend_mat) = names(list_total_enrichments)
  legend_mat[5,] = "T" # set base mutation to T
  
  # reorder matrix
  legend_mat_m =  melt(legend_mat, varnames = c("base position", "names"), value.name = "base")
  test_table = test_table %>% rownames_to_column("names")
  position_enrichment = merge(test_table, legend_mat_m, by = "names")

  bases = ggplot(position_enrichment , aes(x = reorder(names, EcC), y = `base position`, fill = base)) + 
    geom_tile(color = "white", linewidth = 0.9 ) +
    geom_text(aes(label = base ), color = "grey15", size = 2) + theme_minimal() + 
    scale_fill_manual(values = c("green", "grey", 'red')) + 
    theme(axis.title.y = element_text(size = 11),
          axis.text.y = element_text(size=7),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          plot.margin = margin(unit(c(-5, 5.5, 5.5, 20), "points")),
          legend.position = "none") + 
    ylab("base position") +
    scale_y_continuous(breaks = 6:0, labels = c(1:-5)) + 
    labs(fill = element_blank())
  
   F3b_position_enrichment = position_enrichment_plot / plot_spacer() / bases + plot_layout(heights = c(0.75, -0.17, 1))
  F3b_position_enrichment
  
  ###### Perform fisher exact test for selected trinucleotide combinations 
  # See if further selection of AA motif at selected trinucleotides works better
  # do not iterate over all combinations - start wit the most frequently occurring basepair
  EcC_signature = signatures$SBS88
  EcC_signature_ordered = signatures$Type[order(EcC_signature, decreasing = T)]
  EcC_signature_ordered = EcC_signature_ordered[substr(EcC_signature_ordered, 3,3) == "T"] # select only T>N mutations 
  
  # combinations of two
  temp_table = contexts %>% filter(exposure %in% c("Control", "EcC")) %>% as.data.frame()
  colnames(temp_table)[1:21] = 1:21
  test_table = data.frame(names = EcC_signature_ordered,  pval = rep(NA, length(EcC_signature_ordered)))
  
  sim_contexts_ctrl_EcC = sim_contexts %>% filter(name %in% c("Control", "EcC"))
  
  sim_data = list()
  for (j in 1:length(EcC_signature_ordered)) {
  
    trinucs = EcC_signature_ordered[1:j]
    trinucs_index = temp_table$trinucleotide %in% trinucs
    AA_index = substr(temp_table$context, 7,8) == "AA"
    idx = trinucs_index & AA_index
    motif_match = temp_table$exposure[idx] %>% table()
    motif_nomatch = temp_table$exposure[!idx] %>% table()
    mat = rbind(motif_match, motif_nomatch)
    
    # fisher exact test for EcC vs control
    test_table[j,2] = fisher.test(mat[,2:1], alternative = "greater",)$p.value
    
    # perform the same analysis for the simulated data
    Ctrl_sim_muts = sim_contexts_ctrl_EcC %>% filter(name == "Control")
    trinucs_index = Ctrl_sim_muts$trinucleotide %in% trinucs
    AA_index = substr(Ctrl_sim_muts$context, 7,8) == "AA"
    Ctrl_sim_muts$idx = trinucs_index & AA_index
    Ctrl_sim_muts = Ctrl_sim_muts %>% 
      mutate(idx = factor(idx, levels = c(TRUE, FALSE)))
    
    motif_match = Ctrl_sim_muts[  ,c("idx", "bin")] %>% table() %>% 
      as.matrix() %>% t() %>%
      as.data.frame.matrix() %>% as_tibble()
    colnames(motif_match) = c("sim_match", "sim_nomatch")
    motif_match = motif_match %>% 
      mutate(obs_match = mat[1,1], obs_nomatch = mat[2,1])
    
    matrix(motif_match[1,], ncol =2, byrow = TRUE)
    
    pvals = apply(motif_match, 1, \(x) fisher.test(matrix(x, ncol = 2, byrow =TRUE), alternative = "greater")$p.value)  
    nucleotide = test_table[j,"names"]
    sim_data[[nucleotide]] = data.table(pval = pvals)
  }

  
  sim_table = rbindlist(sim_data, idcol = "names") %>% 
    mutate(type = "simulation")
  test_table = test_table %>% 
    mutate(names = factor(names, levels = EcC_signature_ordered),
           type = "observed")
  total_table = rbind(sim_table, test_table)
  
  # save the 17 trinucleotides in the data: 
  write.table(total_table, "trinucs_selection.tsv", sep = "\t", row.names = F)

  F3f_stepwise_trinucs = ggplot(total_table, aes(x = names, y = -log10(pval), color = type)) +
    geom_point() +
    theme_BM() + 
    scale_color_manual(values = c("black", "grey")) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,  size = 7), 
          legend.position = c(0.8,0.4), legend.box.background = element_rect(color = 'black'),
          plot.margin = margin(unit(c(5.5, 5.5, 5.5, 5.5), "points"))) + 
    xlab("Trinucleotide added")  + ylab("-log10 p-value\nenrichment -3-4AA") + 
    labs(color = NULL)
  
  # conclusion: 17 most frequently occuring trinucleotides is the best value. 
  # fisher test for these values for EcN 
  trinucs = EcC_signature_ordered[1:17] # select the 19 most commonly mutated samples 
  contexts17 = ext_context %>% mutate(motif17 = ifelse(trinucleotide %in% trinucs & pos34 == "AA", "motif17", "nomotif"))
  mat = table(contexts17$name, contexts17$motif17) %>% t()
  pval_nissle = fisher.test(mat[,c("EcN","Control")], alternative = "greater")$p.value 
  print(paste("pval_nissle 17:", format(pval_nissle)))
  pval_EcC = fisher.test(mat[,c("EcC","Control")], alternative = "greater")$p.value # p-values = 2.670921e-19
  print(paste("pval_nissle 17:", format(pval_EcC)))
  # step 1: Generate randomly sampled combinations EcC and control data
  # sample 10 times for each percentage point. Sample the same number of mutations as EcN = 983 in total
  # step 2: perform fisher exact testtest on presence %>% of pos34 AA presence 
  # step 3: perform signature extraction 
  
  set.seed(123546)
  
  list_fractions = list()
  list_EcC_sigs = list()
  
  EcC_count = sum(ext_context$name == "EcC")
  for (fraction_EcC in seq(0.0,1, 0.01)) {
    
    n_EcC = round(EcC_count*fraction_EcC)
    n_control = round(EcC_count*(1-fraction_EcC))
    
    data_fraction = data.frame(replicate = 1:1, p_value = rep(NA, 1), odds_ratio = NA, lower_conf = NA, higher_conf = NA)
    set.seed(fraction_EcC * 1000)
  # EcC mutation sampling and counting 
    ctrl_muts = dplyr::slice_sample(TN_contexts$Control, n = n_control, replace = T)
    
    if (fraction_EcC > 0 ) {
      EcC_muts = dplyr::slice_sample(TN_contexts$EcC, n = n_EcC, replace = T)
      total_muts = rbind(EcC_muts, ctrl_muts) 
    } else {total_muts = ctrl_muts}
    
    test_muts_AA = sum(substr(total_muts$context, 7,8) == "AA" & total_muts$trinucleotide %in% trinucs[1:17])
    test_muts_noAA = nrow(total_muts) - test_muts_AA
    
    
    ctrl_muts_AA = sum(substr(TN_contexts$Control$context, 7,8) == "AA"  & TN_contexts$Control$trinucleotide %in% trinucs[1:17])
    ctrl_muts_noAA = nrow(TN_contexts$Control)- ctrl_muts_AA
    
    # fisher test 
    fmat = matrix(c(test_muts_AA, test_muts_noAA, ctrl_muts_AA, ctrl_muts_noAA), ncol = 2)
    
    ftest = fisher.test(fmat)
    data_fraction[i,] = c(i, ftest$p.value, ftest$estimate, ftest$conf.int[1], ftest$conf.int[2])
    list_fractions[[as.character(fraction_EcC)]] = data_fraction
    }
  
  
  total_fractions = rbindlist(list_fractions, idcol = "fraction")
  total_fractions$type = "sampled Control and EcC mixtures 0 - 100%"
  
  # calculate  odds ratio for EcN
  
  ftests = lapply(colnames(mat)[-1], \(x) fisher.test(mat[,c(x, "Control" )]))
  ftests = lapply(ftests, broom::tidy) %>% rbindlist() %>% as.data.frame()
  rownames(ftests) = colnames(mat)[-1]
  colnames(ftests)[1:4] = c("odds_ratio", "p_value", "lower_conf", "higher_conf")
  ftests = ftests %>% 
    dplyr::select(odds_ratio, p_value, lower_conf, higher_conf)
  ftests$replicate = 1
  ftests$fraction = 1
  ftests$type = rownames(ftests)
  fr_total = rbind(total_fractions, ftests)
  fr_total$type = factor(fr_total$type, levels = unique(fr_total$type))
  
  F3g_mixture_plot_EcN = ggplot(fr_total, aes(x = fraction, y = odds_ratio, ymin = lower_conf, ymax = higher_conf )) +
    geom_pointrange(alpha = 0.6, size = 0.3, position = position_dodge(width = 0.1)) + facet_grid( . ~ type, scales = "free_x", space = "free") + 
    ggh4x::force_panelsizes(cols = c(10, 1,1,1,1)) + 
    theme_half_open() + panel_border(size = 0.5) + geom_hline(yintercept = 1) + 
    theme(axis.text.x = element_blank(), 
          panel.spacing.x = unit(0.1, "lines"),
          strip.text.x = element_text(size = 9),
          axis.title.x = element_text(size = 12), 
          axis.title.y = element_text(size = 11),
          plot.margin = margin(unit(c(5.5, 5.5, 5.5, 5.5), "points"))) + 
    xlab("fraction Control/EcC mutations (1% difference/step)                 ") + 
    ylab("-3-4AA enrichment\nodds ratio(Fisher test)") +
    theme()
  
  total_fractions$fraction = as.numeric(total_fractions$fraction)
  model = lm(fraction ~ odds_ratio , data = total_fractions) 
  
  # estimations for relative mutagenicity of the E.coli strains
  ftests$estimate_mean = predict(model, newdata = ftests["odds_ratio"])
  ftests$estimate_low = predict(model, newdata = data.frame(odds_ratio = ftests$lower_conf))
  ftests$estimate_high = predict(model, newdata = data.frame(odds_ratio = ftests$higher_conf))
  
  print(ftests)

  left = ggarrange(F3a_seqlogo_plots, ggarrange(F3c_dinuc_frequencies,  F3d_dinc_enrichment, widths = c(1.5 ,1), labels = c("C", "D")), 
                   nrow = 2,  heights = c(1.8, 1), labels = c("A", ""))
  right = ggarrange(F3b_position_enrichment, F3e_AA_context_profile,nrow = 2, heights = c(1,2.2), labels = c("B", "E"))
  top = ggarrange(left, right, widths = c(1.3,1))
  total_plot = ggarrange(top, ggarrange(F3f_stepwise_trinucs, F3g_mixture_plot_EcN, labels = c("F", "G")), nrow = 2, heights = c(2.5,1))
  
  
  # supplementary plot 
  supp_plot = ggarrange(simulation_plot + theme(legend.position = "right"), histogram_fisher, nrow = 2, heights = c(1.5, 1))
  
  plot_list = list(seqlogo = F3a_seqlogo_plots, dinuc_freqs = F3c_dinuc_frequencies, 
                   dinuc_enrichment = F3d_dinc_enrichment, motif_enrichment = F3b_position_enrichment, 
                   stepwise_trinucs = F3f_stepwise_trinucs, mixture_plot = F3g_mixture_plot_EcN, AA_profile = F3e_AA_context_profile, 
                   simulation_plot = simulation_plot, histogram_fisher = histogram_fisher, supplementary_figure_4 = supplementary_figure_4)
  
  return(list(total_plot = total_plot, supp_plot = supp_plot,  plot_list = plot_list))
}

# plot figure 2
plot_PTA = plot_figures_2(contexts_TN_PTA, cat = cat_PTA, name = 'PTA')
ggsave("Output/Figures/Figure_2.pdf", plot_PTA$total_plot, width = 15, height = 11.5)
ggsave("Output/Figures/Figure_2.png", plot_PTA$total_plot, width = 15, height = 11.5)

# plot figure S3
pCE = plot_figures_2(contexts_TN_CE, cat = cat_CE, name = 'Clonal Expansion')
ggsave("Output/Figures/Fig_S3.pdf", pCE$total_plot, width = 14, height = 11.5)
ggsave("Output/Figures/Fig_S3.png", pCE$total_plot, width = 14, height = 11.5)

supp_figure_4 = ggarrange(plot_PTA$supp_plot, 
          ggarrange(pCE$supp_plot, plot_spacer()  + bgcolor("white")), 
          nrow = 4, labels = c("A", "B", "C", "D"))
ggsave("Output/Figures/Fig_S4.pdf", supp_figure_4, width = 8, height = 11)
ggsave("Output/Figures/Fig_S4.png", supp_figure_4, width = 8, height = 11)

# get the output for the new fused supplementary Figure: 
supp_fig_list = list(plot_PTA$plot_list$simulation_plot, pCE$plot_list$simulation_plot,
                     plot_PTA$plot_list$histogram_fisher, pCE$plot_list$histogram_fisher)
supp_fig_list = lapply(supp_fig_list, \(x) x + theme(plot.margin = unit(c(12, 4, 12, 4), "mm")))

bottom_supp_figure = cowplot::plot_grid(supp_fig_list[[1]], supp_fig_list[[2]],
                                        supp_fig_list[[3]], supp_fig_list[[4]], labels = c("H", "J", "I", "K"), rel_widths = c(1.65,1))

new_combined_fig = ggarrange(pCE$total_plot, bottom_supp_figure, nrow = 2, heights = c(1.5,1))
ggsave("Output/Figures/Fig_S2_new.pdf", new_combined_fig, width = 13, height = 19)
ggsave("Output/Figures/Fig_S2_new.png", new_combined_fig, width = 13, height = 19)