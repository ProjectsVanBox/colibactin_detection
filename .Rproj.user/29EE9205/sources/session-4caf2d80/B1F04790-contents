# New figure 2 plotting functions, to easily performa the analysis for both the PTA and the non-PTA samples
# Nissle Figure 2
# Axel Rosendahl Huber 13-12-2023
setwd("C:/Users/Axel Rosendahl Huber/OneDrive/Nissle_manuscript/Nissle/")
library(gtools)
source("Code//Load_data.R")
source("Code/Functions/Utils.R")

cat_PTA = categories %>%
  filter(method == "PTA")

cat_CE = categories %>% 
  filter(method == "Clonal Expansion")
cat_CE$injection = factor(cat_CE$injection, levels = c("Control", "EcN", "EcC"))

plot_figure_2 = function(cat, mut_mat, name){

  # Boxplots indicating total mutation counts 
  mut_counts = cat  %>% 
    mutate('SBS count' = colSums(mut_mat[, name]))
  
  ktest = kruskal.test(`SBS count` ~ injection, data = mut_counts)
  subtitle = paste0("Kruskal-Wallis p-value: ", round(ktest$p.value, 3))
  
  F2a_sbs_boxplot = ggplot(mut_counts, aes(x = injection, y = `SBS count`, fill = injection)) + 
    geom_boxplot(outlier.shape = NA, width = 0.4) + 
    geom_jitter(shape = 21, width = 0.15) +  
    geom_pwc(aes(group = injection), dodge = 0.2, method = "dunn_test",
             p.adjust.by = "panel", label = "{p.adj.format}", p.adjust.method = "fdr") + 
    scale_fill_manual(values = c("#545863", "#00e8fc", "#f96e46", "#f9c846", "#ffe3e3")) +
    scale_y_continuous(limits = c(0, max(mut_counts$`SBS count`)*1.45)) +
    theme_BM() + 
    theme(plot.title = element_text(hjust = 0.5, size = 11), legend.position = "none", 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    labs(title = name, subtitle = subtitle, x = "")
  
  mut_mat = mut_mat[,cat$name]
  
  # Plot 96-profiles 
  # Make mutation matrices & plot mutation profiles for I3 and PTA-sequenced clones
  mm = t(mut_mat) %>% data.table(keep.rownames = "name")
  mm = merge(mm, mut_counts[,-4])
  mm_cat = mm %>% group_by(injection, method) %>%
    dplyr::select(-name) %>%
    summarize_all(sum) %>%
    dplyr::select(-method) %>% 
    column_to_rownames("injection") %>% 
    t()
  
  F2c_sbs_profile = plot_96_profile3(mm_cat) + 
    theme(axis.text.x = element_text(size = 4.5))
  
  # I3 DYE and I3 Nissle indel analysis
  # Boxplots indicating total mutation counts 
  indel_mut_counts = mut_counts %>% 
    mutate('indel count' = colSums(indel_counts[, name]))
  
  # Save a table with the information per clone
  indel_mut_counts = indel_mut_counts %>% arrange(method, injection)
  sjPlot::tab_df(indel_mut_counts, file = "Output/Supplementary_tables/Supplementary_Table_1.doc")
  
  #  Plot indel counts and mutational profiles
  indel_mut_counts_plot = indel_mut_counts
  
  ktest_indels = kruskal.test(`indel count` ~ injection, indel_mut_counts_plot)
  subtitle = paste0("Kruskal-Wallis p-value: ", round(ktest_indels$p.value, 3))
  
  F2b_boxplot_indels = ggplot(indel_mut_counts_plot, aes(x = injection, y = `indel count`, fill = injection)) + 
    geom_boxplot(outlier.shape = NA, width = 0.4) + 
    geom_jitter(shape = 21, width = 0.15) +
    geom_pwc(aes(group = injection), method = "dunn_test",
             dodge = 0.2, p.adjust.by = "panel", label = "{p.adj.format}", p.adjust.method = "fdr") + 
    scale_fill_manual(values = c("#545863", "#00e8fc", "#f96e46", "#f9c846", "#ffe3e3")) +
    scale_y_continuous(limits = c(0, max(indel_mut_counts$`indel count`)*1.45)) +
    theme_BM() +
    theme(plot.title = element_text(hjust = 0.5, size = 11),
          legend.position = "none",
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    labs(title = "Indels", subtitle = subtitle, x = "")
  
  # Plot indel profiles 
  indel_counts = indel_counts[, cat$name]
  id_mm = t(indel_counts) %>% data.table(keep.rownames = "name")
  id_mm = merge(id_mm, indel_mut_counts[,-4])
  id_mm = id_mm %>% group_by(injection, method) %>%
    dplyr::select(-name, -`indel count`) %>%
    summarize_all(sum) %>% 
    dplyr::select(-method) %>% 
    column_to_rownames("injection") %>% 
    t()
  
  F2d_indel_profile = plot_indel_contexts2(id_mm) +
    theme_minimal_hgrid() +
    scale_y_continuous(n.breaks = 3) + 
    theme(panel.spacing.x = unit(0, "mm"),
          legend.position = "none",
          axis.text.x = element_text(size = 5), 
          strip.text.y = element_text(size = 9), 
          axis.text.y = element_text(size = 8))
  
  # get cosine similarity EcN to Control  
  cos_sim(id_mm[,"Control"], id_mm[, "EcN"])
  
  #  ===== SBS signature re-fitting for each clone ======
  sigs_organoid_culture = signatures[, c("SBS1", "SBS5", "SBS18", "SBS88")]
  fit_res = fit_to_signatures(mut_mat, as.matrix(sigs_organoid_culture))
  
  fit_res_clones_sbs = fit_res$contribution %>% 
    prop.table(2) %>% 
    as.data.frame() %>% rownames_to_column("Signature") %>% 
    pivot_longer(cols = -Signature) %>% 
    filter(Signature == "SBS88")
  fit_res_clones_sbs = merge(fit_res_clones_sbs, categories)
  
  
  ktest_refit = kruskal.test(value ~ injection, data = fit_res_clones_sbs)
  subtitle = paste0("Kruskal-Wallis p-value: ", round(ktest_refit$p.value, 3))
  
  
  F2g_sbs_refit = ggplot(fit_res_clones_sbs, aes(x = injection, y = value, fill = injection,shape = method))  +
    geom_boxplot(aes(shape = NULL), outlier.shape = NA, width = 0.5) + 
    geom_jitter(shape = 21, width = 0.15, alpha = 0.8)  + 
    theme_BM() +
    geom_pwc(aes(group = injection), method = "dunn_test", p.adjust.method =  "fdr", label = "{p.adj.format}") +
    scale_fill_manual(values = c("#545863", "#00e8fc", "#f96e46", "#f9c846", "#ffe3e3")) +
    scale_y_continuous(limits = c(0, max(fit_res_clones_sbs$value)*1.45)) +
    labs( y = "Relative SBS88 contribution", x ="", subtitle = subtitle) +
    theme(legend.position =  "none", 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  F2g_sbs_refit
  
  # ===== Indel signature refitting ====== 
  id_sigs_select = id_signatures[, c("ID1", "ID2", "ID18")]
  fit_res_id = fit_to_signatures(indel_counts, as.matrix(id_sigs_select))
  fit_res_clones = fit_res_id$contribution %>% 
    prop.table(2) %>% 
    as.data.frame() %>% rownames_to_column("Signature") %>% 
    pivot_longer(cols = -Signature) %>% 
    filter(Signature == "ID18")
  fit_res_clones = merge(fit_res_clones, categories)

  
  ktest_refit = kruskal.test(value ~ injection, data = fit_res_clones)
  subtitle = paste0("Kruskal-Wallis p-value: ", round(ktest_refit$p.value, 3))
  
  Fig2h_indel_refit = ggplot(fit_res_clones, aes(x = injection, y = value, fill = injection,shape = method))  +
    geom_boxplot(aes(shape = NULL), outlier.size = 0, width = 0.5) + 
    geom_jitter(shape = 21, width = 0.15, alpha = 0.8)  + 
    theme_BM() +
    geom_pwc(aes(group = injection), method = "dunn_test", p.adjust.method =  "fdr", label = "p.adj.format") +
    scale_fill_manual(values = c("#545863", "#00e8fc", "#f96e46", "#f9c846", "#ffe3e3")) +
    scale_y_continuous(limits = c(0, max(fit_res_clones$value)*1.30)) +
    labs(subtitle = subtitle, y = "Relative ID18 contribution", x = "") +
    theme(legend.position =  "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
  
  # patch all plots together
  figure_2 = F2a_sbs_boxplot + F2b_boxplot_indels + F2c_sbs_profile + F2d_indel_profile + 
    F2g_sbs_refit + Fig2h_indel_refit + 
    plot_layout(guides = "collect", byrow = F, nrow = 2,widths = c(1,1.5,1)) + plot_annotation(tag_levels = "A")
  
  # Rebuttal figure 2: Mutational loads for signatures in absolute counts
  fit_res_clones = fit_res_id$contribution %>% 
    as.data.frame() %>% rownames_to_column("Signature") %>% 
    pivot_longer(cols = -Signature) %>% 
    filter(Signature == "ID18")
  fit_res_clones = merge(fit_res_clones, categories)
  
  rebuttal_2_indel_refit = ggplot(fit_res_clones, aes(x = injection, y = value, fill = injection,shape = method))  +
    geom_boxplot(aes(shape = NULL), outlier.size = 0, width = 0.15, alpha = 0.6) + 
    geom_jitter(width = 0.15, alpha = 0.8)  + 
    theme_BM() +
    geom_pwc(aes(group = injection), ref.group	= "Control", p.adjust.method =  "fdr", label = "p.adj") +
    scale_fill_manual(values = c("#545863", "#00e8fc", "#f96e46", "#f9c846", "#ffe3e3")) +
    ylab("Absolute ID18 contribution") + xlab("") +
    theme(legend.position =  "none", 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    ggtitle("Indels")
  
  fit_res_clones_sbs = fit_res$contribution %>% 
    # prop.table(2) %>% 
    as.data.frame() %>% rownames_to_column("Signature") %>% 
    pivot_longer(cols = -Signature) %>% 
    filter(Signature == "SBS88")
  fit_res_clones_sbs = merge(fit_res_clones_sbs, categories)
  
  rebuttal_2_sbs_refit = ggplot(fit_res_clones_sbs, aes(x = injection, y = value, fill = injection,shape = method))  +
    geom_boxplot(aes(shape = NULL), outlier.shape = NA, width = 0.15, alpha = 0.6) + 
    geom_jitter(width = 0.15, alpha = 0.8)  + 
    theme_BM() +
    geom_pwc(aes(group = injection), ref.group	= "Control", p.adjust.method =  "fdr", label = "p.adj") +
    scale_fill_manual(values = c("#545863", "#00e8fc", "#f96e46", "#f9c846", "#ffe3e3")) +
    ylab("Absolute SBS88 contribution") + xlab("") +
    theme(legend.position =  "none", 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    ggtitle("singe base substitutions") 

  # alternative figure 2
  left = ggarrange(F2a_sbs_boxplot, F2b_boxplot_indels, labels = c("B", "C"), nrow = 2)
  middle = ggarrange(F2c_sbs_profile, F2d_indel_profile, labels = c("D", "E"),  nrow = 2)
  right = ggarrange(F2g_sbs_refit, Fig2h_indel_refit, labels = c("F", "G"), nrow = 2)
  figure_1 =  ggarrange(left, middle, right, widths = c(1,  2, 1), ncol = 3)

  plot_list = list(F2a_sbs_boxplot = F2a_sbs_boxplot, 
                   F2b_boxplot_indels = F2b_boxplot_indels, 
                   F2c_sbs_profile = F2c_sbs_profile, 
                   F2d_indel_profile = F2d_indel_profile, 
                   F2g_sbs_refit = F2g_sbs_refit, 
                   Fig2h_indel_refit = Fig2h_indel_refit)
  
  return(list(plot_list = plot_list, fig = figure_1))
}

fig_1_PTA = plot_figure_2(cat = cat_PTA, mut_mat, name = "PTA")
ggsave("Output/Figures/Figure_1_R.pdf", fig_1_PTA$fig, width = 13, height = 8)
ggsave("Output/Figures/Figure_1_R.png", fig_1_PTA$fig, width = 13, height = 8)

# analyses for all mutations:
mm = t(mut_mat) %>% data.table(keep.rownames = "name")
mm = merge(categories, mm, by = "name")
mm_cat = mm %>% group_by(injection, method) %>%
  dplyr::select(-name) %>%
  summarize_all(sum) %>%
  mutate(name = paste0(method, "-", injection)) %>% 
  mutate(name = gsub("Clonal Expansion-", "CE", name)) %>% 
  ungroup() %>% 
  dplyr::select(-method, -injection) %>% 
  column_to_rownames("name") %>% 
  as.matrix() %>% 
  t()

cossim_all = cos_sim_matrix(mm_cat, mm_cat)
cos_mat_sbs = plot_cosine_heatmap2(cossim_all, cluster_cols = TRUE, plot_values = TRUE)

# Indel cosine similarities: 
id_mm = t(indel_counts) %>% data.table(keep.rownames = "name")
id_mm = merge(id_mm, categories) %>% group_by(injection, method) %>%
  dplyr::select(-name) %>%
  summarize_all(sum) %>%
  mutate(name = paste0(method, "-", injection)) %>% 
  mutate(name = gsub("Clonal Expansion-", "CE", name)) %>% 
  ungroup() %>% 
  dplyr::select(-method, -injection) %>% 
  column_to_rownames("name") %>% 
  as.matrix() %>% 
  t()

cosine_indel_all = cos_sim_matrix(id_mm, id_mm)
cos_mat_indel = plot_cosine_heatmap2(cosine_indel_all, cluster_cols = TRUE, plot_values = TRUE) + 
  guides(fill = "")
CE = plot_figure_2(cat = cat_CE, mut_mat, name = "CE")
CE = CE$plot_list


# plot the different samples:
left_top = ggarrange(CE$F2a_sbs_boxplot, CE$F2b_boxplot_indels, labels = c("B", "C"), nrow = 2)
right_top = ggarrange(CE$F2c_sbs_profile, CE$F2d_indel_profile, labels = c("D", "E"),  nrow = 2)
top = ggarrange(left_top, right_top, ncol = 2, widths = c(1,2))

left_bottom = ggarrange(CE$F2g_sbs_refit, CE$Fig2h_indel_refit, labels = c("F", "G"), nrow = 2)
right_bottom = ggarrange(cos_mat_sbs, cos_mat_indel, labels = c('H', "I"), nrow = 2)
bottom = ggarrange(left_bottom, right_bottom, ncol = 2, widths = c(1,2))
fig_S2 =  ggarrange(top, bottom,  ncol = 2)

ggsave("Output/Figures/Fig_S2_R.pdf", 
       fig_S2, width = 16, height = 9)
ggsave("Output/Figures/Fig_S2_R.png", 
       fig_S2, width = 16, height = 9)