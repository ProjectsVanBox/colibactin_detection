# Figure 5v2 
# POLE mutation status 

# A: HMF WGS pval vs AA (color-coded by mutation load)
# B: HMF WGS highlight hypermutators by color and point out with a line the ones with POLE mutations (as it is in Axel plot right now)
# C: HMF WGS T>N plots with the light color for all and dark color for AA enrichment
# D: HMF AA enrichment in the top 4 trinucleotides of POLE mut sig
# E: HMF AA enrichment in the rest

library(ggrepel)
library(ggExtra)
library(viridis)
library(ggseqlogo)
source("C:/Users/Axel Rosendahl Huber/Documents/Shared/vanBoxtelLab (Groupfolder)/Projects/Axel/Nissle/Manuscript/Figures/Figure_4/Figure_4_v3.R")
setwd("C:/Users/Axel Rosendahl Huber/Documents/Shared/vanBoxtelLab (Groupfolder)/Projects/Axel/Nissle/Analysis/")
source("C:/Users/Axel Rosendahl Huber/Documents/Shared/vanBoxtelLab (Groupfolder)/Projects/Axel/Nissle/Analysis/Scripts/Nissle_functions.R")

# get samples with somatic POLE mutations 
POLE_muts = list.files("New_HMF/POLE_hotspot_muts/", full.names = T)
names(POLE_muts) = gsub("_.*", "", basename(POLE_muts))
POLE_changes = sapply(POLE_muts, function(x) strsplit(read.table(x)[,8], "\\|")[[1]][11])
metadata$POLE_somatic = "WT"
metadata$POLE_somatic[match(names(POLE_changes), metadata$sampleId)] = POLE_changes
metadata = metadata %>% mutate(POLE_somatic = replace(POLE_somatic,!grepl("Val411Leu|Pro286Arg|Ser459Phe", POLE_somatic), "WT"))
metadata$log10muts = log10(metadata$n_sbs)

# Figure 5A - mutation numbers HMF cohort by color and POLE mutations
F5A = ggplot(metadata, aes(x = AA, y = log_p_wgs, color = log10muts)) +
  geom_point(alpha = 0.8, size = 1) + 
  geom_hline(yintercept = 3) + geom_vline(xintercept = 0.22) + 
  scale_color_viridis() + 
  geom_text_repel(data = metadata %>% filter(POLE_somatic != "WT"),
                  mapping =  aes(label = POLE_somatic), min.segment.length = .05, color = "black",size = 3) + 
  scale_shape_discrete(name = "Pks classification", labels = c("pks negative", "pks+ established", "pks+ new")) + 
  xlab("fraction mutations with AA at\n-3-4 position at pks-sites") + ylab("-log10 p-value") + theme_classic() + 
  theme(legend.position=c(.8,0.4),  legend.box.background = element_rect(colour = "black"))
F5A

ggsave("../Manuscript/Figures/Figure_5/Figure_5A.pdf", F5A, width = 4.5, height = 3.5)

# figure 5C: 
# select samples with true pole motifs 
POLE_samples = metadata %>% filter(POLE_somatic != "WT") %>% pull(sampleId)
F5c_96_trinucleotide_matrix = plot_96_profile3(mm_HMF[, POLE_samples])
TRIPLETS_48 = TRIPLETS_96[49:96]
SBS88_TN <- as.data.table(signatures) %>% dplyr::slice(49:96) %>% pull("SBS88")
SBS28_TN <- as.data.table(signatures) %>% dplyr::slice(49:96) %>% pull("SBS28")


# get mutation table for exome data # hpc needs to be connected for this plot
POLE_list = vector(mode = "list")
for ( sample in POLE_samples) {
  print(sample)
  file = paste0("C:/Users/Axel Rosendahl Huber/Documents/Shared/vanBoxtelLab (Groupfolder)/Projects/Axel/Nissle/Analysis/Processed_data/HMF_vcfs/",
                sample, "_10basepair_context_SNVs.txt.gz")
  POLE_list[[sample]]  = data.table::fread(file)
}

POLE_muts = rbindlist(POLE_list, idcol = "name")
POLE_muts = POLE_muts %>% filter(triplet %in% TRIPLETS_96[49:96]) %>% 
  mutate(pos34 =  substr(context, 7,8)) %>% 
  mutate(trinucleotide = factor(triplet, levels = TRIPLETS_96[49:96])) %>% 
  mutate(select =  ifelse(pos34 == "AA", "AA", "other") %>% 
           factor(levels = c("other", "AA")))






label_df_SBS28 = get_profile_labels(POLE_muts, SBS28_TN)
label_df = get_profile_labels(POLE_muts, SBS88_TN)
label_df[, "cosine"] = paste0("cosine similarity to SBS88:\n",label_df$cosine,
                              "\ncosine similarity to SBS28: ", label_df_SBS28$cosine)

plot_profile_absolute = function(mut_list) {
  
  mut_list %>% group_by(trinucleotide, select) %>% 
    ggplot(aes(x = trinucleotide, alpha = select, fill = type)) +
    geom_bar(width = 0.7 ) + facet_grid(name ~ . , scales  = "free_y") + 
    theme_BM() +   theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    scale_alpha_manual(values = c(0.3, 1)) +
    scale_fill_manual(values = COLORS6[4:6]) + 
    theme(legend.position = c(0.4, .98), legend.direction = "horizontal", 
          legend.box = "horizontal",
          axis.text.x = element_text(size= 6.5, angle = 90, vjust = 0, hjust = 1, family = "mono" ),
          legend.box.background = element_rect(colour = "black", fill = NA), 
          legend.background =element_blank(),
          legend.text = element_text(size=7),
          legend.title = element_blank(),
          strip.background = element_blank(), 
          strip.text.y = element_blank(),
          legend.key.size = unit(7, "points")) +
    xlab("") + ylab("Mutation count")
}

F5B = plot_profile_absolute(POLE_muts) +
  labs(alpha = "nucs at pos -3-4", fill = NULL) + 
  ggpp::geom_text_npc(data = label_df, npcx = 0.9, npcy = 0.8,
                      aes(label = cosine), size = 3, hjust = 1) +
  ggpp::geom_text_npc(data = label_df, npcx = 0.02, npcy = 0.8,
                      aes(label = name), size = 3.5, hjust = 0)
  
# Plot contribution TTT
SBS28_tri = TRIPLETS_96[c(96)]
SBS28_counts = fread("C:/Users/Axel Rosendahl Huber/Documents/Shared/vanBoxtelLab (Groupfolder)/Projects/Axel/Nissle/Analysis/New_HMF/SBS_28trinuc_counts.txt")
SBS28_counts = slice(SBS28_counts, match(select_HMF, sample)) # only take the select_HMF samples
SBS28_counts$fraction_AA_pos28 = SBS28_counts$TN_AA_sbs28/SBS28_counts$total_muts
SBS28_counts$fraction_AA_NOTpos28 = SBS28_counts$TN_AA_noSBS28/SBS28_counts$total_muts
SBS_28= cbind(SBS28_counts, metadata)
SBS_28$categories = ifelse(SBS_28$log_p_wgs > 3 & SBS_28$AA > 0.22, "motif+", "motif-")
SBS_28$categories[which(SBS_28$sampleId %in% POLE_samples)] = "POLE hotspot"

F5D = ggplot(SBS_28, aes(y = fraction_AA_pos28, x = categories, fill = categories)) +
  geom_violin() +
  geom_jitter(width = 0.2, alpha = 0.4, size  = 0.3) +
  theme_classic() + ylab("T[T>G]T positions \n with -3-4 AA motif") +
  theme(legend.position = "none") + xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

F5E = ggplot(SBS_28, aes(y = fraction_AA_NOTpos28, x = categories, fill = categories)) +
  geom_violin() +
  geom_jitter(width = 0.2, alpha = 0.4, size = 0.3 ) + 
  theme_classic() + ylab("non T[T>G]T positions \n with -3-4 AA motif") + 
  theme(legend.position = "none") +  xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Generate figure 5
layout = '
AABBB
CDBBB'

fig5 = F5A + F5B + F5D + F5E + 
  plot_layout(design = layout) +
  plot_annotation(tag_levels = "A") +
  theme(plot.tag = element_text(face = 'bold'))
ggsave("C:/Users/Axel Rosendahl Huber/OneDrive/Nissle_manuscript/Nissle_September23/Figures/Fig_S7/Fig_S7.pdf", 
       fig5, width = 10, height = 8)
ggsave("C:/Users/Axel Rosendahl Huber/OneDrive/Nissle_manuscript/Nissle_September23/Figures/Fig_S7/Fig_S7.png", 
       fig5, width = 10, height = 8)

