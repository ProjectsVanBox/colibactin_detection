# Load data script Nissle
library(BSgenome)
library(MutationalPatterns)
library(plyr)
library(data.table)
library(vroom)
library(ggpubr)
library(cowplot)
library(patchwork)
library(BSgenome.Hsapiens.UCSC.hg38)
library(tidyverse)
#source("D:/surfdrive/Shared/pmc_vanboxtel/general/2_Bioinformatics/Scripts/pmc_vanboxtel/Axel_BMseq/Utils.R")
source("Code/Functions/Nissle_functions.R")
# source() Nissle fuctions 
ref_genome = "BSgenome.Hsapiens.UCSC.hg38"
load("Processed_data/Nissle_processed.RData")

# # ------- Load data ------ 
# contexts = readRDS("Processed_data/Contexts/context_list.rds") # load sbs context files 
# 
# id_contexts = readRDS("Processed_data/Contexts/id_contexts_grangeslist.rds")
# id_pks_contexts = readRDS("Processed_data/Contexts/id_pks_contexts_grangeslist.rds")
# indel_loads = lengths(id_contexts) %>% as.data.frame()
# colnames(indel_loads) = "total_indels"
# indel_loads$in_pks_motif = lengths(id_pks_contexts)
# indel_loads$fraction_pksmotif = indel_loads$in_pks_motif/indel_loads$total_indels
# 
# id_signatures = read_delim("https://cancer.sanger.ac.uk/signatures/documents/451/COSMIC_v3.2_ID_GRCh37.txt") 
# signatures = read_delim("https://cancer.sanger.ac.uk/signatures/documents/453/COSMIC_v3.2_SBS_GRCh38.txt") %>% 
#   arrange(match(Type, TRIPLETS_96))
# 
# #  ======== Load sbs data  ========= # 
# vcfs = readRDS("Processed_data/Mutation_data/vcfs_grl.rds")
# vcfs_snv = readRDS("Processed_data/Mutation_data/vcfs_sbs.rds")
# 
# mut_mat = read.delim("Processed_data/Mutation_data/mut_mat_hg38.txt")
# mut_mat_s = read.delim("Processed_data/Mutation_data/mut_mat_s_hg38.txt")
# 
# # ======== Load indel data ======== # 
# id_contexts_grl = readRDS("Processed_data/Contexts/id_contexts_grangeslist.rds")
# id_contexts_grl_pks = readRDS("Processed_data/Contexts/id_pks_contexts_grangeslist.rds")
# indel_counts = read.delim("Processed_data/Mutation_data/indel_counts.txt")
# indel_counts = indel_counts[,colnames(mut_mat)]
# 
# # ===== Load different names ===== # 
# categories = rep("EWT", ncol(mut_mat))
# categories[grepl("EKO", colnames(mut_mat))] = "EKO"
# categories[grepl("ETBF", colnames(mut_mat))] = "ETBF"
# categories[grepl("NIS", colnames(mut_mat))] = "NIS"
# categories[grepl("CDT", colnames(mut_mat))] = "CDT"
# categories[grepl("DYE", colnames(mut_mat))] = "DYE"
# 
# injections = c(rep(1,3), rep(8,2), rep(5,9), rep(8,1), rep(3,10), 
#                rep(8,4), rep(3,3), rep(6,3), rep(3,6))
# 
# # Load dinucleotide categories
# nucs = c("A", "T", "C", "G")
# dinucs = expand.grid(nucs, nucs)
# dinucs = paste0(dinucs[,1],dinucs[,2])
# 
# # ==== select only T>N contexts ======= 
# 
# contexts_TN = list()
# for (type in names(contexts)) {
#   ctx_table= contexts[[type]] %>% 
#     rbindlist() %>% 
#     distinct() %>%
#     filter(grepl("^T", type))
#   contexts_TN[[type]] = ctx_table
# }
# 
# 
# load("Processed_data/select_hmf.Rdata", verbose = T)
# 
# # load HMF mutation matrices
# mm_HMF = read.delim("~/surfdrive/Shared/vanBoxtelLab (Groupfolder)/Projects/Axel/Nissle/Analysis/New_HMF/20210614_HMF_sbs_matrix_somatics.tsv")[,select_HMF]
# mm_HMF_exome = fread("~/surfdrive/Shared/vanBoxtelLab (Groupfolder)/Projects/Axel/Nissle/Analysis/New_HMF/mm_exome.txt", data.table = F )[,select_HMF]
# ids_HMF = read.delim("~/surfdrive/Shared/vanBoxtelLab (Groupfolder)/Projects/Axel/Nissle/Analysis/New_HMF/20210614_HMF_indel_matrix_somatics.tsv")[,select_HMF] 
# 
# # load signature contribution matrices (select HMF samples already removed)
# id_contri = read.delim("~/surfdrive/Shared/vanBoxtelLab (Groupfolder)/Projects/Axel/Nissle/Analysis/Processed_data/Sig_fits/HMF_id_contribution.txt")
# sbs_contri = read.delim("~/surfdrive/Shared/vanBoxtelLab (Groupfolder)/Projects/Axel/Nissle/Analysis/Processed_data/Sig_fits/HMF_sbs_contribution.txt")
# id_contri_exome = read.delim("~/surfdrive/Shared/vanBoxtelLab (Groupfolder)/Projects/Axel/Nissle/Analysis/Processed_data/Sig_fits/HMF_exome_id_contribution.txt")
# sbs_contri_exome = read.delim("~/surfdrive/Shared/vanBoxtelLab (Groupfolder)/Projects/Axel/Nissle/Analysis/Processed_data/Sig_fits/HMF_exome_sbs_contribution.txt")
# 
# # load dinucleotide enrichment results 
# dinuc_mat = read.delim("~/surfdrive/Shared/vanBoxtelLab (Groupfolder)/Projects/Axel/Nissle/Analysis/New_HMF/dinuc_mat_hmf_sbs_peaks.txt")[select_HMF,]
# dinuc_exome = read.delim("~/surfdrive/Shared/vanBoxtelLab (Groupfolder)/Projects/Axel/Nissle/Analysis/New_HMF/dinuc_mat_exome_sbs_peaks.txt")[select_HMF,]
# metadata = read.delim("~/surfdrive/Shared/vanBoxtelLab (Groupfolder)/Projects/Axel/Nissle/Analysis/New_HMF/metadata_extended.txt")
# 
# # # for now, assume that all samples missing in the HMF exonic data do not contain exonic indels (need to check)
# # dif = setdiff(colnames(ids_HMF), colnames(ids_HMF_exome))
# # colSums(ids_HMF[,dif])
# # median(colSums(ids_HMF))
# # mean(colSums(ids_HMF))
# # hist(colSums(ids_HMF), breaks = 2000)
# # max(colSums(ids_HMF))
# # (colSums(ids_HMF))
# # dim(ids_HMF)
# # dim(ids_HMF_exome)
# # ids_HMF_exome = ids_HMF_exome[,select_HMF]
# # 
