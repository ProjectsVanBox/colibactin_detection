# scripts for to analyze all hg38 microinjection data for pks+ bacteria (CCR E.coli and Nissle E. coli)
library(BSgenome)
library(MutationalPatterns)
library(tidyverse)
library(plyr)
library(data.table)
library(vroom)
library(ggpubr)
library(cowplot)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# update script to get the correct working directory
setwd("C:/Users/Axel Rosendahl Huber/OneDrive/Nissle_manuscript/Nissle")

source("Code/Functions/Utils.R")
source("Code/Functions/pks_context_selection_functions.R")
source("Code/Functions/Nissle_functions.R")
ref_genome = "BSgenome.Hsapiens.UCSC.hg38"
genes = genes(TxDb.Hsapiens.UCSC.hg38.knownGene)

# --------------------------------
get_context = function(gr, size_context = 10){
  gr = gr[gr$FILTER == "PASS"]
  gr = gr[which(nchar(gr$REF) == 1 )]
  strand = ifelse(gr$REF == "G" | gr$REF == "A", '-', "+")
  start = start(gr)
  ref = as.character(gr$REF)
  alt = unlist(CharacterList(gr$ALT))
  type = paste0(ref, ">", alt)
  chromosome = as.character(seqnames(gr))
  context = getSeq(Hsapiens, chromosome,  start = start - size_context,
                   end =  start + size_context, 
                   strand = strand)
      
  type = mapvalues(type, c("A>C","A>G","A>T","G>C", "G>A","G>T"), 
                   c("T>G", "T>C", "T>A", "C>G", "C>T", "C>A"), warn_missing	= FALSE)
  trinucleotide = paste0(substr(context, 10, 10), "[", type, "]", substr(context, 12,12))
  
  
  context_table = tibble(chr = chromosome, position = start, type = type, strand = strand, 
                         context = as.character(context), trinucleotide = trinucleotide)
  context_table$id = paste0(context_table$chr, "_", context_table$position, "_", context_table$type)
  return(context_table)
}

vcf_files =  list.files("Data/", recursive = T, full.names = T, pattern = "\\.vcf|.ptato")
vcf_files = vcf_files[!grepl("PTA_indels", vcf_files)]
vcf_indel_files = list.files("Data/PTA_indels/", full.names = T)
names(vcf_indel_files) = gsub(".*_|.vcf|.indels.ptato.|filtered|vcf", "", basename(vcf_indel_files))
names(vcf_indel_files) = gsub("-", ".", names(vcf_indel_files))
####################################
# Process Exposed organoid mutation data 
####################################
names(vcf_files) = gsub(".*_|.vcf|.snvs.ptato.filtered|vcf", "", basename(vcf_files))
names(vcf_files) = gsub("-", ".", names(vcf_files))
vcf_files = vcf_files[!grepl("TM0", names(vcf_files))]

vcfs = read_vcfs_as_granges(vcf_files, names(vcf_files), genome =  ref_genome, type = "all")
vcfs_sbs = get_mut_type(vcfs, "snv")
vcfs_indel = get_mut_type(vcfs, "indel")

vcfs_indel[1:25] = read_vcfs_as_granges(vcf_indel_files, names(vcf_indel_files), ref_genome, type = "indel")

mut_mat = mut_matrix(vcfs_sbs, ref_genome)
mut_mat_s = mut_matrix_stranded(vcfs_sbs, ref_genome, genes)

# get mutation loads 
indel_loads = lengths(vcfs_indel) %>% as.data.frame()
colnames(indel_loads) = "total_indels"
id_contexts = MutationalPatterns::get_indel_context(vcfs_indel, ref_genome)
id_pks_contexts = lapply(id_contexts, select_context_indel, type = "Strelka")
indel_counts = count_indel_contexts(id_contexts)

indel_loads$in_pks_motif = lengths(id_pks_contexts)
indel_loads$fraction_pksmotif = indel_loads$in_pks_motif/lengths(id_contexts)

# ===== Load different names ===== # 
categories = rep("EcC", ncol(mut_mat))
categories[grepl("I3NIS|5N", colnames(mut_mat))] = "EcN"
categories[grepl("19H2", colnames(mut_mat))] = "19H2"
categories[grepl("2F8", colnames(mut_mat))] = "2F8"
categories[grepl("DYE", colnames(mut_mat))] = "Control"
categories[grepl("EKO", colnames(mut_mat))] = "Control"

categories = data.table("injection" = categories)
categories = categories %>% 
  mutate(name = colnames(mut_mat)) %>% 
  mutate(method = ifelse(grepl("NISL", name), "PTA","Clonal Expansion"))
categories$injection = factor(categories$injection, levels = c("Control", "EcN", "EcC", "19H2", "2F8"))

# Load dinucleotide categories
nucs = c("A", "T", "C", "G")
dinucs = expand.grid(nucs, nucs)
dinucs = paste0(dinucs[,1],dinucs[,2])

####################################
# Extended mutational contexts
####################################
contexts = lapply(vcfs_sbs, get_context)


# create the neccasary processed data folders
if (!dir.exists("Processed_data")) { dir.create("Processed_data")}
if (!dir.exists("Processed_data/Contexts")) { dir.create("Processed_data/Contexts")}


# save unique context files 
context_list = list(EcC = contexts[categories$injection == "EcC"],
                    EcN = contexts[categories$injection == "EcN"],
                    `19H2` = contexts[categories$injection == "19H2"],
                    `2F8` = contexts[categories$injection == "2F8"],
                    Control = contexts[categories$injection == "Control"]) 

for (name in names(context_list)) {
  ctx = context_list[[name]]
  ctx = rbindlist(ctx)
  ctx = dplyr::distinct(ctx)
  ctx = ctx[substr(ctx$type, 1,1) == "T","context"]
  write.table(ctx, file = paste0("Processed_data/Contexts/",name, "_unique_contexts.txt"), quote = F, col.names = F, row.names = F)
}

####################################
# Hartwig mutational re-fitting data 
####################################

# NOTE: This part of the script uses patient-level somatic variant and clinical data which have been obtained from the Hartwig Medical Foundation under the data request number DR-084. Somatic variant and clinical data are freely available for academic use from the Hartwig Medical Foundation through standardized procedures. Privacy and publication policies, including co-authorship policies, can be retrieved from: https://www.hartwigmedicalfoundation.nl/en/data-policy/. 
# Data request forms can be downloaded from https://www.hartwigmedicalfoundation.nl/en/applying-for-data/.
# To gain access to the data, this data request form should be emailed to info@hartwigmedicalfoundation.nl., upon which it will be evaluated within 6 weeks by the HMF Scientific Council and an independent Data Access Board.
# When access is granted, the requested data become available through a download link provided by HMF.

# id_signatures = read_delim("https://cancer.sanger.ac.uk/signatures/documents/451/COSMIC_v3.2_ID_GRCh37.txt") 
# signatures = read_delim("https://cancer.sanger.ac.uk/signatures/documents/453/COSMIC_v3.2_SBS_GRCh38.txt") %>% 
#   arrange(match(Type, TRIPLETS_96))
# artefact_signatures = c("SBS27", "SBS43","SBS45","SBS46", "SBS47","SBS48","SBS49","SBS50","SBS51","SBS52",'SBS53',"SBS54",'SBS55','SBS56','SBS57','SBS58','SBS59','SBS60')
# sigs = signatures[,!colnames(signatures) %in% artefact_signatures] # remove signatures marked as artefacts
# 
# # load HMF data (exome and whole genome) 
# HMF mut mat data can be generated using the read_vcfs_as_granges function from the package MutationalPatterns
# mm_HMF = read.delim("C:/Users/Axel Rosendahl Huber/Documents/Shared/vanBoxtelLab (Groupfolder)/Projects/Axel/Nissle/Analysis/New_HMF/20210614_HMF_sbs_matrix_somatics.tsv") 
# mm_HMF_exome = fread("C:/Users/Axel Rosendahl Huber/Documents/Shared/vanBoxtelLab (Groupfolder)/Projects/Axel/Nissle/Analysis/New_HMF/mm_exome.txt", data.table = F)
# select_HMF =  colnames(mm_HMF)[colSums(mm_HMF) > 100 ]  # remove all samples with mutation counts < 100
# mm_HMF = mm_HMF[, select_HMF]
# 
# ids_HMF = read.delim("C:/Users/Axel Rosendahl Huber/Documents/Shared/vanBoxtelLab (Groupfolder)/Projects/Axel/Nissle/Analysis/New_HMF/20210614_HMF_indel_matrix_somatics.tsv")[,select_HMF]
# ids_HMF_exome = fread("C:/Users/Axel Rosendahl Huber/Documents/Shared/vanBoxtelLab (Groupfolder)/Projects/Axel/Nissle/Analysis/New_HMF/indel_counts_HMF_exome.txt", data.table = F)
# 
# # remove all samples with no mutation counts in the exome in both indels and snvs
# sbsExSelect = colnames(mm_HMF_exome)[colSums(mm_HMF_exome) > 0]
# indelExSelect = colnames(ids_HMF_exome)[colSums(ids_HMF_exome) > 0]
# exome_select = intersect(intersect(sbsExSelect, indelExSelect), select_HMF)
# mm_HMF_exome = mm_HMF_exome[,exome_select]
# ids_HMF_exome = ids_HMF_exome[,exome_select]
# 
# # HMF WGS refits 
# sbs_fit = fit_to_signatures(mm_HMF %>% as.matrix(), signatures[,-1] %>% as.matrix() ) 
# sbs_contri = prop.table(sbs_fit$contribution,2)
# 
# id_fit = fit_to_signatures(ids_HMF, id_signatures[,-1] %>% as.matrix())
# id_contri = prop.table(id_fit$contribution,2)
# 
# # HMF exome refits
# sbs_fit_exome = fit_to_signatures(mut_matrix = mm_HMF_exome %>% as.matrix(), signatures = signatures[,-1] %>% as.matrix() ) 
# sbs_contri_exome = prop.table(sbs_fit_exome$contribution,2)
# sbs_contri_exome %>% colSums() %>% is.na() %>% table()
# 
# id_fit_exome = fit_to_signatures(ids_HMF_exome, id_signatures[,-1] %>% as.matrix())
# id_contri_exome = prop.table(id_fit_exome$contribution,2)

#########
# Hartwig dinucleotide enrichments
##########
# combine the results with signature extraction information: 
# below a code block so the signature analysis is performed like the analysis in the Nature paper 
metadata = read.delim("C:/Users/Axel Rosendahl Huber/Documents/Shared/vanBoxtelLab (Groupfolder)/Projects/Axel/Nissle/Analysis/New_HMF/metadata.tsv")
metadata = metadata[match(select_HMF, metadata$sampleId,),] # reorder metadata rows in the same manner as the dinucleotide matrices
metadata$n_sbs =  colSums(mm_HMF)[match(metadata$sampleId, colnames(mm_HMF))] # get total sbs load
metadata$n_indels =  colSums(ids_HMF)[match(metadata$sampleId, colnames(ids_HMF))] # get total indel load

dinuc_mat = read.delim("C:/Users/Axel Rosendahl Huber/Documents/Shared/vanBoxtelLab (Groupfolder)/Projects/Axel/Nissle/Analysis/New_HMF/dinuc_mat_hmf_sbs_peaks.txt")
dinuc_mat = dinuc_mat[select_HMF,]
dinuc_exome = read.delim("C:/Users/Axel Rosendahl Huber/Documents/Shared/vanBoxtelLab (Groupfolder)/Projects/Axel/Nissle/Analysis/New_HMF/dinuc_mat_exome_sbs_peaks.txt")
dinuc_exome  = dinuc_exome[exome_select,]
# p_value exome: fisher test against all other samples 
log_p_wgs_dinucs = get_dinuc_enrichment(dinuc_mat)
metadata$log_p_wgs = log_p_wgs_dinucs
metadata$log_p_all_exome =  0 
dinuc_enrichment = get_dinuc_enrichment(dinuc_exome)
rownames(metadata) = metadata$sampleId
metadata[names(dinuc_enrichment), "log_p_all_exome"] = dinuc_enrichment

save.image(file = paste0("Processed_data/Nissle_processed", ".RData"))

