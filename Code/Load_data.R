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
source("Code/Functions/Nissle_functions.R")
ref_genome = "BSgenome.Hsapiens.UCSC.hg38"
load("Processed_data/Nissle_processed.RData")