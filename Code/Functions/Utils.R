#=== UTILS ===# 
# Axel Rosendahl Huber
# Purpose: Load the commonly used packages and created functions make analysis clear and streamlined 

# ==== Miscellaneous ===== # 
options(stringsAsFactors = F) # Prevent conversion of strings to factors in dataframe
library(tidyverse)
library(reshape2)
library(magrittr)
library(lemon)
# ======== Retrieve and source all functions ===== # 
# Get all available functions in the Scripts/R_functions folder
functions <- list.files("Code/Functions/functions/", pattern = ".R", full.names = TRUE)
for (i in functions) {
  source(i)
  cat("Loading function:", i, "\n" ,  sep = " ")
}