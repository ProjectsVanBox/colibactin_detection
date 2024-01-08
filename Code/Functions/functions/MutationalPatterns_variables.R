# Axel Rosendahl Huber 5-1-2018
# Functions to correctly modify all the mutational signatures variables
# Global package variables

# Default colours for mutation spectrum plotting
COLORS6 = c(
  "#2EBAED", "#000000", "#DE1C14",
  "#D4D2D2", "#ADCC54", "#F0D0CE")

COLORS7 = c(
  "#2EBAED", "#000000", "#DE1C14",
  "#E98C7B", "#D4D2D2", "#ADCC54",
  "#F0D0CE")

SUBSTITUTIONS = c('C>A','C>G','C>T','T>A','T>C','T>G')
SUBSTITUTIONS_96 = rep(SUBSTITUTIONS, each=16)
SUBSTITUTIONS_192 = rep(SUBSTITUTIONS, each=32)
SUBSTITUTIONS_1536 = rep(SUBSTITUTIONS, each=256)


C_TRIPLETS = c(
  "ACA", "ACC", "ACG", "ACT",
  "CCA", "CCC", "CCG", "CCT",
  "GCA", "GCC", "GCG", "GCT",
  "TCA", "TCC", "TCG", "TCT")

T_TRIPLETS = c(
  "ATA", "ATC", "ATG", "ATT",
  "CTA", "CTC", "CTG", "CTT",
  "GTA", "GTC", "GTG", "GTT",
  "TTA", "TTC", "TTG", "TTT")

CONTEXTS_96 = c(rep(C_TRIPLETS, 3), rep(T_TRIPLETS, 3))

# combine substitutions and context in one 
TRIPLETS_96 = paste(substr(CONTEXTS_96,1,1), "[", SUBSTITUTIONS_96, "]", substr(CONTEXTS_96,3,3), sep = "")
STRAND = rep(c("U","T"), 96)
DNA_BASES = c("A", "C", "G", "T")
