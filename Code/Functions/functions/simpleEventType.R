library(VariantAnnotation)
library(StructuralVariantAnnotation)
library(stringr)

#' Simple SV type classifier
simpleEventType <- function(gr) {
  pgr = partner(gr)
  return(ifelse(seqnames(gr) != seqnames(pgr), "CTX", # inter-chromosomosal
                ifelse(strand(gr) == strand(pgr), "INV",
                       ifelse(gr$insLen >= abs(gr$svLen) * 0.7, "INS", # TODO: improve classification of complex events
                              ifelse(xor(start(gr) < start(pgr), strand(gr) == "-"), "DEL",
                                     "DUP")))))
}
