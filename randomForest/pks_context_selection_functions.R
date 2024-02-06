# ----- Functions to determine pks-patterns in whole-genome cancer driver data
library(GenomicRanges)
mutationDataToGranges = function(mutation_data) {
  # function to convert cBioPortal mutation tables to granges format
  mut_table = data.frame(chr = paste0("chr",mutation_data$chr), 
                         start = mutation_data$start_position, 
                         end = mutation_data$end_position,
                         REF = mutation_data$reference_allele, 
                         ALT = mutation_data$variant_allele)
  
  granges = makeGRangesFromDataFrame(mut_table, ignore.strand = T, keep.extra.columns = T)
  seqlevelsStyle(granges) = "UCSC"
  seqlevels(granges, pruning.mode ="coarse") = seqlevels(granges)[seqlevels(granges) != "chrNA"]
  
  return(granges)
}

select_snvs = function(mutation_gr) {
  # select only snvs in the indel
  mutation_gr = remove_mult_alts(mutation_gr)
  mutation_gr_snv = mutation_gr[nchar(mutation_gr$REF) == 1 &
                                  nchar(unlist(mutation_gr$ALT)) == 1 &
                                  !is.na(seqnames(mutation_gr)) & start(ranges(mutation_gr)) != -1 & seqnames(mutation_gr) != "chrNA" &
                                  mutation_gr$REF != "-" &
                                  unlist(mutation_gr$ALT) != "-", ]
  
  strand(mutation_gr_snv) = ifelse(mutation_gr_snv$REF == "C" | mutation_gr_snv$REF == "T", "+", "-")
  
  return(mutation_gr_snv)
}

select_context_snv = function(gr) {
  gr = remove_mult_alts(gr)
  chr = seqnames(gr)
  start = start(gr) - 10 %>% as.integer()
  end = end(gr) + 10 %>% as.integer()
  strand = strand(gr) %>% as.character()
  contexts = getSeq(BSgenome.Hsapiens.UCSC.hg38, names = chr, start =start , end = end ,strand = strand, as.character = T )
  
  gr$context = contexts
  
  gr <- gr[gr$REF == 'T' | gr$REF == 'A']
  
  return(gr)
}

remove_mult_alts = function(gr) {
  
  mult_alts = elementNROWS(gr$ALT) > 1
  nr_mult_alts = sum(mult_alts)
  if (nr_mult_alts > 0){
    gr = gr[!mult_alts]
    warning(paste0("Removed ", nr_mult_alts, " variant(s), because they had multiple alternative alleles."), call. = F)
  }
  return(gr)
}

search_homlen <- function( x ) {
  offset <- 0
  for ( c1 in c(10:1)) {
    if ( substring(x,c1,c1) == "T" ) {
      offset <- offset+1
    } else {
      break
    }
  }
  
  # for ( c2 in c(12:21)) {
  #   if ( substring(x,c2,c2) == "A" ) {
  #     offset <- offset+1
  #   } else {
  #     break
  #   }
  # }
  return( offset )
}

select_context_indel = function(gr, type = 'cBioPortal'){
  
  gr_A = GRanges()
  gr_T = GRanges()
  
  if (type == "cBioPortal") {
    gr_dels = gr[(gr$REF == "T" | gr$REF == "A") & unlist(gr$ALT) == "-" ,]
    if (length(gr_dels) == 0 ) {return()}
    chr = seqnames(gr_dels )
    start = start(gr_dels ) - 10 %>% as.integer()
    end = end(gr_dels ) + 10 %>% as.integer()
    strand = strand(gr_dels ) %>% as.character()
    contexts = getSeq(BSgenome.Hsapiens.UCSC.hg38, names = chr, start =start , end = end ,strand = strand, as.character = T )
    #contextT = contexts[gr_dels$REF == "T"]
    #contextA = contexts[gr_dels$REF == "A"]
    
    #gr_delsT = gr_dels[gr_dels$REF == "T"]
    #gr_delsA = gr_dels[gr_dels$REF == "A"]
  }
  
  else if (type == "Strelka") {
    # Remove locations with multiple alternative alleles
    gr = remove_mult_alts(gr)
    gr_dels = gr[width(gr$REF) == 2 & width(unlist(gr$ALT)) == 1]
    gr_dels <- gr_dels[grepl("[ACG]T", gr_dels$REF) | grepl("[TCG]A", gr_dels$REF)]
    strand(gr_dels[grepl("[ACG]T", gr_dels$REF)]) <- "+"
    strand(gr_dels[grepl("[TCG]A", gr_dels$REF)]) <- "-"
    if (length(gr_dels) == 0 ) {return()}
    chr = seqnames(gr_dels)
    start = end(gr_dels ) - 10 %>% as.integer()
    end = end(gr_dels ) + 10 %>% as.integer()
    strand = strand(gr_dels ) %>% as.character()
    contexts = getSeq(BSgenome.Hsapiens.UCSC.hg38, names = chr, start =start , end = end ,strand = strand, as.character = T )
    offsets <- as.vector(unlist(lapply(contexts, search_homlen)))
    start2  = start
    end2 = end
    while( sum(offsets) > 0 ) {
      #print( offsets )
      start2  = start2+offsets
      end2 = end2+offsets
      contexts2 = getSeq(BSgenome.Hsapiens.UCSC.hg38, names = chr, start =start2 , end = end2 ,strand = strand, as.character = T )
      offsets <- as.vector(unlist(lapply(contexts2, search_homlen)))
    }
    #print( contexts[3] )
    #print( contexts2[3] )
    
    gr_dels$context = contexts2
    
    print( gr_dels) 
    
    
    
    return( gr_dels )
  }
}


select_context_indel_ORI = function(gr, type = 'cBioPortal'){
  
  gr_A = GRanges()
  gr_T = GRanges()
  
  if (type == "cBioPortal") {
    gr_dels = gr[(gr$REF == "T" | gr$REF == "A") & unlist(gr$ALT) == "-" ,]
    if (length(gr_dels) == 0 ) {return()}
    chr = seqnames(gr_dels )
    start = start(gr_dels ) - 10 %>% as.integer()
    end = end(gr_dels ) + 10 %>% as.integer()
    strand = strand(gr_dels ) %>% as.character()
    contexts = getSeq(BSgenome.Hsapiens.UCSC.hg38, names = chr, start =start , end = end ,strand = strand, as.character = T )
    
    contextT = contexts[gr_dels$REF == "T"]
    contextA = contexts[gr_dels$REF == "A"]
    
    gr_delsT = gr_dels[gr_dels$REF == "T"]
    gr_delsA = gr_dels[gr_dels$REF == "A"]
  }
  
  else if (type == "Strelka") {
    # Remove locations with multiple alternative alleles
    gr = remove_mult_alts(gr)
    gr_dels = gr[width(gr$REF) == 2 & width(unlist(gr$ALT)) == 1]
    if (length(gr_dels) == 0 ) {return()}
    chr = seqnames(gr_dels)
    start = end(gr_dels ) - 10 %>% as.integer()
    end = end(gr_dels ) + 10 %>% as.integer()
    strand = strand(gr_dels ) %>% as.character()
    contexts = getSeq(BSgenome.Hsapiens.UCSC.hg38, names = chr, start =start , end = end ,strand = strand, as.character = T )
    
    contextT = contexts[grepl("[ACG]T", gr_dels$REF)]
    contextA = contexts[grepl("[TCG]A", gr_dels$REF)]
    
    gr_delsT = gr_dels[grepl("[ACG]T", gr_dels$REF)]
    gr_delsA = gr_dels[grepl("[TCG]A", gr_dels$REF)]
    }
  else {stop("Enter input type of vcf: 'cBioPortal'(default) or 'Strelka' ")}

  #print( gr_delsT )
  # return 0 when no single-t deletions can be found
  if (length(gr_delsT) > 0 ) {
    # Phase indels to make sure start of T-deletion takes place at position 9. 
    contextT[substring(contextT, 5,9) == "TTTTT"] = paste0("NNNN", contextT[substring(contextT, 5,9) == "TTTTT"])
    contextT[substring(contextT, 6,9) == "TTTT"] = paste0("NNN", contextT[substring(contextT, 6,9) == "TTTT"])
    contextT[substring(contextT, 7,9) == "TTT"] = paste0("NN", contextT[substring(contextT, 7,9) == "TTT"])
    contextT[substring(contextT, 8,9) == "TT"] = paste0("N", contextT[substring(contextT, 8,9) == "TT"])
    
    
    gr_delsT$context <- contextT
    
    # count occurrences of T-deletions with the pks-context
    Tcontexts = gr_delsT[grepl("T[ACG]", substring(contextT, 9,10)) & substring(contextT, 5,8) == "AAAA"]
    TTcontexts =  gr_delsT[grepl("TT[ACG]", substring(contextT, 9,11)) & substring(contextT, 6,8) == "AAA"]
    TTTcontexts =  gr_delsT[grepl("TTT[ACG]", substring(contextT, 9,12)) & substring(contextT, 7,8) == "AA"]
    TTTTcontexts  = gr_delsT[grepl("TTTT[ACG]", substring(contextT, 9,13)) & substring(contextT, 8,8) == "A"]
    T5contexts  = gr_delsT[grepl("TTTTT[ACG]", substring(contextT, 9,14)) & substring(contextT, 8,8) == "A"]
    T6contexts  = gr_delsT[grepl("TTTTTT[ACG]", substring(contextT, 9,15)) & substring(contextT, 8,8) == "A"]
    T7contexts  = gr_delsT[grepl("TTTTTTT[ACG]", substring(contextT, 9,16)) & substring(contextT, 8,8) == "A"]
    T8contexts  = gr_delsT[grepl("TTTTTTTT[ACG]", substring(contextT, 9,17)) & substring(contextT, 8,8) == "A"]
    
    
    gr_T = unlist(GRangesList(Tcontexts, TTcontexts, TTTcontexts, TTTTcontexts, T5contexts, T6contexts, T7contexts, T8contexts ))
  }
  
  # Phase indels to make sure start of T-deletion takes place at 1
  if (length(gr_delsA) > 0 ) {
    contextA[substring(contextA, 5,9) == "AAAAA"] = paste0("NNNN", contextA[substring(contextA, 5,9) == "AAAAA"])
    contextA[substring(contextA, 6,9) == "AAAA"] = paste0("NNN", contextA[substring(contextA, 6,9) == "AAAA"])
    contextA[substring(contextA, 7,9) == "AAA"] = paste0("NN", contextA[substring(contextA, 7,9) == "AAA"])
    contextA[substring(contextA, 8,9) == "AA"] = paste0("N", contextA[substring(contextA, 8,9) == "AA"])
    
    gr_delsA$context <- contextA
    # Count all mutations matching SBS-pks pattern
    information_Adel = integer(length(contextA))
    
    for (i in 1:length(contextA)) {
      context = contextA[i]
      Arepeat = substring(context, 9, 17)
      rle = rle(strsplit(Arepeat,"")[[1]])
      length_As = rle$lengths[1]
      
      if (length_As == 9) {information_Adel[i] = 0 }
      else if (length_As >= 5 & substring(Arepeat, length_As+1, length_As+1) == "T") {
        information_Adel[i] = 5 }
      else if (length_As == 4 & substring(Arepeat, length_As+1, length_As+1) == "T") {
        information_Adel[i] = 4 }
      else if (length_As == 3 & substring(Arepeat, length_As+1, length_As+2) == "TT") {
        information_Adel[i] = 3 }
      else if (length_As == 2 & substring(Arepeat, length_As+1, length_As+3) == "TTT") {
        information_Adel[i] = 2 }
      else if (length_As == 1 & substring(Arepeat, length_As+1, length_As+4) == "TTTT") {
        information_Adel[i] = 1 }
      else information_Adel[i] = 0
    }
  } else (return())
  
  index = ifelse(information_Adel > 0, T, F)
  gr_A = gr_delsA[index]
  
  print(gr_A)
  print(gr_T)
  gr_pks_dels = unlist(GRangesList(gr_A, gr_T))
  
  return(gr_pks_dels)
} 
