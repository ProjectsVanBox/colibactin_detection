indel_strand = function (vcf, ranges) {
  
  TC_positions = substring(vcf$muttype, 1,1)
  vcf = vcf[grep("T|C", TC_positions) ,]
  
  genes = GenomicRanges::reduce(ranges)
  if (!(all(seqlevels(vcf) %in% seqlevels(genes)))) 
    stop(paste("Chromosome names (seqlevels) of vcf and genes Granges", 
               "object do not match. Use the seqlevelsStyle() function", 
               "to rename chromosome names."))
  
  overlap = findOverlaps(vcf, genes)
  overlap = as.data.frame(as.matrix(overlap))
  colnames(overlap) = c("vcf_id", "gene_body_id")
  dup_pos = overlap$vcf_id[duplicated(overlap$vcf_id)]
  dup_idx = which(overlap$vcf_id %in% dup_pos)
  
  if (length(dup_idx) > 0) 
    overlap = overlap[-dup_idx, ]
  
  vcf_overlap = vcf[overlap$vcf_id]
  muttype = vcf_overlap$muttype
  
  # Define which strands are untranscribed 
  strand_muts = rep(0, nrow(overlap))
  for (index in 1:length(strand_muts)) {
    muttype = vcf_overlap$muttype[index]
    if (muttype == "T_deletion" | muttype == "C_deletion") {
      refbase = substring(vcf_overlap$REF[index], 2,2)
      strand_muts[index] = ifelse(refbase == "C" | refbase == "T", "+", "-")
    } 
    if (vcf_overlap$muttype[index] == "T_insertion" | vcf_overlap$muttype[index] == "C_insertion") {
      altbase = substring(vcf_overlap$ALT[index][[1]], 2,2)
      strand_muts[index] = ifelse(altbase == "C" | altbase == "T", "+", "-")
    }
  }
  
  strand_genebodies = as.character(strand(genes)[overlap$gene_body_id])
  same_strand = (strand_muts == strand_genebodies)
  U_index = which(same_strand == TRUE)
  T_index = which(same_strand == FALSE)
  strand = rep(0, nrow(overlap))
  strand[U_index] = "untranscribed"
  strand[T_index] = "transcribed"
  
  strand2 = rep("-", length(vcf))
  strand2[overlap$vcf_id] = strand
  strand2 = factor(strand2, levels = c("untranscribed", 
                                       "transcribed", "-"))
  vcf_overlap$TRANSCRIPTION = strand
  return(vcf_overlap)
  }
  
