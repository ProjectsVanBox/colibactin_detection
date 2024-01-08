count_SVs = function(named_SV_vcfs, cores = 7) { 
  
  vcf_list = mclapply(named_SV_vcfs, FUN = function(x) readVcf(x), mc.cores = cores)

  SVs_bed = mclapply(vcf_list, SimpleSV_bed, mc.cores = 7)
  
  SV_counts_list = sapply(SVs_bed, nrow)
  SV_counts_list[sapply(SV_counts_list, is.null)] = 0 
  SV_counts = as.numeric(SV_counts_list)
  names(SV_counts) = names(SV_counts_list)
  
  SV_type_counts = sapply(SVs_bed, function(x) table(x$name))
  
  SV_table = data.frame(sample = names(SV_counts), patient = patients, SV_count = SV_counts,INS = 0,  INV = 0 , DEL = 0, DUP = 0)
  
  for (sample in names(SV_type_counts)) {
    type_counts= SV_type_counts[[sample]]
    sample_table = SV_table[sample,]
    
    if (dim(type_counts)[1] == 0) {next()}
    else {sample_table[,names(type_counts)] = type_counts}
    
    
    SV_table[sample, ] = sample_table
    
  } 
  
  return(SV_table)
}

  
