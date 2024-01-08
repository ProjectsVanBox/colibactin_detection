# Axel Rosendahl Huber
# 25-03-2019
# Prevent mclapply errors

mut_192_occurrences2 = function(type_context, strand)
{
  # get possible strand values
  values = levels(strand)
  
  idx1 = which(strand == values[1])
  idx2 = which(strand == values[2])
  
  # get type context for both vcf subsets
  type_context_1 = lapply(type_context, function(x) x[idx1])
  type_context_2 = lapply(type_context, function(x) x[idx2])
  
  # make 96-trinucleotide count vector per set
  vector1 = mut_96_occurrences2(type_context_1)
  vector2 = mut_96_occurrences2(type_context_2)
  
  # add names
  names_1 = paste(TRIPLETS_96, values[1], sep = "-")
  names_2 = paste(TRIPLETS_96, values[2], sep = "-")
  
  # combine vectors in alternating fashion
  vector = c(rbind(vector1, vector2))
  names = c(rbind(names_1, names_2))
  names(vector) = names

  return(vector)
}




mut_matrix_stranded2 = function (vcf_list, ref_genome, ranges, mode = "transcription") 
{
    df = data.frame()
    if (mode == "transcription") {
        rows <- lapply(as.list(vcf_list), function(vcf) {
            type_context = type_context(vcf, ref_genome)
            strand = mut_strand(vcf, ranges, mode = "transcription")
            row = mut_192_occurrences2(type_context, strand)
            return(row)
        })
        for (row in rows) df = rbind(df, row)
    }
    if (mode == "replication") {
        rows <- lapply(as.list(vcf_list), function(vcf) {
            type_context = type_context(vcf, ref_genome)
            strand = mut_strand(vcf, ranges, mode = "replication")
            row = mut_192_occurrences2(type_context, strand)
            return(row)
        })
        for (row in rows) {
            if (class(row) == "try-error") 
                stop(row)
            df = rbind(df, row)
        }
    }
    names(df) = names(row)
    row.names(df) = names(vcf_list)
    return(t(df))
}
