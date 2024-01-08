# Axel Rosendahl Huber 
# 25-03-2019
# Modification of the mut_matrix function to deal with the mclapply errors. 

mut_96_occurrences2 = function(type_context)
{
    vector = rep(0,96)
    names(vector) = TRIPLETS_96
    
    # if type_context is empty, return vector with zeroes
    if (isEmpty(type_context$types) || isEmpty(type_context$context))
        return(vector)
    
    # for all mutations in this sample
    for (i in 1:length(type_context[[1]]))
    {
        # Find mutation type
        type = which(SUBSTITUTIONS == type_context[[1]][i])

        # Find triplet
        if(type < 4)
            context = which(C_TRIPLETS == type_context[[2]][i])
        else
            context = which(T_TRIPLETS == type_context[[2]][i])

        pos = (type - 1)*16 + context
        vector[pos] = vector[pos] + 1
    }

    return(vector)
}

mut_matrix2 <- function (vcf_list, ref_genome) 
{
    df = data.frame()
    rows = lapply(as.list(vcf_list), function(vcf) {
        type_context = type_context(vcf, ref_genome)
        row = mut_96_occurrences2(type_context)
        return(row)})
    for (row in rows) {
        if (class(row) == "try-error") 
            stop(row)
        df = rbind(df, row)
    }
    names(df) = names(row)
    row.names(df) = names(vcf_list)
    return(t(df))
}
