intersect_with_region = function(vcf, surveyed, region)
{
  # Number of mutations in vcf file
  n_muts = length(vcf)
  
  # Number of base pairs that were surveyed
  surveyed_length = sum(as.numeric(width(surveyed)))
  
  # Check if chromosome names are the same in the objects
  if (seqlevelsStyle(vcf) != seqlevelsStyle(surveyed))
    stop(paste("The chromosome names (seqlevels) of the VCF and the",
               "surveyed GRanges object do not match."))
  
  if (seqlevelsStyle(region) != seqlevelsStyle(surveyed))
    stop(paste("The chromosome names (seqlevels) of the surveyed and",
               "the region GRanges object do not match."))
  
  # Intersect genomic region and surveyed region
  surveyed_region = intersect(surveyed, region, ignore.strand = TRUE)
  surveyed_region_length = sum(width(surveyed_region))
  
  # Find which mutations lie in surveyed genomic region
  overlap = findOverlaps(vcf, surveyed_region)
  muts_in_region = as.data.frame(as.matrix(overlap))$queryHits
  
  observed = length(muts_in_region)
  prob = n_muts / surveyed_length
  expected = prob * surveyed_region_length
  
  res = data.frame(n_muts,
                   surveyed_length,
                   prob, surveyed_region_length,
                   expected,
                   observed)
  return(res)
}
