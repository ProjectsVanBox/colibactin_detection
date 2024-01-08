# Annotate SV granges (modified from https://github.com/PapenfussLab/gridss/blob/master/example/simple-event-annotation.R)
SimpleSV_bed = function(vcf) {
  gr <- breakpointRanges(vcf)
  if (length(gr) > 0) {
    svtype <- simpleEventType(gr)
    gr$SIMPLE_TYPE <- NA_character_
    gr$SIMPLE_TYPE <- svtype
    gr$svLen 
    
    
    # TODO: perform event filtering here
    # By default, GRIDSS is very sensitive but this comes at the cost of a high false discovery rate
    gr <- gr[gr$FILTER == "PASS" & partner(gr)$FILTER == "PASS"] # Remove low confidence calls
    
    simplegr <- gr[simpleEventType(gr) %in% c("INS", "INV", "DEL", "DUP")]
    if (length(simplegr)> 0 ) {
      simplebed <- data.frame(chrom=seqnames(simplegr),
                              # call the centre of the homology/inexact interval
                              start=as.integer((start(simplegr) + end(simplegr)) / 2),
                              end=as.integer((start(partner(simplegr)) + end(partner(simplegr))) / 2),
                              name=simpleEventType(simplegr),
                              score=simplegr$QUAL,
                              length = simplegr$svLen,
                              strand=".")
      # Only select the lower of the two breakends so we don't output everything twice
      simplebed <- simplebed[simplebed$start < simplebed$end,]
      return(simplebed)
      
    } else {  return() }
    
  } else {return()}
}







