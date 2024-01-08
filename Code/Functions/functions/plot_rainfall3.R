plot_rainfall3 = function (vcf, chromosomes, title = "", colors, cex = 2.5, cex_text = 3, 
          ylim = 1e+08) 
{
  if (missing(colors)) 
    colors = COLORS6
  if (length(colors) != 6) 
    stop("colors vector length not 6")
  chr_length = seqlengths(vcf)
  chr_length = chr_length[names(chr_length) %in% chromosomes]
  chr_cum = c(0, cumsum(as.numeric(chr_length)))
  names(chr_cum) = names(chr_length)
  labels = gsub("chr", "", names(chr_length))
  m = c()
  for (i in 2:length(chr_cum)) m = c(m, (chr_cum[i - 1] + chr_cum[i])/2)
  type = loc = dist = chrom = c()
  for (i in 1:length(chromosomes)) {
    chr_subset = vcf[seqnames(vcf) == chromosomes[i]]
    n = length(chr_subset)
    if (n <= 1) {
      next
    }
    type = c(type, mut_type(chr_subset)[-1])
    loc = c(loc, (start(chr_subset) + chr_cum[i])[-1])
    dist = c(dist, diff(start(chr_subset)))
    chrom = c(chrom, rep(chromosomes[i], n - 1))
  }
  data = data.frame(type = type, location = loc, distance = dist, 
                    chromosome = chrom)
  typesin = SUBSTITUTIONS %in% unique(data$type)
  colors = colors[typesin]
  location = NULL
  plot = ggplot(data, aes(x = location, y = distance)) + geom_point(aes(colour = factor(type)), 
                                                                    cex = cex) + geom_vline(xintercept = as.vector(chr_cum), 
                                                                                            linetype = "dotted") + annotate("text", x = m, y = ylim, 
                                                                                                                            label = labels, cex = cex_text) + xlab("Genomic Location") + 
    ylab("Genomic Distance") + scale_y_log10() + scale_colour_manual(values = colors) + 
    scale_x_continuous(expand = c(0, 0), limits = c(0, max(chr_cum))) + 
    ggtitle(title) + theme_bw() + theme(legend.position = "bottom", 
                                        legend.title = element_blank(), legend.key = element_blank(), 
                                        panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(), 
                                        axis.ticks.x = element_blank(), axis.text.x = element_blank()) + 
    guides(colour = guide_legend(nrow = 1))
  return(plot)
}
