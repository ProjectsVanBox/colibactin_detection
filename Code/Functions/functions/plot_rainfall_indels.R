# Written by Axel Rosendahl Huber
# 22-07-19 
# Plot a rainfall plot for indels 

plot_rainfall_indels = function (gr_exp, chromosomes, title = "", colors, cex = 2.5, cex_text = 3, 
                                 ylim = 1e+08) 
{
  if (missing(colors)) 
    colors = c("#FDBE6F", "#FF8001", "#B0DD8B", "#36A12E", "#FDCAB5", "#FC8A6A", "#F14432", "#BC141A", "#D0E1F2", "#94C4DF", "#4A98C9", "#1764AB", "#E2E2EF", "#B6B6D8", "#8683BD", "#61409B")
  if (length(colors) != 16) 
    stop("colors vector length not 16")
  
  chr_length = seqlengths(gr_exp)
  chr_length = chr_length[names(chr_length) %in% chromosomes]
  chr_cum = c(0, cumsum(as.numeric(chr_length)))
  names(chr_cum) = names(chr_length)
  labels = gsub("chr", "", names(chr_length))
  m = c()
  for (i in 2:length(chr_cum)) m = c(m, (chr_cum[i - 1] + chr_cum[i])/2)
  type = loc = dist = chrom = c()
  for (i in 1:length(chromosomes)) {
    chr_subset = gr_exp[seqnames(gr_exp) == chromosomes[i]]
    n = length(chr_subset)
    if (n <= 1) {
      next
    }
    type = c(type, chr_subset$muttype_sub[-1])
    loc = c(loc, (start(chr_subset) + chr_cum[i])[-1])
    dist = c(dist, diff(start(chr_subset)))
    chrom = c(chrom, rep(chromosomes[i], n - 1))
  }
  
  muttype = c("C_deletion","T_deletion","C_insertion","T_insertion","2bp_deletion","3bp_deletion","4bp_deletion","5+bp_deletion","2bp_insertion","3bp_insertion","4bp_insertion","5+bp_insertion","2bp_deletion_with_microhomology","3bp_deletion_with_microhomology","4bp_deletion_with_microhomology","5+bp_deletion_with_microhomology")
  
  data = data.frame(type = as.factor(muttype[type + 1]), location = loc, distance = dist, chromosome = chrom)
  
  typesin = muttype %in% levels(data$type)
  colors = colors[typesin]
  location = NULL
  
  plot = ggplot(data, aes(x = location, y = distance)) + geom_point(aes(colour = factor(type)), 
                                                                    cex = cex) + geom_vline(xintercept = as.vector(chr_cum), 
                                                                                            linetype = "dotted") + annotate("text", x = m, y = ylim, 
                                                                                                                            label = labels, cex = cex_text) + xlab("Genomic Location") + 
    ylab("Genomic Distance") + scale_y_log10() + scale_colour_manual(values = colors) + 
    scale_x_continuous(expand = c(0, 0), limits = c(0, max(chr_cum))) + 
    ggtitle(title) + theme_bw() + theme(legend.title = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(), 
                                        axis.ticks.x = element_blank(), axis.text.x = element_blank()) + 
    guides(colour = guide_legend(nrow = 16))
  return(plot)
}
