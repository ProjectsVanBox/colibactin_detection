plot_192_profile4 = function(mut_matrix,colors, ymax = 5, condensed = TRUE)
{
  number_sigs = dim(mut_matrix)[2]
  if (is.null(number_sigs)) {
    strand = sapply(names(mut_matrix), function(x) strsplit(x, "-")[[1]][2])
  } else if (number_sigs >= 1) {
    strand = sapply(rownames(mut_matrix), function(x) strsplit(x, "-")[[1]][2])
  }
  if (missing(colors)) {
    colors = COLORS6
  }
  if (length(colors) != 6) {
    stop("Provide colors vector with length 6")
  }
  context = rep(CONTEXTS_96, each = 2)
  substitution = rep(SUBSTITUTIONS, each = 32)
  substring(context, 2, 2) = "."
  df = data.frame(substitution = substitution, context = context, strand = strand)
  rownames(mut_matrix) = NULL
  df2 = cbind(df, as.data.frame(mut_matrix))
  seq  = 1:192 + rep(c(1,-1), times = 96)
  df2 = df2[seq,]
  df3 = melt(df2, id.vars = c("substitution", "context", "strand"))  
  df3$strand = factor(df3$strand, levels = c("transcribed", "untranscribed"))
  
  value = NULL
  if (condensed) {
    plot = ggplot(data = df3, aes(x = context, y = value, fill = substitution, width = 1, alpha = strand)) +
      geom_bar(position = "dodge", stat = "identity", colour = "black", size = 0.2) +
      scale_fill_manual(values = colors) + facet_grid(variable ~ substitution) + ylab("Relative contribution") + coord_cartesian(ylim = c(0,  ymax)) +
      scale_y_continuous(expand = c(0, 0), breaks = seq(0, ymax, 0.5)) +
      expand_limits(y = 0) +  scale_alpha_manual(values = c(1,.2)) +
      guides(fill = FALSE) + theme_classic() +
      theme(axis.title.y = element_text(size = 12, vjust = 1),
            axis.text.y = element_text(size = 8),
            axis.title.x = element_text(size = 12),
            axis.text.x = element_text(size = 5, angle = 90, vjust = 0.4), 
            strip.text.x = element_text(size = 9),
            strip.text.y = element_text(size = 9),
            panel.grid.major.x = element_blank(),
            panel.spacing.x = unit(0, "lines")
      )
  }
  else {
    plot = ggplot(data = df3, aes(x = context, y = value, fill = substitution, width = 0.6, alpha = strand)) +
      geom_bar(stat = "identity", colour = "black", size = 0.2) +
      scale_fill_manual(values = colors) + facet_grid(variable ~ substitution) +
      ylab("Relative contribution") + coord_cartesian(ylim = c(0, ymax)) + 
      scale_y_continuous(breaks = seq(0, ymax, 0.1)) +
      guides(fill = FALSE) + theme_bw() + theme( axis.title.y = element_text(size = 12, vjust = 1),
                                                 axis.text.y = element_text(size = 8),
                                                 axis.title.x = element_text(size = 12),
                                                 axis.text.x = element_text(size = 5, angle = 90, vjust = 0.4),
                                                 strip.text.x = element_text(size = 9),
                                                 strip.text.y = element_text(size = 9),
                                                 panel.grid.major.x = element_blank()
      )
  }
  return(plot)
}
