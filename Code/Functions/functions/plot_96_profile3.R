plot_96_profile3 = function (mut_matrix, colors, relative_values = FALSE) {
  number_sigs = dim(mut_matrix)[2]
  
  if(relative_values) {
    mut_matrix = as.matrix(mut_matrix)
    mut_matrix = prop.table(mut_matrix, 2)
  }
    
  
  if (missing(colors)) {
    colors = COLORS6
  }
  if (length(colors) != 6) {
    stop("Provide colors vector with length 6")
  }
  context = CONTEXTS_96
  substitution = rep(SUBSTITUTIONS, each = 16)
  substring(context, 2, 2) = "."
  df = data.frame(substitution = substitution, context = context)
  rownames(mut_matrix) = NULL
  df2 = cbind(df, as.data.frame(mut_matrix))
  df3 = reshape2::melt(df2, id.vars = c("substitution", "context"))
  value = NULL
 
   plot = ggplot2::ggplot(data = df3, aes(x = context, y = value, 
                                  fill = substitution, width = 0.8)) + 
      geom_bar(stat = "identity", size = 0.3) + scale_fill_manual(values = colors) + 
      facet_grid(variable ~ substitution, scales = "free_y") + 
      guides(fill = "none") + theme_BM() + 
      theme(axis.title.y = element_text(size = 12, vjust = 1), 
            axis.text.y = element_text(size = 8), axis.title.x = element_text(size = 12), 
            axis.text.x = element_text(size = 6, angle = 90, vjust = 0.3),
            strip.text.x = element_text(size = 9), 
            strip.text.y = element_text(size = 9),
            panel.grid.major.x = element_blank(),
            panel.spacing.x = unit(0.1, "lines")) +  
      xlab("96-trinucleotide context") + xlab("")
  
  if (relative_values) {
    plot = plot + ylab("Relative contribution")
  } 
  else{ 
    plot = plot + ylab("Number of single base substitutions")
  }
  
  return(plot)
}
