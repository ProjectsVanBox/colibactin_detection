# spectrum_96 
# Same function as plot_96_profile but with Error bars 
# Axel Rosendahl Huber 25-05-2018

spectrum_96 <- function(mut_mat, ymax = 0.2) {
  norm_mut_mat <- apply(mut_mat, 2, function(x) x/sum(x))
  context <- CONTEXTS_96
  substring(context, 2, 2) = "."
  
  # Create list of total contributions + list SDs
  sd <- list(apply(as.matrix(norm_mut_mat), 1, FUN = sd))
  total <- data.frame(rowSums(norm_mut_mat), sd = sd)
  rownames(total) <- rownames(norm_mut_mat)
  colnames(total) <- c("values", "sd")
  value = NULL
  total$context <- context
  total$fraction <- total$values/sum(total$values)
  total$sdfraction <- total$sd
  total$top <- total$fraction + total$sdfraction
  total$bottom <- total$fraction - total$sdfraction
  total$bottom[total$bottom < 0] <- 0
  total$substitution <- SUBSTITUTIONS_96
  context = CONTEXTS_96
  colors = COLORS6
  ggplot(data = total, aes(x = context, y = fraction, fill = substitution, width = 1, ymax = top, ymin = bottom)) + 
    geom_bar(stat = "identity",colour = "black", size = 0.4) + 
    scale_fill_manual(values = colors) + 
    facet_grid( ~ substitution) + 
    ylab("Relative contribution") + 
    coord_cartesian(ylim = c(0, ymax)) + 
    scale_y_continuous(breaks = seq(0, ymax, 0.02), expand = c(0, 0)) + 
    guides(fill = FALSE) + theme_classic() + 
    expand_limits(y = 0) +  
    theme(axis.title.y = element_text(size = 12, vjust = 1), 
          axis.text.y = element_text(size = 8), axis.title.x = element_text(size = 12), 
          axis.text.x = element_text(size = 5, angle = 90, 
                                     vjust = 0.4), strip.text.x = element_text(size = 9), 
          strip.text.y = element_text(size = 9), panel.grid.major.x = element_blank(), 
          panel.spacing.x = unit(0, "lines")) + 
    geom_errorbar(color = "black", width = 0.5, size = 0.3)
}
