# plot_contribution2 
# Axel Rosendahl Huber 5 february 2018
# Aim: Modify plot_contribution so that the y axis crosses x axis at zero

plot_contribution2 <- function (contribution, signatures, index = c(), coord_flip = FALSE, 
          mode = "relative", palette = c()) 
{
  if (!(mode == "relative" | mode == "absolute")) 
    stop("mode parameter should be either 'relative' or 'absolute'")
  if (length(index > 0)) {
    contribution = contribution[, index]
  }
  Sample = NULL
  Contribution = NULL
  Signature = NULL
  if (mode == "relative") {
    m_contribution = melt(contribution)
    colnames(m_contribution) = c("Signature", "Sample", "Contribution")
    plot = ggplot(m_contribution, aes(x = factor(Sample), 
                                      y = Contribution, fill = factor(Signature), order = Sample)) + 
      geom_bar(position = "fill", stat = "identity", colour = "black") + 
      labs(x = "", y = "Relative contribution") + theme_BM() + 
      expand_limits(y = 0) +  scale_y_continuous(expand = c(0, 0)) +
      theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank()) + 
      theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank())
  }
  else {
    if (missing(signatures)) 
      stop(paste("For contribution plotting in mode 'absolute':", 
                 "also provide signatures matrix"))
    total_signatures = colSums(signatures)
    abs_contribution = contribution * total_signatures
    m_contribution = melt(abs_contribution)
    colnames(m_contribution) = c("Signature", "Sample", "Contribution")
    plot = ggplot(m_contribution, aes(x = factor(Sample),y = Contribution, fill = factor(Signature), order = Sample)) + 
      geom_bar(stat = "identity", colour = "black") +
      labs(x = "",y = "Absolute contribution (no. mutations)") + 
      expand_limits(y = 0) +  scale_y_continuous(expand = c(0, 0)) +
      theme_BM() + 
      theme(panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank()) +
      theme(panel.grid.minor.y = element_blank(),panel.grid.major.y = element_blank())
  }
  if (length(palette) > 0) 
    plot = plot + scale_fill_manual(name = "Signature", values = palette)
  else plot = plot + scale_fill_discrete(name = "Signature")
  if (coord_flip) 
    plot = plot + coord_flip() + xlim(rev(levels(factor(m_contribution$Sample))))
  else plot = plot + xlim(levels(factor(m_contribution$Sample)))
  return(plot)
}
