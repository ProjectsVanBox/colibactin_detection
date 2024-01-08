# Function plot_strand2
# Axel Rosendahl Huber 23-01-2018
# Modified function of MutationalPatterns function plot_strand 
# Modifications: No floating bargraph, theme_classic, 

plot_strand2 <- function (strand_bias_df, mode = "relative", colors, ymax = 0.3) 
{
  if (missing(colors)) 
    colors = COLORS6
  type = NULL
  relative_contribution = NULL
  no_mutations = NULL
  if (mode == "relative") {
    plot = ggplot(strand_bias_df, aes(x = type, y = relative_contribution, 
            fill = type, alpha = strand)) + geom_bar(stat = "identity", 
            position = "dodge", colour = "black", cex = 0.5) + 
            scale_fill_manual(values = colors) + scale_alpha_discrete(range = c(1, 
            0.4)) + ylab("Relative contribution") + facet_grid(. ~ group) +
            theme_classic() + 
            expand_limits(y = 0) +  scale_y_continuous(expand = c(0, 0), limits = c(0, ymax), breaks = seq(0,ymax, 0.1)) + 
            scale_x_discrete(breaks = NULL) + xlab("")
  }
  else if (mode == "absolute") {
    plot = ggplot(strand_bias_df, aes(x = type, y = no_mutations, 
              fill = type, alpha = strand)) + geom_bar(stat = "identity", 
              position = "dodge", colour = "black", cex = 0.5) + 
              scale_fill_manual(values = colors) + 
              scale_alpha_discrete(range = c(1, 0.4)) + 
              expand_limits(y = 0) +  scale_y_continuous(expand = c(0, 0)) + 
              ylab("Total number of mutations") + facet_grid(. ~ group) + 
              theme_bw() + scale_x_discrete(breaks = NULL) + 
              xlab("")
  }
  return(plot)
}
