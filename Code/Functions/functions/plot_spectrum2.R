# plot_spectrum2
# 
# Axel Rosendahl Huber 2-8-2018
# Modified version of the plot_spectrum function in the MutationalPatterns package.
# Additions: Ability to add individual data points to the plot
# Aestetic changes: Black outline added to the plot. (theme_BM() = TRUE)

plot_spectrum2 <- function (type_occurrences, CT_at_CpG = TRUE, by, colors, legend = TRUE, 
          dot_plot = TRUE ,theme_BM = TRUE) 
{
  value = NULL
  nmuts = NULL
  sub_type = NULL
  variable = NULL
  error_pos = NULL
  stdev = NULL
  
  if (missing(colors)) 
    colors = COLORS7
  if (length(colors) != 7) 
    stop("Colors parameter: supply color vector with length 7")
  if (CT_at_CpG == FALSE) 
    type_occurrences2 = type_occurrences[, 1:6]
  else type_occurrences2 = type_occurrences[, c(1:2, 8, 7, 4:6)]
  df2 = type_occurrences2/rowSums(type_occurrences2)
  if (missing(by)) 
    by = "all"
  df2$by = by
  df3 = melt(df2, id.vars = "by")
  df3$subtype <- substring(df3$variable,1,3)
  counts = melt(type_occurrences2, measure.vars = colnames(type_occurrences2))
  df4 = cbind(df3, counts$value)
  colnames(df4)[5] = "nmuts"
  x = ddply(df4, c("by", "variable"), summarise, mean = mean(value), 
            stdev = sd(value))
  info_x = ddply(df4, c("by"), summarise, total_individuals = sum(value), 
                 total_mutations = sum(nmuts))
  x = merge(x, info_x)
  info_type = data.frame(sub_type = c("C>A", "C>G", "C>T","C>T", "C>T", "T>A", "T>C", "T>G"), 
                         variable = c("C>A","C>G", "C>T", "C>T at CpG", "C>T other", "T>A", "T>C","T>G"))
  x = merge(x, info_type)
  x$total_mutations = prettyNum(x$total_mutations, big.mark = ",")
  x$total_mutations = paste("n =", as.character(x$total_mutations), "substitutions")
  x$error_pos = x$mean
  if (CT_at_CpG == FALSE) 
    colors = colors[c(1, 2, 3, 5:7)]
  else {
    x = x[order(x$by), ]
    CpG = which(x$variable == "C>T at CpG")
    other = which(x$variable == "C>T other")
    x$error_pos[CpG] = x$error_pos[other] + x$error_pos[CpG]
    order = order(factor(x$variable, 
                         levels = c("C>A", "C>G", "C>T other", "C>T at CpG", "T>A", "T>C", "T>G")))
    x = x[order, ]
    
  }
  plot = ggplot(data = x, aes(x = sub_type, y = mean, fill = variable, group = sub_type)) + 
    geom_bar(stat = "identity", color = "black", alpha = 0.7) + 
    scale_fill_manual(values = colors,name = "Point mutation type") + 
    xlab("") + theme_BM() + 
    ylab("Relative contribution") + 
    theme(axis.text.x = element_blank(), 
          panel.grid.major.x = element_blank()) +
    scale_y_continuous(expand = c(0, 0))
  
    if (sum(is.na(x$stdev)) > 0) 
    warning("No standard deviation error bars can be plotted, because there is only one sample per mutation spectrum")
  else plot = plot + geom_errorbar(aes(ymin = error_pos - stdev, 
                                       ymax = error_pos + stdev), width = 0.2)
  if (dot_plot) {
    df_dotplot <- df3
    if (CT_at_CpG) {
      for (i in row.names(df_dotplot[df_dotplot$variable == "C>T at CpG",])) {
        i = as.numeric(i)
        new_value <- df_dotplot$value[i] + x$mean[x$variable == "C>T other" & x$by == df_dotplot$by[i]]
        df_dotplot$value[i] <- new_value
      }
    }
    plot <- plot + geom_dotplot(data =  df_dotplot, aes(x = subtype, y = value, group = variable), fill= "black", binaxis = "y", stackdir = "center", binwidth = 0.01, dotsize =  .5)
  }
  
  if (length(by) == 1) 
    plot = plot + facet_wrap(~total_mutations)
  else plot = plot + facet_wrap(by ~ total_mutations)
  if (legend == FALSE) 
    plot = plot + theme(legend.position = "none")
  return(plot)
}

