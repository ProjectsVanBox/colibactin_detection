table_spectrum <- function (type_occurrences, CT = FALSE, by, colors, legend = TRUE) 
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
  if (CT == FALSE) 
    type_occurrences = type_occurrences[, 1:6]
  else type_occurrences = type_occurrences[, c(1:2, 8, 7, 4:6)]
  df2 = type_occurrences/rowSums(type_occurrences)
  if (missing(by)) 
    by = "all"
  df2$by = by
  df3 = melt(df2, id.vars = "by")
  counts = melt(type_occurrences, measure.vars = colnames(type_occurrences))
  df4 = cbind(df3, counts$value)
  colnames(df4)[4] = "nmuts"
  x = ddply(df4, c("by", "variable"), summarise, mean = mean(value), 
            stdev = sd(value))
  info_x = ddply(df4, c("by"), summarise, total_individuals = sum(value), 
                 total_mutations = sum(nmuts))
  x = merge(x, info_x)
  info_type = data.frame(sub_type = c("C>A", "C>G", "C>T", 
                                      "C>T", "C>T", "T>A", "T>C", "T>G"), variable = c("C>A", 
                                                                                       "C>G", "C>T", "C>T at CpG", "C>T other", "T>A", "T>C", 
                                                                                       "T>G"))
  x = merge(x, info_type)
  x$total_mutations = prettyNum(x$total_mutations, big.mark = ",")
  x$total_mutations = paste("No. mutations =", as.character(x$total_mutations))
  x$error_pos = x$mean
  if (CT == FALSE) 
    colors = colors[c(1, 2, 3, 5:7)]
  else {
    x = x[order(x$by), ]
    CpG = which(x$variable == "C>T at CpG")
    other = which(x$variable == "C>T other")
    x$error_pos[CpG] = x$error_pos[other] + x$error_pos[CpG]
    order = order(factor(x$variable, levels = c("C>A", "C>G", 
                                                "C>T other", "C>T at CpG", "T>A", "T>C", "T>G")))
    x = x[order, ]
  }
  plot = ggplot(data = x, aes(x = sub_type, y = mean, fill = variable, 
                              group = sub_type)) + geom_bar(stat = "identity") + scale_fill_manual(values = colors, 
                                                                                                   name = "Point mutation type") + theme_bw() + xlab("") + 
    ylab("Relative contribution") + theme(axis.ticks = element_blank(), 
                                          axis.text.x = element_blank(), panel.grid.major.x = element_blank())
  if (sum(is.na(x$stdev)) > 0) 
    warning("No standard deviation error bars can be plotted, because there is only one sample per mutation spectrum")
  else plot = plot + geom_errorbar(aes(ymin = error_pos - stdev, 
                                       ymax = error_pos + stdev), width = 0.2)
  if (length(by) == 1) 
    plot = plot + facet_wrap(~total_mutations)
  else plot = plot + facet_wrap(by ~ total_mutations)
  if (legend == FALSE) 
    plot = plot + theme(legend.position = "none")
  return(x)
}
