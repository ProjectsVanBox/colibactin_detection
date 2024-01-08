# plot_compare_profiles192
plot_compare_profiles192 <- function (profile1, profile2, profile_names = c("profile 1","profile 2"), 
                                      profile_ymax = 0.05, diff_ylim = c(-0.02, 0.02), condensed = TRUE) { 
  colors = COLORS6
  s1_relative = profile1/sum(profile1)
  s2_relative = profile2/sum(profile2)
  diff = s1_relative - s2_relative
  RSS = sum(diff^2)
  RSS = format(RSS, scientific = TRUE, digits = 3)
  cosine_sim = cos_sim(profile1, profile2)
  cosine_sim = round(cosine_sim, 3)
  x = cbind(s1_relative, s2_relative, diff)
  colnames(x) = c(profile_names, "Difference")
  substitutions = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  index = c(rep(1, 1, 32), rep(2, 1, 32), rep(3, 1, 32), rep(4, 
                                                             1, 32), rep(5, 1, 32), rep(6, 1, 32))
  load("~/surfdrive/Shared/vanBoxtelLab (Groupfolder)/Boxtel_General/Scripts/pmc_vanboxtel/Axel_BMseq/data/CONTEXTS_192")
  context = CONTEXTS_192
  context <- gsub(x = context,pattern =  "untranscribed", "U")
  context <- gsub(x = context,pattern =  "transcribed", "T")
  context <- gsub(x = context, pattern = "\\[.*\\]", ".")
  df = data.frame(substitution = substitutions[index], context = context)
  rownames(x) = NULL
  df2 = cbind(df, as.data.frame(x))
  df3 = melt(df2, id.vars = c("substitution", "context"))
  df3$substitution <- as.factor(df3$substitution)
  df3$strand <- "T"
  df3$strand[grep(df3$context, pattern = "U")] <- "U"
  df3$strand <- factor(df3$strand, levels = c("U", "T"))
  value = NULL
  substitution = NULL
  Sample = NULL
  Contribution = NULL
  Signature = NULL
  df4 = data.frame(substitution = rep("C>A", 4), context = rep("A.A", 4), variable = c(profile_names, "Difference", "Difference"), 
                   value = c(profile_ymax, profile_ymax, diff_ylim[1], diff_ylim[2]))
  
  plot = ggplot(data = df3, aes(x = context, y = value, fill = substitution, alpha = strand)) +
    geom_bar(stat = "identity", colour = "black", size = 0.2) +
    scale_fill_manual(values = colors) + facet_grid(variable ~ substitution , scales = "free") + ylab("Relative contribution") + 
    guides(fill = FALSE) + theme_BM() + 
    theme(axis.title.y = element_text(size = 12, vjust = 1), axis.text.y = element_text(size = 8), 
    axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 5, angle = 90, vjust = 0.4), strip.text.x = element_text(size = 14), 
    strip.text.y = element_text(size = 14), panel.grid.major.x = element_blank(), panel.spacing.x = unit(0, "lines")) +
    ggtitle(paste0("Cosine similarity = ", cosine_sim)) 
  return(plot)
}

