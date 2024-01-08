plot_indel_contexts2 = function (counts, same_y = FALSE, extra_labels = FALSE, condensed = FALSE) 
{
  
  
  INDEL_COLORS <- c(
    "#FDBE6F", "#FF8001", "#B0DD8B", "#36A12E", "#FDCAB5", "#FC8A6A",
    "#F14432", "#BC141A", "#D0E1F2", "#94C4DF", "#4A98C9", "#1764AB",
    "#E2E2EF", "#B6B6D8", "#8683BD", "#61409B"
  )
  
  
  count <- muttype <- muttype_sub <- muttype_total <- sample <- NULL
  counts <- counts %>% as.data.frame() %>% tibble::rownames_to_column("muttype_total") %>% 
    tidyr::separate(muttype_total, c("muttype", "muttype_sub"), 
                    sep = "_(?=[0-9])") %>% dplyr::mutate(muttype = factor(muttype, 
                                                                           levels = unique(muttype))) %>% tidyr::gather(key = "sample", 
                                                                                                                        value = "count", -muttype, -muttype_sub) %>% dplyr::mutate(sample = factor(sample, 
                                                                                                                                                                                                   levels = unique(sample)))
  nr_muts <- counts %>% dplyr::group_by(sample) %>% dplyr::summarise(nr_muts = round(sum(count)))
  facet_labs_y <- nr_muts$sample
  names(facet_labs_y) <- nr_muts$sample
  facet_labs_x <- c("1: C", "1: T", "1: C", "1: T", 2, 3, 4, 
                    "5+", 2, 3, 4, "5+", 2, 3, 4, "5+")
  names(facet_labs_x) <- levels(counts$muttype)
  if (same_y) {
    facet_scale <- "free_x"
  }
  else {
    facet_scale <- "free"
  }
  if (extra_labels) {
    title <- stringr::str_c("Deletion           ", "Insertion          ", 
                            "Deletion                                   ", "Insertion                                  ", 
                            "Deletion (MH)")
    x_lab <- stringr::str_c("Homopolymer length                            ", 
                            "Number of repeat units                                                                               ", 
                            "Microhomology length")
  }
  else {
    title <- x_lab <- ""
  }
  if (condensed == TRUE) {
    width <- 1
    spacing <- 0
  }
  else {
    width <- 0.6
    spacing <- 0.5
  }
  fig <- ggplot(counts, aes(x = muttype_sub, y = count, fill = muttype, 
                            width = width)) + geom_bar(stat = "identity") + facet_grid(sample ~ 
                                                                                         muttype, scales = facet_scale, space = "free_x", labeller = labeller(muttype = facet_labs_x, 
                                                                                                                                                              sample = facet_labs_y)) + scale_fill_manual(values = INDEL_COLORS) + 
    theme_bw() + labs(fill = "Mutation type", title = title, 
                      y = "Nr of indels", x = x_lab) + theme(panel.grid.major.x = element_blank(), 
                                                             panel.grid.minor.y = element_blank(), panel.spacing.x = unit(spacing, 
                                                                                                                          "lines"))
  return(fig)
}