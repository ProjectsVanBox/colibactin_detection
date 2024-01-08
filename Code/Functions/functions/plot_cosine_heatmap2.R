plot_cosine_heatmap2 = function (cos_sim_matrix, col_order = NA, row_order = NA, cluster_rows = TRUE, 
                                 cluster_cols = FALSE, method = "complete", plot_values = FALSE) 
{
  if (!inherits(cos_sim_matrix, "matrix")) {
    stop("cos_sim_matrix must be a matrix")
  }
  if (length(colnames(cos_sim_matrix)) == 0) {
    stop("cos_sim_matrix is missing colnames")
  }
  if (length(rownames(cos_sim_matrix)) == 0) {
    stop("cos_sim_matrix is missing rownames")
  }
  Cosine.sim <- Signature <- Sample <- x <- y <- xend <- yend <- NULL
  if (!is.na(row_order) & cluster_rows == TRUE) {
    stop("row_order can only be provided when cluster_rows is FALSE", 
         call. = FALSE)
  }
  else if (!is.na(row_order)) {
    if (!inherits(row_order, "character")) {
      stop("row_order must be a character vector", call. = FALSE)
    }
    if (length(row_order) != nrow(cos_sim_matrix)) {
      stop("row_order must have the same length as the number of\n          samples in the explained matrix", 
           call. = FALSE)
    }
  }
  else if (cluster_rows == TRUE) {
    hc.sample <- hclust(dist(cos_sim_matrix), method = method)
    row_order <- rownames(cos_sim_matrix)[hc.sample$order]
    dhc <- as.dendrogram(hc.sample)
    ddata <- ggdendro::dendro_data(dhc, type = "rectangle")
    dendrogram_rows <- ggplot(ggdendro::segment(ddata)) + 
      geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
      coord_flip() + scale_y_reverse(expand = c(0.2, 0)) + 
      ggdendro::theme_dendro()
  }
  else {
    row_order <- rownames(cos_sim_matrix)
  }
  if (!is.na(col_order) & cluster_cols == TRUE) {
    stop("col_order can only be provided when cluster_cols is FALSE", 
         call. = FALSE)
  }
  else if (!is.na(col_order)) {
    if (!inherits(col_order, "character")) {
      stop("col_order must be a character vector", call. = FALSE)
    }
    if (length(col_order) != ncol(cos_sim_matrix)) {
      stop("col_order must have the same length as the number of \n          signatures in the explained matrix", 
           call. = FALSE)
    }
  }
  else if (cluster_cols == TRUE) {
    hc.sample2 <- cos_sim_matrix %>% t() %>% dist() %>% hclust(method = method)
    col_order <- colnames(cos_sim_matrix)[hc.sample2$order]
    dhc <- as.dendrogram(hc.sample2)
    ddata <- ggdendro::dendro_data(dhc, type = "rectangle")
    dendrogram_cols <- ggplot(ggdendro::segment(ddata)) + 
      geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
      ggdendro::theme_dendro() + scale_y_continuous(expand = c(0.2, 
                                                               0))
  }
  else {
    col_order <- colnames(cos_sim_matrix)
  }
  cos_sim_matrix.m <- cos_sim_matrix %>% as.data.frame() %>% 
    tibble::rownames_to_column("Sample") %>% tidyr::pivot_longer(-Sample, 
                                                                 names_to = "Signature", values_to = "Cosine.sim") %>% 
    dplyr::mutate(Signature = factor(Signature, levels = col_order), 
                  Sample = factor(Sample, levels = row_order))
  heatmap <- ggplot(cos_sim_matrix.m, aes(x = Signature, y = Sample, 
                                          fill = Cosine.sim, order = Sample)) + geom_raster() + 
    scale_fill_distiller(palette = "YlGnBu", direction = 1, 
                         name = "Cosine \nsimilarity", limits = c(0, 1.000000001)) + 
    theme_bw() + theme(axis.text.x = element_text(angle = 90, 
                                                  hjust = 1, vjust = 0.5, size = 8), panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank()) + labs(x = NULL, 
                                                                  y = NULL)
  if (plot_values == TRUE) {
    heatmap <- heatmap + geom_text(aes(label = round(Cosine.sim, 
                                                     2)), size = 2)
  }
  if (cluster_rows == TRUE & cluster_cols == TRUE) {
    plot_final <- cowplot::plot_grid(dendrogram_rows, heatmap, 
                                     align = "h", rel_widths = c(0.3, 1))
  }
  else if (cluster_rows == TRUE & cluster_cols == FALSE) {
    plot_final <- cowplot::plot_grid(dendrogram_rows, heatmap, 
                                     align = "h", rel_widths = c(0.3, 1))
  }
  else if (cluster_rows == FALSE & cluster_cols == TRUE) {
    plot_final <- cowplot::plot_grid(dendrogram_cols, heatmap, 
                                     align = "v", rel_heights = c(0.3, 1)) + ylim(rev(levels(factor(cos_sim_matrix.m$Sample))))
  }
  else {
    plot_final <- heatmap + ylim(rev(levels(factor(cos_sim_matrix.m$Sample))))
  }
  return(plot_final)
}
