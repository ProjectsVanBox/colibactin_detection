plot_compare_profiles2 = function (profile1,
                                   profile2,
                                   profile_names = c("profile 1",
                                                     "profile 2"),
                                   profile_ymax = 0.2,
                                   diff_ylim = c(-0.02, 0.02),
                                   colors = NA,
                                   condensed = FALSE)
{
  value <-
    substitution <-
    Sample <- Contribution <- Signature <- variable <- NULL
  full_context <- context <- NULL
  if (is.na(colors)) {
    colors <- COLORS6
  }
  comp <-
    .create_profile_comparison(profile1, profile2, profile_names)
  df <-
    comp$matrix %>% as.data.frame() %>% tibble::rownames_to_column("full_context") %>%
    dplyr::mutate(
      substitution = stringr::str_replace(full_context,
                                          "\\w\\[(.*)\\]\\w", "\\1"),
      context = stringr::str_replace(full_context,
                                     "\\[.*\\]", "\\.")
    ) %>% dplyr::select(-full_context) %>%
    tidyr::pivot_longer(c(-substitution,-context),
                        names_to = "sample",
                        values_to = "value") %>% dplyr::mutate(sample = factor(sample,
                                                                               levels = unique(sample)))
  df_blank <-
    .create_dummy_limits(df[, c("substitution", "context")],
                         profile_names, profile_ymax, diff_ylim)
  if (condensed == TRUE) {
    width <- 1
    spacing <- 0
  }
  else {
    width <- 0.6
    spacing <- 0.5
  }
  plot <-
    ggplot(data = df,
           aes(
             x = context,
             y = value,
             fill = substitution,
             width = width
           )) + geom_bar(stat = "identity",
                         position = "identity",
                         size = 0.2)
  + geom_blank(data = df_blank,
               aes(x = context, y = value)) + scale_fill_manual(values = colors) +
    facet_grid(sample ~ substitution, scales = "free_y") +
    labs(y = "Relative contribution", title = comp$title) +
    guides(fill = FALSE) + theme_bw() + theme(
      axis.title.y = element_text(size = 12,
                                  vjust = 1),
      axis.text.y = element_text(size = 8),
      axis.title.x = element_text(size = 12),
      axis.text.x = element_text(
        size = 5,
        angle = 90,
        vjust = 0.5
      ),
      strip.text.x = element_text(size = 14),
      strip.text.y = element_text(size = 14),
      panel.grid.major.x = element_blank(),
      panel.spacing.x = unit(spacing,
                             "lines")
    )
  return(plot)
}
