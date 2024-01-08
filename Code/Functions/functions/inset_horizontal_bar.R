inset_horizontal_bar = function (tree_view, insets, width, height, hjust = 0, vjust = 0, 
          x = "node", reverse_x = FALSE, reverse_y = FALSE) 
{
#  if (width < 0 || width > 1) 
#    stop("width should be in range of (0,1)")
#  if (height < 0 || height > 1) 
#    stop("height should be in range of (0,1)")
  df <- tree_view$data[as.numeric(names(insets)), ]
  x <- match.arg(x, c("node", "branch", "edge"))
  if (x == "node") {
    xx <- df$x
  }
  else {
    xx <- df$branch
  }
  yy <- df$y
  xx <- xx - hjust
  yy <- yy - vjust
  if (reverse_x) 
    xx <- -xx
  if (reverse_y) 
    yy <- -yy
  #width <- width * diff(range(tree_view$data$x, na.rm = TRUE))
  height <- height * diff(range(tree_view$data$y, na.rm = TRUE))
  subview_barplot =   tree_view + ggimage::geom_subview(subview = insets, width = width, 
                           height = height, x = xx, y = yy) + coord_flip()
  return(subview_barplot)
}

