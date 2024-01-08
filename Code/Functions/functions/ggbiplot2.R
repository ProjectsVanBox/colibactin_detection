# ggbiplot2 
# Axel Rosendahl Huber 
# 30-4-2018
# modification of ggbiplot: https://github.com/vqv/ggbiplot
# Enables selection of the most significant contributors to the 

ggbiplot2 <- function(pcobj, choices = 1:2, scale = 1, pc.biplot = TRUE, 
                       obs.scale = 1 - scale, var.scale = scale, groups = NULL, 
                       ellipse = FALSE, ellipse.prob = 0.68, labels = NULL, labels.size = 3, 
                       alpha = 1, var.axes = TRUE, circle = FALSE, circle.prob = 0.69, 
                       varname.size = 3, varname.adjust = 1.5, varname.abbrev = FALSE, selected_vectors = NULL,
                       repel_text, selected_labels = NULL, point_size = 2,
                       ...) {
  library(ggplot2)
  library(plyr)
  library(scales)
  library(grid)
  library(ggrepel)
  
  stopifnot(length(choices) == 2)
  if (inherits(pcobj, "prcomp")) {
    nobs.factor <- sqrt(nrow(pcobj$x) - 1)
    d <- pcobj$sdev
    u <- sweep(pcobj$x, 2, 1/(d * nobs.factor), FUN = "*")
    v <- pcobj$rotation
  }
  else if (inherits(pcobj, "princomp")) {
    nobs.factor <- sqrt(pcobj$n.obs)
    d <- pcobj$sdev
    u <- sweep(pcobj$scores, 2, 1/(d * nobs.factor), FUN = "*")
    v <- pcobj$loadings
  }
  else if (inherits(pcobj, "PCA")) {
    nobs.factor <- sqrt(nrow(pcobj$call$X))
    d <- unlist(sqrt(pcobj$eig)[1])
    u <- sweep(pcobj$ind$coord, 2, 1/(d * nobs.factor), FUN = "*")
    v <- sweep(pcobj$var$coord, 2, sqrt(pcobj$eig[1:ncol(pcobj$var$coord), 
                                                  1]), FUN = "/")
  }
  else if (inherits(pcobj, "lda")) {
    nobs.factor <- sqrt(pcobj$N)
    d <- pcobj$svd
    u <- predict(pcobj)$x/nobs.factor
    v <- pcobj$scaling
    d.total <- sum(d^2)
  }
  else {
    stop("Expected a object of class prcomp, princomp, PCA, or lda")
  }
  choices <- pmin(choices, ncol(u))
  df.u <- as.data.frame(sweep(u[, choices], 2, d[choices]^obs.scale, 
                              FUN = "*"))
  v <- sweep(v, 2, d^var.scale, FUN = "*")
  df.v <- as.data.frame(v[, choices])
  names(df.u) <- c("xvar", "yvar")
  names(df.v) <- names(df.u)
  if (pc.biplot) {
    df.u <- df.u * nobs.factor
  }
  r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(df.u^2))^(1/4)
  v.scale <- rowSums(v^2)
  df.v <- r * df.v/sqrt(max(v.scale))
  if (obs.scale == 0) {
    u.axis.labs <- paste("standardized PC", choices, sep = "")
  }
  else {
    u.axis.labs <- paste("PC", choices, sep = "")
  }
  u.axis.labs <- paste(u.axis.labs, sprintf("(%0.1f%% explained var.)", 
                                            100 * pcobj$sdev[choices]^2/sum(pcobj$sdev^2)))
  if (!is.null(labels)) {
    df.u$labels <- labels
  }
  if (!is.null(groups)) {
    df.u$groups <- groups
  }
  if (varname.abbrev) {
    df.v$varname <- abbreviate(rownames(v))
  }
  else {
    df.v$varname <- rownames(v)
  }
  
  if (!is.null(selected_vectors)) {
    size_vdist <- sqrt(df.v$xvar^2 + df.v$yvar^2) 
    names(size_vdist) <- TRIPLETS_96
    sort_vdist <- sort(size_vdist, decreasing = T)
    select_vdist <- names(sort_vdist[1:selected_vectors])
    df.v <- df.v[select_vdist,]
  }
  
  
  
  df.v$angle <- with(df.v, (180/pi) * atan(yvar/xvar))
  df.v$hjust = with(df.v, (1 - varname.adjust * sign(xvar))/2)
  g <- ggplot(data = df.u, aes(x = xvar, y = yvar)) + xlab(u.axis.labs[1]) + 
    ylab(u.axis.labs[2]) + coord_equal()
  if (var.axes) {
    if (circle) {
      theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, 
                                                length = 50))
      circle <- data.frame(xvar = r * cos(theta), yvar = r * 
                             sin(theta))
      g <- g + geom_path(data = circle, color = muted("white"), 
                         size = 1/2, alpha = 1/3)
    }
    g <- g + geom_segment(data = df.v, aes(x = 0, y = 0, 
                                           xend = xvar, yend = yvar), arrow = arrow(length = unit(1/2, 
                                                                                                  "picas")), color = muted("red"))
  }
  if (!is.null(df.u$labels)) {
    if (!is.null(df.u$groups)) {
      g <- g + geom_text(aes(label = labels, color = groups), 
                         size = labels.size)
    }
    else {
      g <- g + geom_text(aes(label = labels), size = labels.size)
    }
  }
  else {
    if (!is.null(df.u$groups)) {
      g <- g + geom_point(aes(fill = groups), color = "black",alpha = alpha, shape = 21, size = point_size)
    }
    else {
      g <- g + geom_point(alpha = alpha)
    }
  }
  if (!is.null(df.u$groups) && ellipse) {
    theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
    circle <- cbind(cos(theta), sin(theta))
    ell <- ddply(df.u, "groups", function(x) {
      if (nrow(x) <= 2) {
        return(NULL)
      }
      sigma <- var(cbind(x$xvar, x$yvar))
      mu <- c(mean(x$xvar), mean(x$yvar))
      ed <- sqrt(qchisq(ellipse.prob, df = 2))
      data.frame(sweep(circle %*% chol(sigma) * ed, 2, 
                       mu, FUN = "+"), groups = x$groups[1])
    })
    names(ell)[1:2] <- c("xvar", "yvar")
    g <- g + geom_path(data = ell, aes(color = groups, group = groups))
  }
  if (var.axes) {
    g <- g + geom_text(data = df.v, aes(label = varname, 
                                        x = xvar, y = yvar, angle = angle, hjust = hjust), 
                       color = "darkred", size = varname.size)
  }
  if (repel_text) {
    if (is.null(selected_labels)){
      g <- g + geom_text_repel(data = df.u, aes(label = groups), 
                               box.padding = unit(1.75, 'lines'))
    } else {
      index <- groups %in% selected_labels
      g <- g + geom_text_repel(data = df.u[index,], 
                               aes(label = as.character(df.u$groups[index])),
                               box.padding = unit(1.75, 'lines'))
      }
   
  }
  return(g)
}
