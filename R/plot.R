plot_diag_ <- function(path, dim, t.lim) {
  p <- ggplot()

  jd <- c(1, which(path$W_t[,dim] != path$W_tm[,dim]), length(path$t))
  for(j in 2:length(jd)) {
    if(j != length(jd))
      p <- p + geom_vline(xintercept = path$t[jd[j]], linetype = "dotted")
    if(jd[j]-jd[j-1] > 0)
      p <- p + geom_line(aes(x = .data$t, y = .data$x), tibble(t = path$t[jd[j-1]:jd[j]], x = c(path$W_t[jd[j-1]:(jd[j]-1),dim], path$W_tm[jd[j],dim])), colour = "grey")
  }

  p <- p + geom_point(aes(x = .data$t, y = .data$x), tibble(t = path$t, x = path$W_t[,dim]), colour = "black", size = 0.1)
  jd <- jd[c(-1, -length(jd))]
  if(length(jd) > 0) {
    p <- p + geom_point(aes(x = .data$t, y = .data$x), tibble(t = path$t[jd], x = path$W_tm[jd,dim]), colour = "black", size = 1, shape = 1)
  }

  lyrs <- data.frame(t.l = path$layers$t.l,
                     t.u = path$layers$t.u,
                     L = path$layers$outer.L[,dim],
                     U = path$layers$outer.U[,dim])

  if(nrow(path$layers) > 0) {
    p <- p +
      geom_segment(aes(x = .data$t.l, xend = .data$t.u, y = .data$L, yend = .data$L), lyrs, colour = "black") +
      geom_segment(aes(x = .data$t.l, xend = .data$t.u, y = .data$U, yend = .data$U), lyrs, colour = "black")
  }

  p + coord_cartesian(xlim = t.lim) + xlab("t") + ylab(paste0("W(",dim,")")) #, ylim = c(ymin,ymax))
}

plot_offdiag_ <- function(path, dim, t.lim, layers.2d, show.layers) {
  p <- ggplot()

  lyrs <- data.frame(t.l = path$layers$t.l,
                     t.u = path$layers$t.u,
                     L1 = path$layers$outer.L[,dim[1]],
                     U1 = path$layers$outer.U[,dim[1]],
                     L2 = path$layers$outer.L[,dim[2]],
                     U2 = path$layers$outer.U[,dim[2]])

  if(show.layers && nrow(path$layers) > 0) {
    if(layers.2d == "detailed") {
      for(i in 1:nrow(path$layers)) {
        cube <- as.data.frame(path$layers$inner.cube[[i]][,dim])
        p <- p +
          geom_polygon(aes(x = .data$V2, y = .data$V1), cube[chull(cube),], alpha = 0, colour = "red")
      }

      p <- p +
        geom_segment(aes(x = .data$L2, xend = .data$L2, y = .data$L1, yend = .data$U1), lyrs, colour = "black") +
        geom_segment(aes(x = .data$U2, xend = .data$U2, y = .data$L1, yend = .data$U1), lyrs, colour = "black") +
        geom_segment(aes(x = .data$L2, xend = .data$U2, y = .data$L1, yend = .data$L1), lyrs, colour = "black") +
        geom_segment(aes(x = .data$L2, xend = .data$U2, y = .data$U1, yend = .data$U1), lyrs, colour = "black")
    } else if(layers.2d == "union") {
      if(!("sf" %in% installed.packages()[,"Package"])) {
        stop("The union layering option requires the suggested package sf. NB this option is more computationally intensive than default option.")
      }

      union.parts <- sort(unique(c(min(path$layers$t.l), max(path$layers$t.u), path$labels$seg.end)))
     # union.parts <- sort(unique(c(t.lim, path$labels$seg.end)))
      union.parts <- rbind(union.parts[1:(length(union.parts)-1)],union.parts[2:length(union.parts)])

      for(i in 1:ncol(union.parts)) {
        poly.lyrs <- which(path$layers$t.l >= union.parts[1,i] & path$layers$t.u <= union.parts[2,i])
        poly <- lapply(path$layers[poly.lyrs,]$inner.cube, function(poly) {
          poly <- poly[chull(poly[,dim[2:1]]),dim[2:1]]
          poly <- rbind(poly, poly[1,])
          sf::st_polygon(list(poly))
        })

        cube <- sf::st_union(sf::st_buffer(sf::st_sfc(poly), 0.000001));
        # p <- p + geom_sf(data=cube, fill=NA, colour="red")
        cube2 <- as.data.frame(cube[[1]][[1]]); names(cube2) <- c("V1", "V2")
        p <- p + geom_polygon(aes(x = .data$V1, y = .data$V2), cube2, alpha = 0, colour = "red")
      }

      lyrs <- data.frame(L1 = min(lyrs$L1),
                         L2 = min(lyrs$L2),
                         U1 = max(lyrs$U1),
                         U2 = max(lyrs$U2))

      p <- p +
        geom_segment(aes(x = .data$L2, xend = .data$L2, y = .data$L1, yend = .data$U1), lyrs, colour = "black") +
        geom_segment(aes(x = .data$U2, xend = .data$U2, y = .data$L1, yend = .data$U1), lyrs, colour = "black") +
        geom_segment(aes(x = .data$L2, xend = .data$U2, y = .data$L1, yend = .data$L1), lyrs, colour = "black") +
        geom_segment(aes(x = .data$L2, xend = .data$U2, y = .data$U1, yend = .data$U1), lyrs, colour = "black")
    } else if(layers.2d == "convex") {
      cube <- as.data.frame(do.call("rbind", lapply(path$layers$inner.cube, function(x) { x[,dim] })))
      p <- p +
          geom_polygon(aes(x = .data$V2, y = .data$V1), cube[chull(cube),], alpha = 0, colour = "red")

      lyrs <- data.frame(L1 = min(lyrs$L1),
                         L2 = min(lyrs$L2),
                         U1 = max(lyrs$U1),
                         U2 = max(lyrs$U2))

      p <- p +
        geom_segment(aes(x = .data$L2, xend = .data$L2, y = .data$L1, yend = .data$U1), lyrs, colour = "black") +
        geom_segment(aes(x = .data$U2, xend = .data$U2, y = .data$L1, yend = .data$U1), lyrs, colour = "black") +
        geom_segment(aes(x = .data$L2, xend = .data$U2, y = .data$L1, yend = .data$L1), lyrs, colour = "black") +
        geom_segment(aes(x = .data$L2, xend = .data$U2, y = .data$U1, yend = .data$U1), lyrs, colour = "black")
    } else {
      stop("Unrecognised setting for layers.2d argument. Valid choices are 'convex' (default), 'union' and 'detailed'.")
    }
  }

  jd <- c(1, which(sapply(1:nrow(path$W_t), function(i) { any(path$W_t[i,dim]!=path$W_tm[i,dim]) })), length(path$t))
  for(j in 2:length(jd)) {
    if(jd[j]-jd[j-1] > 0)
      p <- p + geom_path(aes(x = .data$x2, y = .data$x1), tibble(x1 = c(path$W_t[jd[j-1]:(jd[j]-1),dim[1]], path$W_tm[jd[j],dim[1]]), x2 = c(path$W_t[jd[j-1]:(jd[j]-1),dim[2]], path$W_tm[jd[j],dim[2]])), colour = "grey")
  }
  for(t.segend in path$labels$seg.end) {
    p <- p + geom_path(aes(x = .data$x2, y = .data$x1), tibble(x1 = c(path$W_tm[path$t %in% t.segend,dim[1]], path$W_t[path$t %in% t.segend,dim[1]]), x2 = c(path$W_tm[path$t %in% t.segend,dim[2]], path$W_t[path$t %in% t.segend,dim[2]])), colour = "grey", linetype = "dotted")
  }

  p <- p + geom_point(aes(x = .data$x2, y = .data$x1), tibble(x1 = path$W_t[,dim[1]], x2 = path$W_t[,dim[2]]), colour = "black", size = 0.1)
  jd <- jd[c(-1, -length(jd))]
  if(length(jd) > 0) {
    p <- p + geom_point(aes(x = .data$x2, y = .data$x1), tibble(x1 = path$W_tm[jd,dim[1]], x2 = path$W_tm[jd,dim[2]]), colour = "black", size = 1, shape = 1)
  }

  p + xlab("") + ylab("")
}

#' @rdname plot.BrownianMotion
#' @export
plot.BrownianMotionNd <- function(x, y, ...) {
  opts <- list(...)
  bm <- x

  # What dimensions to plot?
  if(!is.null(opts[["dims"]])) {
    dims <- opts[["dims"]]
  } else {
    dims <- 1:bm$dim
  }

  # Figure out plotting limits and extract visible path
  if(!is.null(opts[["t.lim"]])) {
    t.lim <- opts[["t.lim"]]
  } else {
    t.lim <- c(bm$Z.bm[[1]]$labels$start, bm$Z.bm[[1]]$labels$end)
  }
  t <- bm[,"t"]
  l <- max(min(t), t[t < t.lim[1]])
  r <- min(max(t), t[t > t.lim[2]])
  path.1d <- bm[t=l:r,]
  #path.2d <- bm[t=max(t[t <= t.lim[1]]):min(t[t > t.lim[2]]),]
  path.2d <- bm[t=t.lim[1]:t.lim[2],]
  # FIXME
  for(iii in seq_along(path.2d$layers$inner.cube)) {
    path.2d$layers$inner.cube[[iii]] <- unname(path.2d$layers$inner.cube[[iii]])
  }
  not.na <- !is.na(path.2d$W_t[,1])
  path.2d$t <- path.2d$t[not.na]
  path.2d$W_t <- path.2d$W_t[not.na,,drop = FALSE]
  path.2d$W_tm <- path.2d$W_tm[not.na,,drop = FALSE]
  # Correct inclusion of layers at the right end of interval
  # path.2d$layers <- path.2d$layers[path.2d$layers$t.l < t.lim[2],]

  # Layer display
  if(!is.null(opts[["layers.2d"]])) {
    layers.2d <- opts[["layers.2d"]]
  } else {
    layers.2d <- "convex"
  }

  # Layers visible?
  if(nrow(path.1d$layers) > 0) {
    if(t.lim[2] <= max(path.1d$layers$t.u) &&
       t.lim[1] >= min(path.1d$layers$t.l) &&
       all(head(sort(path.1d$layers$t.u), -1) == tail(sort(path.1d$layers$t.l), -1))) {
      show.layers <- TRUE
    } else {
      warning("an interval within the plot time limits requested has no layer information, so all layers suppressed on 2D plots.")
      show.layers <- FALSE
    }
  } else {
    show.layers <- FALSE
  }

  # Pairs or not
  if(!is.null(opts[["pairs"]])) {
    pairs <- opts[["pairs"]]
  } else {
    pairs <- TRUE
  }

  if(length(dims) == 1) {
    return(plot_diag_(path.1d, dims[1], t.lim))
  } else if(pairs) {
    p <- plot_diag_(path.1d, dims[1], t.lim)
    for(i in 1:length(dims)) {
      for(j in 1:length(dims)) {
        if(i == 1 && j == 1) {
          next
        }
        if(i != j) {
          p <- p + plot_offdiag_(path.2d, dims[c(i,j)], t.lim, layers.2d, show.layers)
        } else {
          p <- p + plot_diag_(path.1d, dims[i], t.lim)
        }
      }
    }
  } else {
    combs <- combn(dims, 2)
    p <- plot_offdiag_(path.2d, combs[,1], t.lim, layers.2d, show.layers) +
      xlab(paste0("W(",combs[2,1],")")) +
      ylab(paste0("W(",combs[1,1],")"))
    if(length(dims) > 2) {
      for(i in 2:ncol(combs)) {
        p <- p + plot_offdiag_(path.2d, combs[,i], t.lim, layers.2d, show.layers) +
          xlab(paste0("W(",combs[2,i],")")) +
          ylab(paste0("W(",combs[1,i],")"))
      }
    }
  }

  p
}



#' @rdname plot.BrownianMotion
#' @export
plot.BrownianMotionNdZ <- function(x, y, ...) {
  opts <- list(...)
  bm <- x[[1]]

  # What dimensions to plot?
  if(!is.null(opts[["dims"]])) {
    dims <- opts[["dims"]]
  } else {
    dims <- 1:bm$dim
  }

  # Figure out plotting limits and extract visible path
  if(!is.null(opts[["t.lim"]])) {
    t.lim <- opts[["t.lim"]]
  } else {
    t.lim <- c(bm$Z.bm[[1]]$labels$start, bm$Z.bm[[1]]$labels$end)
  }
  t <- bm[,"t"]
  l <- max(min(t), t[t < t.lim[1]])
  r <- min(max(t), t[t > t.lim[2]])
  path.1d <- bm$Z[t=l:r,]
  #path.2d <- bm[t=max(t[t <= t.lim[1]]):min(t[t > t.lim[2]]),]
  path.2d <- bm$Z[t=t.lim[1]:t.lim[2],]
  not.na <- !is.na(path.2d$Z_t[,1])
  path.2d$t <- path.2d$t[not.na]
  path.2d$W_t <- path.2d$Z_t[not.na,,drop = FALSE]
  path.2d$W_tm <- path.2d$Z_tm[not.na,,drop = FALSE]
  # Correct inclusion of layers at the right end of interval
  # path.2d$layers <- path.2d$layers[path.2d$layers$t.l < t.lim[2],]

  # Layer display
  if(!is.null(opts[["layers.2d"]])) {
    layers.2d <- opts[["layers.2d"]]
  } else {
    layers.2d <- "convex"
  }

  # Layers visible?
  if(nrow(path.1d$layers) > 0) {
    if(t.lim[2] <= max(path.1d$layers$t.u) &&
       t.lim[1] >= min(path.1d$layers$t.l) &&
       all(head(sort(path.1d$layers$t.u), -1) == tail(sort(path.1d$layers$t.l), -1))) {
      show.layers <- TRUE
    } else {
      warning("an interval within the plot time limits requested has no layer information, so all layers suppressed on 2D plots.")
      show.layers <- FALSE
    }
  } else {
    show.layers <- FALSE
  }

  # Pairs or not
  if(!is.null(opts[["pairs"]])) {
    pairs <- opts[["pairs"]]
  } else {
    pairs <- TRUE
  }

  if(length(dims) == 1) {
    return(plot(bm$Z.bm[[dims]], t.lim = t.lim) + ylab(paste0("Z(",dims[1],")")))
  } else if(pairs) {
    p <- plot(bm$Z.bm[[dims[1]]], t.lim = t.lim) + ylab(paste0("Z(",dims[1],")"))
    for(i in 1:length(dims)) {
      for(j in 1:length(dims)) {
        if(i == 1 && j == 1) {
          next
        }
        if(i != j) {
          p <- p + plot_offdiag_(path.2d, dims[c(i,j)], t.lim, layers.2d, show.layers)
        } else {
          p <- p + (plot(bm$Z.bm[[dims[i]]], t.lim = t.lim) + ylab(paste0("Z(",dims[i],")")))
        }
      }
    }
  } else {
    combs <- combn(dims, 2)
    p <- plot_offdiag_(path.2d, combs[,1], t.lim, layers.2d, show.layers) +
      xlab(paste0("Z(",combs[2,1],")")) +
      ylab(paste0("Z(",combs[1,1],")"))
    if(length(dims) > 2) {
      for(i in 2:ncol(combs)) {
        p <- p + plot_offdiag_(path.2d, combs[,i], t.lim, layers.2d, show.layers) +
          xlab(paste0("Z(",combs[2,i],")")) +
          ylab(paste0("Z(",combs[1,i],")"))
      }
    }
  }

  p
}



#' Plot a Brownian motion path and its layers
#'
#' Plots the skeleton of a Brownian motion path with \pkg{ggplot2}. Black dots
#' are the simulated points, joined by grey lines; these lines are only a guide
#' and are not the path itself. Layers are drawn as horizontal bounds coloured
#' by type (purple for Bessel, blue for intersection, red for localised), with
#' solid lines for bounds the path stays within and dashed or dotted lines for
#' inner bounds. Bounds given by the user to \code{\link{first_passage}} are
#' drawn in green. Jumps (see \code{\link{add_segment}}) are shown with open
#' circles at the left limit and a dotted vertical line.
#'
#' A multi-dimensional path is drawn as a grid of plots (combined with
#' \pkg{patchwork}): each component against time on the diagonal, and each pair
#' of components against each other off the diagonal, where the red region
#' contains the path over the plotted time range and the black box bounds it.
#' \code{plot(bm$Z)} draws the transformed, independent components instead.
#'
#' @param x a Brownian motion object.
#' @param y unused.
#' @param ... plotting options:
#'   \describe{
#'     \item{\code{t.lim}}{a time range \code{c(from, to)} to plot.}
#'     \item{\code{hide.user}}{(one dimension) \code{TRUE} to hide the green
#'       user-specified bounds.}
#'     \item{\code{dims}}{(several dimensions) which components to plot.}
#'     \item{\code{pairs}}{(several dimensions) \code{FALSE} to draw only the
#'       pairwise plots, without the plots against time.}
#'     \item{\code{layers.2d}}{(several dimensions) how to draw layers in the
#'       pairwise plots: \code{"convex"} (default) for the convex hull over the
#'       time range, \code{"detailed"} for each layer separately, or
#'       \code{"union"} for their union (requires the \pkg{sf} package).}
#'   }
#'
#' @return A \pkg{ggplot2} plot (a \pkg{patchwork} of plots for several
#'   dimensions), which is drawn when printed.
#'
#' @examples
#' set.seed(1)
#' bm <- create_bm()
#' bm |> sim(1:4) |> layers(0, 2, type = "bessel") |> first_passage(l = -3, u = 3)
#' plot(bm)
#' plot(bm, t.lim = c(0, 2))
#' plot(bm, hide.user = TRUE)
#'
#' bm2 <- create_bm(dim = 2, cov = matrix(c(1, 0.5, 0.5, 1), 2))
#' bm2 |> sim(1:2) |> layers(0, 2)
#' plot(bm2)
#' plot(bm2, layers.2d = "detailed", pairs = FALSE)
#' plot(bm2$Z)
#' @export
plot.BrownianMotion <- function(x, y, ...) {
  opts <- list(...)

  bm <- x

  p <- ggplot()

  jd <- c(1, which(bm$W_t != bm$W_tm), length(bm$t))
  for(j in 2:length(jd)) {
    if(j != length(jd))
      p <- p + geom_vline(xintercept = bm$t[jd[j]], linetype = "dotted")
    if(jd[j]-jd[j-1] > 0)
      p <- p + geom_line(aes(x = .data$t, y = .data$x), tibble(t = bm$t[jd[j-1]:jd[j]], x = c(bm$W_t[jd[j-1]:(jd[j]-1)], bm$W_tm[jd[j]])), colour = "grey")
  }

  p <- p + geom_point(aes(x = .data$t, y = .data$x), tibble(t = bm$t, x = bm$W_t), colour = "black", size = 0.1)
  jd <- jd[c(-1, -length(jd))]
  if(length(jd) > 0) {
    p <- p + geom_point(aes(x = .data$t, y = .data$x), tibble(t = bm$t[jd], x = bm$W_tm[jd]), colour = "black", size = 1, shape = 1)
  }

  localised <- bm$layers[bm$layers$type == "localised",]
  intersection <- bm$layers[bm$layers$type == "intersection",]
  bessel <- bm$layers[bm$layers$type == "bessel",]
  localised.bb <- bm$layers[bm$layers$type == "localised-bb",]
  intersection.bb <- bm$layers[bm$layers$type == "intersection-bb",]
  bessel.bb <- bm$layers[bm$layers$type == "bessel-bb",]

  if(is.null(opts[["hide.user"]]) || !opts[["hide.user"]]) {
    if(nrow(bm$user.layers) > 0) {
      p <- p +
        geom_segment(aes(x = .data$t.l, xend = .data$t.u, y = .data$L, yend = .data$L), bm$user.layers, colour = "green", linewidth = 1.3) +
        geom_segment(aes(x = .data$t.l, xend = .data$t.u, y = .data$U, yend = .data$U), bm$user.layers, colour = "green", linewidth = 1.3)
    }
  }

  if(nrow(localised) > 0) {
    p <- p +
      geom_segment(aes(x = .data$t.l, xend = .data$t.u, y = .data$Ld, yend = .data$Ld), localised, colour = "red") + #, linetype = "dashed") +
      geom_segment(aes(x = .data$t.l, xend = .data$t.u, y = .data$Uu, yend = .data$Uu), localised, colour = "red") + #, linetype = "dashed") +
      geom_segment(aes(x = .data$t.l, xend = .data$t.u, y = .data$Ud, yend = .data$Ud), localised, colour = "red", linetype = ifelse(localised$Ud.hard, "longdash", "dotted")) +
      geom_segment(aes(x = .data$t.l, xend = .data$t.u, y = .data$Lu, yend = .data$Lu), localised, colour = "red", linetype = ifelse(localised$Lu.hard, "longdash", "dotted"))
  }

  if(nrow(intersection) > 0) {
    p <- p +
      geom_segment(aes(x = .data$t.l, xend = .data$t.u, y = .data$Ld, yend = .data$Ld), intersection, colour = "blue") +
      geom_segment(aes(x = .data$t.l, xend = .data$t.u, y = .data$Uu, yend = .data$Uu), intersection, colour = "blue") +
      geom_segment(aes(x = .data$t.l, xend = .data$t.u, y = .data$Ud, yend = .data$Ud), intersection, colour = "blue", linetype = ifelse(intersection$Ud.hard, "longdash", "dotted")) +
      geom_segment(aes(x = .data$t.l, xend = .data$t.u, y = .data$Lu, yend = .data$Lu), intersection, colour = "blue", linetype = ifelse(intersection$Lu.hard, "longdash", "dotted"))
  }

  if(nrow(bessel) > 0) {
    p <- p +
      geom_segment(aes(x = .data$t.l, xend = .data$t.u, y = .data$Ld, yend = .data$Ld), bessel, colour = "purple") +
      geom_segment(aes(x = .data$t.l, xend = .data$t.u, y = .data$Uu, yend = .data$Uu), bessel, colour = "purple") +
      geom_segment(aes(x = .data$t.l, xend = .data$t.u, y = .data$Ud, yend = .data$Ud), bessel, colour = "purple", linetype = ifelse(bessel$Ud.hard, "longdash", "dotted")) +
      geom_segment(aes(x = .data$t.l, xend = .data$t.u, y = .data$Lu, yend = .data$Lu), bessel, colour = "purple", linetype = ifelse(bessel$Lu.hard, "longdash", "dotted"))
  }

  if(nrow(localised.bb) > 0) {
    localised.bb <- dplyr::left_join(localised.bb,
                                     bm$bb.local$layers,
                                     by = c("t.l", "t.u"),
                                     suffix = c("",".bb"))

    p <- p +
      geom_segment(aes(x = .data$t.l, xend = .data$t.u, y = .data$Ld.bb+.data$Ld, yend = .data$Ld.bb+.data$Uu), localised.bb, colour = "red") +
      geom_segment(aes(x = .data$t.l, xend = .data$t.u, y = .data$Uu.bb+.data$Ld, yend = .data$Uu.bb+.data$Uu), localised.bb, colour = "red") +
      geom_segment(aes(x = .data$t.l, xend = .data$t.u, y = .data$Ud.bb+.data$Ld, yend = .data$Ud.bb+.data$Uu), localised.bb, colour = "red", linetype = ifelse(localised.bb$Ud.hard, "longdash", "dotted")) +
      geom_segment(aes(x = .data$t.l, xend = .data$t.u, y = .data$Lu.bb+.data$Ld, yend = .data$Lu.bb+.data$Uu), localised.bb, colour = "red", linetype = ifelse(localised.bb$Lu.hard, "longdash", "dotted"))
    if(is.null(opts[["hide.user"]]) || !opts[["hide.user"]]) {
      p <- p +
        geom_segment(aes(x = .data$t.l, xend = .data$t.u, y = pmin(.data$Ld.bb+.data$Ld,.data$Ld.bb+.data$Uu), yend = pmin(.data$Ld.bb+.data$Ld,.data$Ld.bb+.data$Uu)), localised.bb, colour = "red", alpha = 0.75) +
        geom_segment(aes(x = .data$t.l, xend = .data$t.u, y = pmax(.data$Uu.bb+.data$Ld,.data$Uu.bb+.data$Uu), yend = pmax(.data$Uu.bb+.data$Ld,.data$Uu.bb+.data$Uu)), localised.bb, colour = "red", alpha = 0.75)
    }
  }

  if(nrow(intersection.bb) > 0) {
    intersection.bb <- dplyr::left_join(intersection.bb,
                                        bm$bb.local$layers,
                                        by = c("t.l", "t.u"),
                                        suffix = c("",".bb"))

    p <- p +
      geom_segment(aes(x = .data$t.l, xend = .data$t.u, y = .data$Ld.bb+.data$Ld, yend = .data$Ld.bb+.data$Uu), intersection.bb, colour = "blue") +
      geom_segment(aes(x = .data$t.l, xend = .data$t.u, y = .data$Uu.bb+.data$Ld, yend = .data$Uu.bb+.data$Uu), intersection.bb, colour = "blue") +
      geom_segment(aes(x = .data$t.l, xend = .data$t.u, y = .data$Ud.bb+.data$Ld, yend = .data$Ud.bb+.data$Uu), intersection.bb, colour = "blue", linetype = ifelse(intersection.bb$Ud.hard, "longdash", "dotted")) +
      geom_segment(aes(x = .data$t.l, xend = .data$t.u, y = .data$Lu.bb+.data$Ld, yend = .data$Lu.bb+.data$Uu), intersection.bb, colour = "blue", linetype = ifelse(intersection.bb$Lu.hard, "longdash", "dotted"))
    if(is.null(opts[["hide.user"]]) || !opts[["hide.user"]]) {
      p <- p +
        geom_segment(aes(x = .data$t.l, xend = .data$t.u, y = pmin(.data$Ld.bb+.data$Ld,.data$Ld.bb+.data$Uu), yend = pmin(.data$Ld.bb+.data$Ld,.data$Ld.bb+.data$Uu)), intersection.bb, colour = "blue", alpha = 0.75) +
        geom_segment(aes(x = .data$t.l, xend = .data$t.u, y = pmax(.data$Uu.bb+.data$Ld,.data$Uu.bb+.data$Uu), yend = pmax(.data$Uu.bb+.data$Ld,.data$Uu.bb+.data$Uu)), intersection.bb, colour = "blue", alpha = 0.75)
    }
  }

  if(nrow(bessel.bb) > 0) {
    bessel.bb <- dplyr::left_join(bessel.bb,
                                  bm$bb.local$layers,
                                  by = c("t.l", "t.u"),
                                  suffix = c("",".bb"))

    p <- p +
      geom_segment(aes(x = .data$t.l, xend = .data$t.u, y = .data$Ld.bb+.data$Ld, yend = .data$Ld.bb+.data$Uu), bessel.bb, colour = "purple") +
      geom_segment(aes(x = .data$t.l, xend = .data$t.u, y = .data$Uu.bb+.data$Ld, yend = .data$Uu.bb+.data$Uu), bessel.bb, colour = "purple") +
      geom_segment(aes(x = .data$t.l, xend = .data$t.u, y = .data$Ud.bb+.data$Ld, yend = .data$Ud.bb+.data$Uu), bessel.bb, colour = "purple", linetype = ifelse(bessel.bb$Ud.hard, "longdash", "dotted")) +
      geom_segment(aes(x = .data$t.l, xend = .data$t.u, y = .data$Lu.bb+.data$Ld, yend = .data$Lu.bb+.data$Uu), bessel.bb, colour = "purple", linetype = ifelse(bessel.bb$Lu.hard, "longdash", "dotted"))
    if(is.null(opts[["hide.user"]]) || !opts[["hide.user"]]) {
      p <- p +
        geom_segment(aes(x = .data$t.l, xend = .data$t.u, y = pmin(.data$Ld.bb+.data$Ld,.data$Ld.bb+.data$Uu), yend = pmin(.data$Ld.bb+.data$Ld,.data$Ld.bb+.data$Uu)), bessel.bb, colour = "purple", alpha = 0.75) +
        geom_segment(aes(x = .data$t.l, xend = .data$t.u, y = pmax(.data$Uu.bb+.data$Ld,.data$Uu.bb+.data$Uu), yend = pmax(.data$Uu.bb+.data$Ld,.data$Uu.bb+.data$Uu)), bessel.bb, colour = "purple", alpha = 0.75)
    }
  }

  if(!is.null(opts[["t.lim"]])) {
    t.lim <- opts[["t.lim"]]

    lyrs <- bm$layers[!(bm$layers$type %in% c("localised-bb","intersection-bb","bessel-bb")),]
    ymin <- min(lyrs[lyrs$t.u>t.lim[1] & lyrs$t.l<t.lim[2],]$Ld, bm$W_t[bm$t>=t.lim[1] & bm$t<=t.lim[2]])
    ymax <- max(lyrs[lyrs$t.u>t.lim[1] & lyrs$t.l<t.lim[2],]$Uu, bm$W_t[bm$t>=t.lim[1] & bm$t<=t.lim[2]])

    bb.lyrs <- rbind(localised.bb, intersection.bb, bessel.bb)
    if(nrow(bb.lyrs) > 0) {
      bb.lyrs$min <- pmin(bb.lyrs$Ld.bb+bb.lyrs$Ld,bb.lyrs$Ld.bb+bb.lyrs$Uu)
      bb.lyrs$max <- pmax(bb.lyrs$Uu.bb+bb.lyrs$Ld,bb.lyrs$Uu.bb+bb.lyrs$Uu)
      ymin <- min(ymin,
                  bb.lyrs[bb.lyrs$t.u>t.lim[1] & bb.lyrs$t.l<t.lim[2],]$min)
      ymax <- max(ymax,
                  bb.lyrs[bb.lyrs$t.u>t.lim[1] & bb.lyrs$t.l<t.lim[2],]$max)
    }

    return(p + coord_cartesian(xlim = opts[["t.lim"]], ylim = c(ymin,ymax)))
  } else {
    return(p)
  }
}
