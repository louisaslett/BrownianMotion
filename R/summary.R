#' Summarise a Brownian motion object
#'
#' Summarises the skeleton of a Brownian motion path, its layers, jumps and
#' labels, and the default settings used for further simulation.
#'
#' The values of the path are summarised at the skeleton points only (including
#' the left limits at any jumps), so the minimum and maximum are those of the
#' simulated points, not of the whole path. Over the parts of the time range
#' covered by layers, however, the whole path is known to lie within the layer
#' bounds, and the summary reports the overall bounds this gives. For a
#' multi-dimensional path these are the bounds of the smallest box, with sides
#' parallel to the axes, which contains every layer.
#'
#' @param object a Brownian motion object, or the transformed view \code{bm$Z}
#'   of a multi-dimensional one.
#' @param x a summary, as returned by \code{summary()}.
#' @param digits the minimum number of significant digits to show, as in
#'   \code{\link[base]{format}}. Defaults to \code{max(3, getOption("digits") -
#'   3)}.
#' @param width the width of the output in characters. Defaults to the width of
#'   the console.
#' @param ... unused.
#'
#' @return \code{summary()} returns an object of class
#'   \code{"summary.BrownianMotion"} (one dimension) or
#'   \code{"summary.BrownianMotionNd"} (several dimensions, and \code{bm$Z}),
#'   which is a list with elements:
#'   \describe{
#'     \item{\code{dim}}{the dimension of the path.}
#'     \item{\code{view}}{(several dimensions) \code{"W"} for the path, or
#'       \code{"Z"} for its transformed view \code{bm$Z}.}
#'     \item{\code{cov}}{(several dimensions) the covariance matrix of the
#'       components summarised.}
#'     \item{\code{n.points}}{the number of skeleton points.}
#'     \item{\code{t.range}}{the first and last times of the skeleton.}
#'     \item{\code{origin}}{the number of skeleton points of each origin:
#'       \code{"fixed"}, \code{"sim"}, \code{"passage"} and \code{"other"} (see
#'       \code{\link[=print.BrownianMotion]{print}}).}
#'     \item{\code{spacing}}{the minimum, median and maximum time between
#'       consecutive skeleton points, or \code{NULL} for a single point.}
#'     \item{\code{path}}{(one dimension) a data frame giving the time and value
#'       of the path at its start, end, minimum and maximum.}
#'     \item{\code{components}}{(several dimensions) a data frame with a row per
#'       component giving its value at the start and end, its minimum and
#'       maximum, and the lower (\code{L}) and upper (\code{U}) layer bounds.}
#'     \item{\code{jumps}}{a data frame with a row per jump, giving its time
#'       \code{t} and \code{size} (in one dimension also the left limit
#'       \code{W_tm} and value \code{W_t}; in several dimensions the size is
#'       the Euclidean length of the jump).}
#'     \item{\code{layers}}{(one dimension) a data frame with a row per layer
#'       type, giving the number of intervals, the total time they cover, and
#'       the mean distance between their bounds. In several dimensions, a matrix
#'       of the number of intervals of each type (rows) in each component of
#'       \code{bm$Z} (columns).}
#'     \item{\code{n.layers}}{the number of layer intervals.}
#'     \item{\code{coverage}}{the proportion of the time range covered by
#'       layers (\code{NA} for a single point).}
#'     \item{\code{gaps}}{a data frame of the intervals, \code{from} and
#'       \code{to}, not covered by any layer.}
#'     \item{\code{bounds}}{the bounds which contain the whole path over the
#'       layered time: a vector \code{c(L = , U = )} in one dimension, and part
#'       of \code{components} in several. \code{NULL} if there are no layers.}
#'     \item{\code{labels}}{the number of points with each user label.}
#'     \item{\code{passages}}{(one dimension) the number of first passages
#'       requested with \code{\link{first_passage}}.}
#'     \item{\code{settings}}{the defaults \code{refine}, \code{mult} and
#'       \code{prefer} used by \code{\link{sim}} and \code{\link{layers}} (one
#'       value per component in several dimensions).}
#'   }
#'   \code{print()} returns \code{x} invisibly.
#'
#' @seealso \code{\link[=print.BrownianMotion]{print}} for a table of the
#'   skeleton, and \code{\link[=BrownianMotion-indexing]{indexing}} and
#'   \code{\link{layers}} for the details of each point and layer.
#'
#' @examples
#' set.seed(1)
#' bm <- create_bm()
#' bm |> sim(c(1, 2, peak = 3)) |> layers(0, 2, type = "bessel") |>
#'   first_passage(delta = 1)
#' summary(bm)
#'
#' # The summary can also be used programmatically
#' s <- summary(bm)
#' s$coverage
#' s$bounds
#'
#' bm2 <- create_bm(dim = 2, cov = matrix(c(1, 0.5, 0.5, 1), 2))
#' bm2 |> sim(1:2) |> layers(0, 2)
#' summary(bm2)
#' summary(bm2$Z)
#' @export
summary.BrownianMotion <- function(object, ...) {
  bm <- object
  t <- bm$t
  lyrs <- layer_bounds_(bm)
  gaps <- layer_gaps_(lyrs$t.l, lyrs$t.u, min(t), max(t))

  # Extremes of the skeleton, including left limits at jumps
  values <- c(bm$W_t, bm$W_tm)
  times <- c(t, t)
  path <- data.frame(t = c(t[1], t[length(t)], times[which.min(values)], times[which.max(values)]),
                     W_t = c(bm$W_t[1], bm$W_t[length(t)], min(values), max(values)),
                     row.names = c("start", "end", "min", "max"))

  jumps <- which(bm$W_t != bm$W_tm)
  jumps <- data.frame(t = t[jumps],
                      W_tm = unname(bm$W_tm[jumps]),
                      W_t = unname(bm$W_t[jumps]),
                      size = unname(bm$W_t[jumps] - bm$W_tm[jumps]))

  types <- layer_type_order_(lyrs$type)
  layer.types <- data.frame(type = types,
                            intervals = vapply(types, function(ty) { sum(lyrs$type == ty) }, integer(1)),
                            time = vapply(types, function(ty) { with(lyrs[lyrs$type == ty, ], sum(t.u - t.l)) }, numeric(1)),
                            width = vapply(types, function(ty) { with(lyrs[lyrs$type == ty, ], mean(U - L)) }, numeric(1)),
                            row.names = NULL, stringsAsFactors = FALSE)

  res <- list(dim = 1L,
              n.points = length(t),
              t.range = range(t),
              origin = origin_counts_(bm$labels, t),
              spacing = spacing_(t),
              path = path,
              jumps = jumps,
              layers = layer.types,
              n.layers = nrow(lyrs),
              coverage = coverage_(gaps, t),
              gaps = gaps,
              bounds = if(nrow(lyrs) > 0) c(L = min(lyrs$L), U = max(lyrs$U)) else NULL,
              labels = label_counts_(bm$labels),
              passages = nrow(bm$user.layers),
              settings = list(refine = bm$refine, mult = bm$mult, prefer = bm$prefer))
  class(res) <- "summary.BrownianMotion"
  res
}

#' @rdname summary.BrownianMotion
#' @export
summary.BrownianMotionNd <- function(object, ...) {
  summary_nd_(object, view = "W")
}

#' @rdname summary.BrownianMotion
#' @export
summary.BrownianMotionNdZ <- function(object, ...) {
  summary_nd_(object[[1]], view = "Z")
}

summary_nd_ <- function(bm, view) {
  t <- bm$Z.bm[[1]]$t
  i <- seq_along(t)
  if(view == "W") {
    X <- getWtNd_(bm, i)
    Xm <- getWtmNd_(bm, i)
    cov <- bm$cov
  } else {
    X <- getZtNd_(bm, i)
    Xm <- getZtmNd_(bm, i)
    cov <- diag(bm$dim)
  }
  comps <- paste0(view, "[", seq_len(bm$dim), "]")

  # Outer bounding box over all layers, in the view summarised
  lyrs <- nd_layer_types_(bm)
  if(nrow(lyrs) > 0) {
    box <- nd_layer_box_(bm, lyrs$t.l, transform = (view == "W"))
    L <- apply(box$L, 2, min)
    U <- apply(box$U, 2, max)
  } else {
    L <- U <- rep(NA_real_, bm$dim)
  }
  components <- data.frame(start = X[1, ],
                           end = X[length(t), ],
                           min = apply(rbind(X, Xm), 2, min),
                           max = apply(rbind(X, Xm), 2, max),
                           L = L,
                           U = U,
                           row.names = comps)

  jumps <- which(rowSums(X != Xm) > 0)
  jumps <- data.frame(t = t[jumps],
                      size = sqrt(rowSums((X[jumps, , drop = FALSE] - Xm[jumps, , drop = FALSE])^2)))

  # Layer types are those of the independent components in Z
  types <- layer_type_order_(unlist(lapply(bm$Z.bm, function(Z.bm) { Z.bm$layers$type })))
  layer.types <- vapply(bm$Z.bm, function(Z.bm) {
    vapply(types, function(ty) { sum(Z.bm$layers$type == ty) }, integer(1))
  }, integer(length(types)))
  layer.types <- matrix(layer.types, nrow = length(types),
                        dimnames = list(types, paste0("Z[", seq_len(bm$dim), "]")))

  gaps <- layer_gaps_(lyrs$t.l, lyrs$t.u, min(t), max(t))

  res <- list(dim = bm$dim,
              view = view,
              cov = cov,
              n.points = length(t),
              t.range = range(t),
              origin = origin_counts_(bm$Z.bm[[1]]$labels, t),
              spacing = spacing_(t),
              components = components,
              jumps = jumps,
              layers = layer.types,
              n.layers = nrow(lyrs),
              coverage = coverage_(gaps, t),
              gaps = gaps,
              labels = label_counts_(bm$Z.bm[[1]]$labels),
              settings = data.frame(refine = vapply(bm$Z.bm, function(Z.bm) { Z.bm$refine }, logical(1)),
                                    mult = vapply(bm$Z.bm, function(Z.bm) { Z.bm$mult }, numeric(1)),
                                    prefer = vapply(bm$Z.bm, function(Z.bm) { Z.bm$prefer }, character(1)),
                                    row.names = paste0("Z[", seq_len(bm$dim), "]"),
                                    stringsAsFactors = FALSE))
  class(res) <- "summary.BrownianMotionNd"
  res
}

#' @rdname summary.BrownianMotion
#' @export
print.summary.BrownianMotion <- function(x, digits = NULL, ...) {
  cat(format(x, digits = digits, ...), sep = "\n")
  invisible(x)
}

#' @rdname summary.BrownianMotion
#' @export
print.summary.BrownianMotionNd <- function(x, digits = NULL, ...) {
  cat(format(x, digits = digits, ...), sep = "\n")
  invisible(x)
}

#' @rdname summary.BrownianMotion
#' @export
format.summary.BrownianMotion <- function(x, digits = NULL, width = cli::console_width(), ...) {
  digits <- default_digits_(digits)

  overview <- overview_lines_(x, digits, width)

  path <- x$path
  path.lines <- align_cols_(list(" " = rownames(path),
                                 t = format_num_(path$t, digits),
                                 W_t = format_num_(path$W_t, digits)),
                            left = " ", indent = "  ")

  c(header_rule_("Brownian motion summary", dim_text_(1), 72, width),
    overview,
    "",
    section_("Path at the skeleton points"),
    path.lines,
    "",
    layer_lines_(x, digits, width))
}

#' @rdname summary.BrownianMotion
#' @export
format.summary.BrownianMotionNd <- function(x, digits = NULL, width = cli::console_width(), ...) {
  digits <- default_digits_(digits)
  view <- x$view

  title <- if(view == "W") "Brownian motion summary" else "Brownian motion summary, transformed view Z"
  if(view == "W") {
    cov.lines <- c(section_("Covariance"),
                   paste0("  ", cov_text_(x$cov, digits)))
    if(x$dim <= 6) {
      cov.lines <- c(cov.lines, paste0("  ", utils::capture.output(print(x$cov, digits = digits))))
    }
  } else {
    cov.lines <- c(section_("Covariance"),
                   "  independent standard components, with W_t = Z_t %*% chol(cov)")
  }

  comp <- x$components
  cols <- list(" " = rownames(comp),
               start = format_num_(comp$start, digits),
               end = format_num_(comp$end, digits),
               min = format_num_(comp$min, digits),
               max = format_num_(comp$max, digits))
  if(x$n.layers > 0) {
    cols$`layer bounds` <- paste0("[", format_num_(comp$L, digits), ", ", format_num_(comp$U, digits), "]")
  }
  comp.lines <- align_cols_(cols, left = " ", indent = "  ")

  c(header_rule_(title, dim_text_(x$dim), 72, width),
    overview_lines_(x, digits, width),
    "",
    cov.lines,
    "",
    section_("Components at the skeleton points"),
    comp.lines,
    "",
    layer_lines_(x, digits, width))
}



# Helpers for summary() -------------------------------------------------------

# Lower (L) and upper (U) bounds of the axis-parallel box containing each of the
# layers starting at times t.l, as matrices with a row per layer. Each layer is
# a box in Z, whose image under W = Z %*% chol is bounded in each component by
# taking the smaller (or larger) end of every side.
nd_layer_box_ <- function(bm, t.l, transform = TRUE) {
  bounds <- lapply(bm$Z.bm, function(Z.bm) {
    lyrs <- layer_bounds_(Z.bm)
    lyrs[match(t.l, lyrs$t.l), c("L", "U")]
  })
  Lz <- sapply(bounds, function(b) { b$L })
  Uz <- sapply(bounds, function(b) { b$U })
  Lz <- matrix(Lz, nrow = length(t.l))
  Uz <- matrix(Uz, nrow = length(t.l))
  if(!transform) {
    return(list(L = Lz, U = Uz))
  }
  L <- U <- matrix(0, nrow = length(t.l), ncol = bm$dim)
  for(d in seq_len(bm$dim)) {
    a <- sweep(Lz, 2, bm$chol[, d], "*")
    b <- sweep(Uz, 2, bm$chol[, d], "*")
    L[, d] <- rowSums(pmin(a, b))
    U[, d] <- rowSums(pmax(a, b))
  }
  list(L = L, U = U)
}

layer_type_order_ <- function(types) {
  all.types <- c("bessel", "intersection", "localised", "bessel-bb", "intersection-bb", "localised-bb")
  c(all.types[all.types %in% types], sort(setdiff(unique(types), all.types)))
}

origin_counts_ <- function(labels, t) {
  table(factor(point_origin_(labels, t), levels = c("fixed", "sim", "passage", "other")),
        dnn = NULL)
}

spacing_ <- function(t) {
  if(length(t) < 2) {
    return(NULL)
  }
  dt <- diff(t)
  c(min = min(dt), median = stats::median(dt), max = max(dt))
}

coverage_ <- function(gaps, t) {
  span <- max(t) - min(t)
  if(span == 0) {
    return(NA_real_)
  }
  (span - sum(gaps$to - gaps$from)) / span
}

label_counts_ <- function(labels) {
  nms <- user_labels_(labels)
  vapply(labels[nms], function(x) { length(unique(x)) }, integer(1))
}

section_ <- function(title) {
  as.character(cli::style_bold(title))
}

# Lines of "key  value" pairs, with the keys aligned
kv_lines_ <- function(keys, values, width) {
  keys <- format(keys)
  lines <- paste0(keys, "  ", values)
  as.character(cli::ansi_strtrim(lines, width))
}

intervals_text_ <- function(from, to, digits) {
  paste0("[", fmt_(from, digits), ", ", fmt_(to, digits), "]")
}

# The overview lines shared by the one and multi-dimensional summaries
overview_lines_ <- function(x, digits, width) {
  origin.names <- c(fixed = "fixed", sim = "simulated", passage = "first passage", other = "other")
  origin <- x$origin[x$origin > 0]
  points <- paste0(x$n.points, if(length(origin) > 0) {
    paste0(" (", paste(origin, origin.names[names(origin)], collapse = ", "), ")")
  })

  if(x$n.points == 1) {
    time <- paste0("t = ", fmt_(x$t.range[1], digits))
  } else {
    time <- paste0(intervals_text_(x$t.range[1], x$t.range[2], digits),
                   ", length ", fmt_(diff(x$t.range), digits))
  }

  spacing <- if(is.null(x$spacing)) {
    "-"
  } else {
    paste0("min ", fmt_(x$spacing[["min"]], digits),
           ", median ", fmt_(x$spacing[["median"]], digits),
           ", max ", fmt_(x$spacing[["max"]], digits))
  }

  n.jumps <- nrow(x$jumps)
  jumps <- if(n.jumps == 0) {
    "none"
  } else if(n.jumps <= 3) {
    paste0(n.jumps, ", at t = ", paste(fmt_(x$jumps$t, digits), collapse = ", "),
           if(n.jumps == 1) " (size " else " (sizes ",
           paste(fmt_(x$jumps$size, digits), collapse = ", "), ")")
  } else {
    paste0(n.jumps, ", with sizes from ", fmt_(min(x$jumps$size), digits),
           " to ", fmt_(max(x$jumps$size), digits))
  }

  labels <- if(length(x$labels) == 0) {
    "none"
  } else {
    paste0(names(x$labels), " (", vapply(x$labels, count_text_, character(1), what = "point"), ")",
           collapse = ", ")
  }

  settings <- x$settings
  defaults <- paste0("refine = ", collapse_values_(settings$refine),
                     ", mult = ", collapse_values_(settings$mult),
                     ", prefer = ", collapse_values_(paste0('"', settings$prefer, '"')))

  kv_lines_(c("Time", "Points", "Spacing", "Jumps", "Labels", "Defaults"),
            c(time, points, spacing, jumps, labels, defaults),
            width)
}

# A single value if all are equal, otherwise all of them, as in c(1, 2)
collapse_values_ <- function(x) {
  if(length(unique(x)) == 1) {
    as.character(x[1])
  } else {
    paste0("c(", paste(x, collapse = ", "), ")")
  }
}

# The layers section shared by the one and multi-dimensional summaries
layer_lines_ <- function(x, digits, width) {
  if(x$n.layers == 0) {
    return(c(section_("Layers"), "  none"))
  }
  span.text <- if(is.na(x$coverage)) "" else paste0(", covering ", format_pct_(x$coverage), " of the time")
  lines <- section_(paste0("Layers: ", count_text_(x$n.layers, "interval"), span.text))

  if(x$dim == 1) {
    lt <- x$layers
    lines <- c(lines,
               align_cols_(list(" " = rep("", nrow(lt)),
                                type = style_layer_type_(lt$type),
                                intervals = as.character(lt$intervals),
                                time = format_num_(lt$time, digits),
                                `mean width` = format_num_(lt$width, digits)),
                           left = "type", indent = ""))
  } else {
    lt <- x$layers
    cols <- c(list(" " = rep("", nrow(lt)), type = style_layer_type_(rownames(lt))),
              lapply(seq_len(ncol(lt)), function(d) { as.character(lt[, d]) }))
    names(cols)[-(1:2)] <- colnames(lt)
    max.cols <- 2 + max(1, floor((width - 20) / 7))
    if(length(cols) > max.cols) {
      cols <- cols[seq_len(max.cols)]
    }
    lines <- c(lines,
               "  Intervals of each type, in each independent component of Z:",
               align_cols_(cols, left = "type", indent = ""))
  }

  gaps <- x$gaps
  if(nrow(gaps) > 0) {
    shown <- utils::head(seq_len(nrow(gaps)), 3)
    gap.text <- paste(intervals_text_(gaps$from[shown], gaps$to[shown], digits), collapse = ", ")
    if(nrow(gaps) > 3) {
      gap.text <- paste0(gap.text, " and ", nrow(gaps) - 3, " more")
    }
    lines <- c(lines, kv_lines_("  Not covered", gap.text, width))
  }

  if(x$dim == 1) {
    lines <- c(lines,
               wrap_(paste0("Over the layered time the whole path lies in ",
                            intervals_text_(x$bounds[["L"]], x$bounds[["U"]], digits), "."),
                     width, indent = 2, prefix = "  "))
    if(x$passages > 0) {
      lines <- c(lines, paste0("  ", count_text_(x$passages, "first passage"), " requested with first_passage()."))
    }
  } else {
    lines <- c(lines, wrap_("Over the layered time the whole path lies within the layer bounds above.",
                            width, indent = 2, prefix = "  "))
  }
  lines
}
