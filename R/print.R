#' Print a Brownian motion object
#'
#' Prints an overview of a Brownian motion path followed by its skeleton as a
#' table, with one row per simulated point. \code{format()} returns the same
#' output as a character vector, one element per line.
#'
#' The first lines give the number of skeleton points and the time range they
#' span, the number of layers and the proportion of that time range they cover,
#' and the number of jumps and user labels. For a multi-dimensional path the
#' covariance is also described.
#'
#' The table then lists, for each skeleton point:
#' \describe{
#'   \item{\code{t}, \code{W_t}}{the time and the value of the path (one column
#'     per component in several dimensions, or \code{Z_t} for the transformed
#'     view \code{bm$Z}).}
#'   \item{\code{W_tm}}{(one dimension, paths with jumps only) the left limit of
#'     the path at each jump. In several dimensions a \code{jump} column gives
#'     the (Euclidean) size of the jump instead.}
#'   \item{\code{origin}}{how the point arose: \code{"fixed"} for points given to
#'     \code{\link{create_bm}} or \code{\link{add_points}}, \code{"sim"} for
#'     points simulated with \code{\link{sim}} (or as the ends of a layer
#'     interval), \code{"passage"} for first passage times (see
#'     \code{\link{first_passage}}, which are also used to build localised
#'     layers) and \code{"other"} for anything else.}
#'   \item{\code{layer}, \code{bounds}}{the type of layer, if any, over the
#'     interval from this point to the next, and (in one dimension) the bounds
#'     \code{[L, U]} which contain the whole path over that interval. For a
#'     multi-dimensional path the layer types are those of the independent
#'     components of \code{bm$Z}; see \code{bm[, "layers"]} for the bounds.}
#'   \item{\code{labels}}{any user labels on the point (see
#'     \code{\link{add_labels}}).}
#' }
#' Long paths show only the first and last few points; use \code{n} to see
#' more. Layer types are coloured as in \code{\link[=plot.BrownianMotion]{plot}}
#' when the console supports colour. Columns which do not fit in the console
#' width are dropped or truncated.
#'
#' @param x a Brownian motion object, or the transformed view \code{bm$Z} of a
#'   multi-dimensional one.
#' @param n the number of skeleton points to show. By default every point is
#'   shown when there are at most 20, and 10 otherwise. Use \code{Inf} to show
#'   all of them.
#' @param digits the minimum number of significant digits to show, as in
#'   \code{\link[base]{format}}. Defaults to \code{max(3, getOption("digits") -
#'   3)}.
#' @param width the width of the output in characters. Defaults to the width of
#'   the console.
#' @param ... further arguments passed on from \code{print()} to
#'   \code{format()}, and otherwise ignored.
#'
#' @return \code{print()} returns \code{x} invisibly. \code{format()} returns a
#'   character vector with one element per line of output.
#'
#' @seealso \code{\link[=summary.BrownianMotion]{summary}} for a fuller
#'   summary, and \code{\link[=BrownianMotion-indexing]{indexing}} to extract
#'   the values shown.
#'
#' @examples
#' set.seed(1)
#' bm <- create_bm()
#' bm |> sim(c(1, 2, peak = 3)) |> layers(0, 2, type = "bessel")
#' bm
#'
#' # Long paths show only the first and last few points
#' bm |> first_passage(delta = 0.5) |> sim(seq(3.1, 6, by = 0.1))
#' bm
#' print(bm, n = 4, digits = 2)
#'
#' # Multi-dimensional paths, and their transformed view
#' bm2 <- create_bm(dim = 2, cov = matrix(c(1, 0.5, 0.5, 1), 2))
#' bm2 |> sim(1:2) |> layers(0, 2)
#' bm2
#' bm2$Z
#' @export
print.BrownianMotion <- function(x, n = NULL, digits = NULL, ...) {
  cat(format(x, n = n, digits = digits, ...), sep = "\n")
  invisible(x)
}

#' @rdname print.BrownianMotion
#' @export
print.BrownianMotionNd <- function(x, n = NULL, digits = NULL, ...) {
  cat(format(x, n = n, digits = digits, ...), sep = "\n")
  invisible(x)
}

#' @rdname print.BrownianMotion
#' @export
print.BrownianMotionNdZ <- function(x, n = NULL, digits = NULL, ...) {
  cat(format(x, n = n, digits = digits, ...), sep = "\n")
  invisible(x)
}

#' @rdname print.BrownianMotion
#' @export
format.BrownianMotion <- function(x, n = NULL, digits = NULL, width = cli::console_width(), ...) {
  bm <- x
  digits <- default_digits_(digits)
  n.pts <- length(bm$t)
  rows <- print_rows_(n.pts, n)

  lyrs <- layer_bounds_(bm)
  jumps <- which(bm$W_t != bm$W_tm)
  user.labels <- user_labels_(bm$labels)

  # Overview
  overview <- c(points_text_(bm$t, digits),
                layers_text_(lyrs, bm$t),
                if(length(jumps) > 0) count_text_(length(jumps), "jump"),
                if(length(user.labels) > 0) count_text_(length(user.labels), "user label"))

  # Skeleton table, for the rows shown
  t <- bm$t[rows]
  cols <- list(t = format_num_(t, digits),
               W_t = format_num_(bm$W_t[rows], digits))
  if(length(jumps) > 0) {
    W_tm <- bm$W_tm[rows]
    W_tm[!(rows %in% jumps)] <- NA
    cols$W_tm <- format_num_(W_tm, digits)
  }
  cols$origin <- point_origin_(bm$labels, t)
  if(nrow(lyrs) > 0) {
    j <- interval_layer_(lyrs, t)
    cols$layer <- style_layer_type_(lyrs$type[j])
    bounds <- paste0("[", format_num_(lyrs$L[j], digits), ", ", format_num_(lyrs$U[j], digits), "]")
    bounds[is.na(j)] <- ""
    cols$bounds <- bounds
  }
  if(length(user.labels) > 0) {
    cols$labels <- point_labels_(bm$labels, t)
  }

  tbl <- format_table_(cols, rows, n.pts, width,
                       left = c("origin", "layer", "labels"))

  notes <- c(if(nrow(lyrs) > 0) "Layers and bounds apply from each point to the next.",
             tbl$notes)

  c(header_rule_("Brownian motion", "1 dimension", max(tbl$width, 60), width),
    wrap_(paste(overview, collapse = sep_()), width),
    "",
    tbl$lines,
    footer_(notes, width))
}

#' @rdname print.BrownianMotion
#' @export
format.BrownianMotionNd <- function(x, n = NULL, digits = NULL, width = cli::console_width(), ...) {
  format_nd_(x, n = n, digits = digits, width = width, view = "W")
}

#' @rdname print.BrownianMotion
#' @export
format.BrownianMotionNdZ <- function(x, n = NULL, digits = NULL, width = cli::console_width(), ...) {
  format_nd_(x[[1]], n = n, digits = digits, width = width, view = "Z")
}

# Shared by the multi-dimensional path (view = "W") and its transformed view
# bm$Z (view = "Z"), which differ only in the values shown.
format_nd_ <- function(bm, n, digits, width, view) {
  digits <- default_digits_(digits)
  t.all <- bm$Z.bm[[1]]$t
  n.pts <- length(t.all)
  rows <- print_rows_(n.pts, n)
  labels <- bm$Z.bm[[1]]$labels
  user.labels <- user_labels_(labels)

  lyrs <- nd_layer_types_(bm)
  jumps <- which(rowSums(getZtNd_(bm, seq_len(n.pts)) != getZtmNd_(bm, seq_len(n.pts))) > 0)

  overview <- c(points_text_(t.all, digits),
                layers_text_(lyrs, t.all),
                if(length(jumps) > 0) count_text_(length(jumps), "jump"),
                if(length(user.labels) > 0) count_text_(length(user.labels), "user label"))
  if(view == "W") {
    cov.text <- paste0("Covariance: ", cov_text_(bm$cov, digits))
  } else {
    cov.text <- paste0("Independent standard components, with W_t = Z_t %*% chol(cov)")
  }

  t <- t.all[rows]
  if(view == "W") {
    X <- getWtNd_(bm, rows)
    Xm <- getWtmNd_(bm, rows)
  } else {
    X <- getZtNd_(bm, rows)
    Xm <- getZtmNd_(bm, rows)
  }
  cols <- list(t = format_num_(t, digits))
  comp.cols <- paste0(view, "_t[", seq_len(bm$dim), "]")
  for(d in seq_len(bm$dim)) {
    cols[[comp.cols[d]]] <- format_num_(X[, d], digits)
  }
  if(length(jumps) > 0) {
    size <- sqrt(rowSums((X - Xm)^2))
    size[!(rows %in% jumps)] <- NA
    cols$jump <- format_num_(size, digits)
  }
  cols$origin <- point_origin_(labels, t)
  if(nrow(lyrs) > 0) {
    cols$layer <- style_layer_type_(lyrs$type[interval_layer_(lyrs, t)])
  }
  if(length(user.labels) > 0) {
    cols$labels <- point_labels_(labels, t)
  }

  # Components which do not fit are dropped from the right, keeping the first
  tbl <- format_table_(cols, rows, n.pts, width,
                       left = c("origin", "layer", "labels"),
                       droppable = rev(comp.cols[-1]))

  notes <- c(if(nrow(lyrs) > 0) {
               if(view == "W") {
                 "Layers are those of the components of bm$Z, from each point to the next."
               } else {
                 "Layers apply from each point to the next."
               }
             },
             tbl$notes)

  title <- if(view == "W") "Brownian motion" else "Brownian motion, transformed view Z"
  c(header_rule_(title, dim_text_(bm$dim), max(tbl$width, 60), width),
    wrap_(paste(overview, collapse = sep_()), width),
    wrap_(cov.text, width),
    "",
    tbl$lines,
    footer_(notes, width))
}



# Helpers shared by format() and summary() ------------------------------------

# Labels used internally by the package, which are not shown as user labels
bm_internal_labels_ <- c("start", "end", "user", "seg.start", "seg.end", "fpt", "internal", "forced")

default_digits_ <- function(digits) {
  if(is.null(digits)) {
    digits <- max(3L, getOption("digits") - 3L)
  }
  if(!checkmate::test_int(digits, lower = 1, upper = 22)) {
    stop("digits must be a whole number between 1 and 22")
  }
  as.integer(digits)
}

# The names of user labels, ignoring internal and empty labels
user_labels_ <- function(labels) {
  nms <- names(labels)
  if(is.null(nms)) {
    return(character())
  }
  nms <- nms[!(nms %in% bm_internal_labels_) & nms != "" & !is.na(nms)]
  nms[vapply(labels[nms], length, integer(1)) > 0]
}

# The user labels at each of the times t, as comma separated strings
point_labels_ <- function(labels, t) {
  out <- character(length(t))
  for(l in user_labels_(labels)) {
    has <- t %in% labels[[l]]
    out[has] <- ifelse(out[has] == "", l, paste0(out[has], ", ", l))
  }
  out
}

# How each of the points at times t arose, from the internal labels
point_origin_ <- function(labels, t) {
  origin <- rep("other", length(t))
  origin[t %in% labels[["user"]]] <- "sim"
  origin[t %in% labels[["fpt"]]] <- "passage"
  origin[t %in% labels[["forced"]]] <- "fixed"
  origin
}

# Layers of a 1-dimensional path sorted by time, with the bounds L and U which
# contain the whole path over each interval.
layer_bounds_ <- function(bm) {
  lyrs <- bm$layers[order(bm$layers$t.l), ]
  L <- lyrs$Ld
  U <- lyrs$Uu
  bb <- grepl("-bb$", lyrs$type)
  if(any(bb)) {
    # Bridge layers bound the auxiliary path, and store the offset of the path
    # from it at either end (Ld on the left, Uu on the right), which is
    # interpolated linearly in between
    if(is.null(bm$bb.local)) {
      L[bb] <- U[bb] <- NA
    } else {
      aux <- bm$bb.local$layers
      idx <- match(lyrs$t.l[bb], aux$t.l)
      L[bb] <- aux$Ld[idx] + pmin(lyrs$Ld[bb], lyrs$Uu[bb])
      U[bb] <- aux$Uu[idx] + pmax(lyrs$Ld[bb], lyrs$Uu[bb])
    }
  }
  data.frame(type = lyrs$type, t.l = lyrs$t.l, t.u = lyrs$t.u, L = L, U = U,
             stringsAsFactors = FALSE)
}

# Layer intervals of a multi-dimensional path sorted by time, with the layer
# type of each (types which differ between components are joined with "/")
nd_layer_types_ <- function(bm) {
  lyrs <- bm$Z.bm[[1]]$layers
  lyrs <- lyrs[order(lyrs$t.l), ]
  if(nrow(lyrs) == 0) {
    return(data.frame(type = character(), t.l = numeric(), t.u = numeric(),
                      stringsAsFactors = FALSE))
  }
  types <- sapply(bm$Z.bm, function(Z.bm) {
    Z.bm$layers$type[match(lyrs$t.l, Z.bm$layers$t.l)]
  })
  types <- matrix(types, nrow = nrow(lyrs))
  type <- apply(types, 1, function(x) { paste(sort(unique(x[!is.na(x)])), collapse = "/") })
  data.frame(type = as.character(type), t.l = lyrs$t.l, t.u = lyrs$t.u,
             stringsAsFactors = FALSE)
}

# Index of the layer over the interval starting at each time t (NA if none),
# for layers sorted by t.l
interval_layer_ <- function(lyrs, t) {
  j <- findInterval(t, lyrs$t.l)
  j[j == 0] <- NA
  j[!is.na(j) & t >= lyrs$t.u[pmax(j, 1)]] <- NA
  j
}

# The parts of [start, end] which are not covered by any layer
layer_gaps_ <- function(t.l, t.u, start, end) {
  o <- order(t.l)
  from <- c(start, cummax(t.u[o]))
  to <- c(t.l[o], end)
  keep <- to > from
  data.frame(from = from[keep], to = to[keep])
}

# Text describing a covariance matrix, eg "unit variances, correlation 0.5"
cov_text_ <- function(cov, digits) {
  v <- diag(cov)
  if(all(v == 1)) {
    var.text <- "unit variances"
  } else if(all(v == v[1])) {
    var.text <- paste0("variances ", fmt_(v[1], digits))
  } else if(length(v) <= 4) {
    var.text <- paste0("variances ", paste(fmt_(v, digits), collapse = ", "))
  } else {
    var.text <- paste0("variances from ", fmt_(min(v), digits), " to ", fmt_(max(v), digits))
  }
  r <- stats::cov2cor(cov)[upper.tri(cov)]
  if(all(r == 0)) {
    return(paste0("independent, ", var.text))
  }
  if(length(r) == 1) {
    cor.text <- paste0("correlation ", fmt_(r, digits))
  } else if(all(r == r[1])) {
    cor.text <- paste0("all correlations ", fmt_(r[1], digits))
  } else {
    cor.text <- paste0("correlations from ", fmt_(min(r), digits), " to ", fmt_(max(r), digits))
  }
  paste0(var.text, ", ", cor.text)
}



# Text formatting helpers -----------------------------------------------------

fmt_ <- function(x, digits) {
  vapply(x, function(x) { format(x, digits = digits) }, character(1))
}

# Format a column of numbers to a common number of decimal places, leaving
# NA blank
format_num_ <- function(x, digits) {
  out <- rep("", length(x))
  ok <- !is.na(x)
  if(any(ok)) {
    out[ok] <- format(x[ok], digits = digits)
  }
  out
}

format_pct_ <- function(x) {
  paste0(format(signif(100 * x, 3)), "%")
}

count_text_ <- function(n, what, plural = paste0(what, "s")) {
  paste(n, if(n == 1) what else plural)
}

dim_text_ <- function(dim) {
  count_text_(dim, "dimension")
}

points_text_ <- function(t, digits) {
  if(length(t) == 1) {
    paste0("1 point at t = ", fmt_(t, digits))
  } else {
    paste0(count_text_(length(t), "point"), " on [", fmt_(min(t), digits), ", ", fmt_(max(t), digits), "]")
  }
}

layers_text_ <- function(lyrs, t) {
  if(nrow(lyrs) == 0) {
    return("no layers")
  }
  span <- max(t) - min(t)
  gaps <- layer_gaps_(lyrs$t.l, lyrs$t.u, min(t), max(t))
  covered <- (span - sum(gaps$to - gaps$from)) / span
  paste0(count_text_(nrow(lyrs), "layer"), " (", format_pct_(covered), " of time)")
}

sep_ <- function() {
  if(cli::is_utf8_output()) " \u00b7 " else ", "
}

vellipsis_ <- function() {
  if(cli::is_utf8_output()) "\u22ee" else ":"
}

# Colour layer types to match plot.BrownianMotion
style_layer_type_ <- function(type) {
  out <- ifelse(is.na(type), "", type)
  base <- sub("-bb$", "", out)
  out[base == "bessel"] <- cli::col_magenta(out[base == "bessel"])
  out[base == "intersection"] <- cli::col_blue(out[base == "intersection"])
  out[base == "localised"] <- cli::col_red(out[base == "localised"])
  as.character(out)
}

header_rule_ <- function(title, right, rule.width, width) {
  as.character(cli::rule(left = cli::style_bold(title), right = right,
                         width = min(rule.width, width)))
}

footer_ <- function(notes, width) {
  if(length(notes) == 0) {
    return(character())
  }
  as.character(cli::col_grey(wrap_(paste(cli::symbol$info, notes), width, indent = 2)))
}

# Wrap each of x to the width, with prefix before every line and continuation
# lines indented further
wrap_ <- function(x, width, indent = 0, prefix = "") {
  unlist(lapply(x, function(x) {
    lines <- as.character(cli::ansi_strwrap(x, width = width - nchar(prefix), exdent = indent))
    paste0(prefix, if(length(lines) == 0) "" else lines)
  }))
}

# Which of n.pts rows to show: all of them, or the first and last few
print_rows_ <- function(n.pts, n) {
  if(is.null(n)) {
    n <- if(n.pts <= 20) n.pts else 10
  }
  if(!checkmate::test_number(n, lower = 1)) {
    stop("n must be a number of at least 1")
  }
  if(n >= n.pts) {
    return(seq_len(n.pts))
  }
  head.n <- ceiling(n / 2)
  tail.n <- floor(n / 2)
  c(seq_len(head.n), seq.int(n.pts - tail.n + 1, length.out = tail.n))
}

# Lay out a table of character columns (named by their headers) for the given
# rows of n.pts, with row numbers and a marker where rows are omitted. Columns
# in left are left aligned, and others right aligned. If the table is wider than
# width, the labels column is truncated, and then columns in droppable are
# removed in turn.
format_table_ <- function(cols, rows, n.pts, width, left = character(), droppable = character()) {
  cols <- c(list(" " = as.character(rows)), cols)
  notes <- character()

  col_width <- function(nm) {
    max(cli::ansi_nchar(c(nm, cols[[nm]]), type = "width"))
  }
  table_width <- function() {
    sum(vapply(names(cols), col_width, numeric(1))) + 2 * (length(cols) - 1)
  }

  dropped <- character()
  repeat {
    over <- table_width() - width
    if(over <= 0) {
      break
    }
    if(!is.null(cols$labels)) {
      room <- col_width("labels") - over
      if(room >= 10) {
        cols$labels <- as.character(cli::ansi_strtrim(cols$labels, room))
        next
      }
    }
    if(length(droppable) > length(dropped)) {
      nm <- droppable[length(dropped) + 1]
      cols[[nm]] <- NULL
      dropped <- c(dropped, nm)
    } else if(!is.null(cols$labels)) {
      cols$labels <- NULL
      notes <- c(notes, "Labels not shown: use bm[, \"labels\"] to see them.")
    } else {
      break
    }
  }
  if(length(dropped) > 0) {
    notes <- c(notes, paste0(count_text_(length(dropped), "component"),
                             " not shown: use a wider console to see ",
                             if(length(dropped) == 1) "it." else "them."))
  }

  lines <- align_cols_(cols, left)

  # Mark omitted rows in the middle of the table
  gap <- which(diff(rows) > 1)
  if(length(gap) > 0) {
    n.more <- n.pts - length(rows)
    marker <- paste0(cli::ansi_align(vellipsis_(), col_width(" "), "right"), "  ",
                     cli::col_grey(paste0(count_text_(n.more, "more point"),
                                          " (use print(n = ...) to see more)")))
    lines <- append(lines, as.character(marker), after = gap + 1)
  }

  list(lines = lines, width = table_width(), notes = notes)
}

# Align character columns (named by their headers) into lines of text, with a
# bold header. The first column holds row names (named " "), which are shown in
# grey and have no header. Columns in left are left aligned, and others right
# aligned.
align_cols_ <- function(cols, left = character(), indent = "") {
  aligned <- lapply(names(cols), function(nm) {
    x <- c(if(nm == names(cols)[1]) "" else nm, cols[[nm]])
    width <- max(cli::ansi_nchar(x, type = "width"))
    as.character(cli::ansi_align(x, width, if(nm %in% left) "left" else "right"))
  })
  header <- vapply(aligned, `[`, character(1), 1)
  header <- style_padded_(header, cli::style_bold)
  body <- lapply(aligned, `[`, -1)
  body[[1]] <- style_padded_(body[[1]], cli::col_grey)
  lines <- c(paste(header, collapse = "  "),
             do.call("paste", c(body, sep = "  ")))
  # Left aligned final columns are padded, which leaves trailing spaces
  paste0(indent, sub(" +$", "", lines))
}

# Apply a style to the text of padded cells, leaving the padding unstyled so
# that trailing spaces can be removed
style_padded_ <- function(x, style) {
  text <- trimws(x)
  out <- x
  has <- text != ""
  out[has] <- paste0(sub("^( *).*$", "\\1", x[has]),
                     style(text[has]),
                     sub("^.*[^ ]( *)$", "\\1", x[has]))
  out
}
