#' Extract parts of a Brownian motion object
#'
#' Brownian motion objects are indexed like a data frame, \code{bm[i, var]},
#' where \code{i} selects skeleton points and \code{var} selects what to return.
#'
#' Points can be selected in three ways:
#' \itemize{
#'   \item by position, e.g. \code{bm[1, "W_t"]} for the first point;
#'   \item by time, using the named argument \code{t}, e.g.
#'     \code{bm[t = 1.5, "W_t"]}. A range \code{t = a:b} selects \code{a},
#'     \code{b} and every skeleton time in between (it is not a sequence of
#'     integers). Times not in the skeleton give \code{NA} values, but the
#'     layers containing them are still returned;
#'   \item by label, e.g. \code{bm["user", "t"]} (see \code{\link{add_labels}}).
#' }
#' Leaving \code{i} empty, as in \code{bm[, "W_t"]}, selects every point.
#'
#' \code{var} is one of \code{"t"}, \code{"W_t"} (the path), \code{"W_tm"}
#' (the left limit of the path, which differs from \code{W_t} only at jumps
#' created by \code{\link{add_segment}} or \code{\link{concat_bm}}),
#' \code{"layers"} or \code{"labels"}. If it is omitted, a list of all of these
#' is returned. For the transformed view \code{bm$Z} of a multi-dimensional
#' path, use \code{"Z_t"} and \code{"Z_tm"} in place of \code{"W_t"} and
#' \code{"W_tm"}.
#'
#' @param bm a Brownian motion object.
#' @param i the points to select, as indices or labels. Leave empty to select
#'   every point.
#' @param var what to return (see Details). Omit to return a list of
#'   everything.
#' @param t times to select, in place of \code{i}; must be given by name, as
#'   in \code{bm[t = 1.5, "W_t"]}.
#'
#' @return Depending on \code{var}: a numeric vector (\code{"t"}, and
#'   \code{"W_t"}/\code{"W_tm"} in one dimension); a matrix with one column per
#'   dimension (\code{"W_t"}/\code{"W_tm"} in several dimensions); a tibble of
#'   layers (\code{"layers"}); a named list of label times (\code{"labels"});
#'   or a list of all of these if \code{var} is omitted.
#'
#' @examples
#' set.seed(1)
#' bm <- create_bm(t = c(-1, 0, 1), W_t = c(100, 101, 99))
#' bm |> sim(c(0.5, 2)) |> layers(0, 1, type = "bessel")
#'
#' bm[, "t"]
#' bm[1, "W_t"]           # first point
#' bm[t = 0.5, "W_t"]     # by time
#' bm[t = 0:1, "W_t"]     # 0, 1 and all skeleton times in between
#' bm["user", "t"]        # by label
#' bm[, "labels"]
#' bm[t = 0.25, "layers"] # layer containing a time not in the skeleton
#' bm[2]                  # everything about the second point
#' what <- "W_tm"
#' bm[t = 0:1, what]      # var can be held in a variable
#'
#' # Multi-dimensional paths, and their transformed view
#' bm2 <- create_bm(dim = 2, cov = matrix(c(1, 0.5, 0.5, 1), 2))
#' bm2 |> sim(1:2)
#' bm2[, "W_t"]
#' bm2$Z[, "Z_t"]
#' @aliases BrownianMotion-indexing
#' @rdname indexing
#' @export
`[.BrownianMotion` <- function(bm, i, var, t) {
  # var is read from the call below; i and t are optional
  x <- if(missing(i)) NULL else i
  if(missing(t)) t <- NULL
  # Determine what variable has been requested (eg bm[1,"W_t"])
  # or is NA if not specified (eg bm[1])
  var <- index_var_(sys.call(), parent.frame(), c("t", "W_t", "W_tm", "layers", "labels"))

  if(is.null(x)) {
    x <- 1:length(bm$t)
  }

  if(!is.null(t)) {
    t2 <- substitute(t)
    env <- parent.frame()
    if(as.character(t2)[1] == ":") {
      t <- c(eval(t2[[2]], envir = env),
             bm$t[bm$t > eval(t2[[2]], envir = env) &
                  bm$t < eval(t2[[3]], envir = env)],
             eval(t2[[3]], envir = env))
    } else {
      t <- eval(t, envir = env)
    }
  }

  # Do we have label in x, and maybe variable to retrieve in t?
  # If so, get time for each label (and dump variable from t, as we now have in var)
  if(!is.numeric(t) && is.character(x)) {
    x <- x[x %in% names(bm$labels)]
    t <- c(bm$labels[x], recursive = TRUE)
  }

  # Create dummy to hold times that are not in the skeleton for which
  # we will only return the relevant layer info (if any)
  missing <- NULL

  # Do we have a time in t?
  # If so, find index and put in x,
  # otherwise, x must be an index already
  if(is.numeric(t)) {
    t <- unique(t)
    reorder.t <- order(order(t))
    missing <- t[!(t %in% bm$t)]
    x <- which(bm$t %in% t)
  } else {
    reorder.t <- order(order(x))
    t <- bm$t[x]
  }

  # Deal with the case where none of the times are skeleton points
  if(length(x) == 0) {
    if(is.na(var)) {
      # Get everything
      if(nrow(bm$layers) == 0) {
        lyrs.x <- NULL
      } else {
        # 20/11/23: If there is exactly one layer, the next line is a bug, fixed (I think) by the line after
        #lyrs.x <- (rowSums(sapply(c(bm$t[x], missing), function(t) { bm$layers$t.l <= t & bm$layers$t.u > t })) > 0)
        lyrs.x <- (rowSums(simplify2array(lapply(c(bm$t[x], missing), function(t) { bm$layers$t.l <= t & bm$layers$t.u > t }), except = 0)) > 0)
      }
      return(list(t = t,
                  W_t = rep(NA, length(t)),
                  W_tm = rep(NA, length(t)),
                  layers = bm$layers[lyrs.x,],
                  labels = list()))
    } else if(var == "t") {
      # Get just t
      return(unname(t))
    } else if(var == "W_t") {
      # Get just W_t
      return(rep(NA, length(t)))
    } else if(var == "W_tm") {
      # Get just W_tm
      return(rep(NA, length(t)))
    } else if(var == "layers") {
      # Get just layers
      if(nrow(bm$layers) == 0) {
        lyrs.x <- NULL
      } else {
        # 20/11/23: If there is exactly one layer, the next line is a bug, fixed (I think) by the line after
        #lyrs.x <- (rowSums(sapply(c(bm$t[x], missing), function(t) { bm$layers$t.l <= t & bm$layers$t.u > t })) > 0)
        lyrs.x <- (rowSums(simplify2array(lapply(c(bm$t[x], missing), function(t) { bm$layers$t.l <= t & bm$layers$t.u > t }), except = 0)) > 0)
      }
      return(bm$layers[lyrs.x,])
    } else if(var == "labels") {
      return(list())
    }
  }

  if(length(x)>0 && any(as.integer(x)!=x)) {
    stop("must be an integer index")
  }
  if(length(x)>0 && any(x < 1 | x > length(bm$t))) {
    stop("index out of bounds")
  }

  if(is.na(var) || var == "") {
    # Get everything

    if(nrow(bm$layers) == 0) {
      lyrs.x <- NULL
    } else {
      # Extract layers at both skeleton times and non-skeleton times
      # 20/11/23: If there is exactly one layer, the next line is a bug, fixed (I think) by the line after
      #lyrs.x <- (rowSums(as.matrix(sapply(c(bm$t[x], missing), function(t) { bm$layers$t.l <= t & bm$layers$t.u > t }))) > 0)
      lyrs.x <- (rowSums(simplify2array(lapply(c(bm$t[x], missing), function(t) { bm$layers$t.l <= t & bm$layers$t.u > t }), except = 0)) > 0)
    }

    res <- matrix(c(unname(bm$t[x]), missing,
                    unname(bm$W_t[x]), rep(NA, length(missing)),
                    unname(bm$W_tm[x]), rep(NA, length(missing))), ncol = 3)
    res <- res[order(res[,1]),,drop=FALSE][reorder.t,,drop=FALSE]
    res.labels <- lapply(bm$labels, function(ts) { ts[ts %in% res[,1]] })

    return(list(t = res[,1],
                W_t = res[,2],
                W_tm = res[,3],
                layers = bm$layers[lyrs.x,],
                labels = res.labels[lapply(res.labels, length)>0]))
  } else if(var == "t") {
    # Get just t
    return(unname(t))
  } else if(var == "W_t") {
    # Get just W_t
    res <- matrix(c(unname(bm$t[x]), missing,
                    unname(bm$W_t[x]), rep(NA, length(missing))), ncol = 2)
    res <- res[order(res[,1]),,drop=FALSE][reorder.t,,drop=FALSE]
    return(unname(res[,2]))
  } else if(var == "W_tm") {
    # Get just W_tm
    res <- matrix(c(unname(bm$t[x]), missing,
                    unname(bm$W_tm[x]), rep(NA, length(missing))), ncol = 2)
    res <- res[order(res[,1]),,drop=FALSE][reorder.t,,drop=FALSE]
    return(unname(res[,2]))
  } else if(var == "layers") {
    # Get just layers
    if(nrow(bm$layers) == 0) {
      lyrs.x <- NULL
    } else {
      # 20/11/23: If there is exactly one layer, the next line is a bug, fixed (I think) by the line after
      #lyrs.x <- (rowSums(sapply(c(bm$t[x], missing), function(t) { bm$layers$t.l <= t & bm$layers$t.u > t })) > 0)
      lyrs.x <- (rowSums(simplify2array(lapply(c(bm$t[x], missing), function(t) { bm$layers$t.l <= t & bm$layers$t.u > t }), except = 0)) > 0)
    }
    return(bm$layers[lyrs.x,])
  } else if(var == "labels") {
    res <- matrix(c(unname(bm$t[x]), missing,
                    unname(bm$W_t[x]), rep(NA, length(missing))), ncol = 2)

    res.labels <- lapply(bm$labels, function(ts) { ts[ts %in% res[,1]] })

    return(res.labels[lapply(res.labels, length)>0])
  }
}



#' @rdname indexing
#' @export
`[.BrownianMotionNd` <- function(bm, i, var, t) {
  # var is read from the call below; i and t are optional
  x <- if(missing(i)) NULL else i
  if(missing(t)) t <- NULL
  # Determine what variable has been requested (eg bm[1,"W_t"])
  # or is NA if not specified (eg bm[1])
  var <- index_var_(sys.call(), parent.frame(), c("t", "W_t", "W_tm", "layers", "labels"))

  if(is.null(x)) {
    x <- 1:length(bm$Z.bm[[1]]$t)
  }

  if(!is.null(t)) {
    t2 <- substitute(t)
    env <- parent.frame()
    if(as.character(t2)[1] == ":") {
      t <- c(eval(t2[[2]], envir = env),
             bm$Z.bm[[1]]$t[bm$Z.bm[[1]]$t > eval(t2[[2]], envir = env) &
                            bm$Z.bm[[1]]$t < eval(t2[[3]], envir = env)],
             eval(t2[[3]], envir = env))
    } else {
      t <- eval(t, envir = env)
    }
  }

  # Do we have label in x, and maybe variable to retrieve in t?
  # If so, get time for each label (and dump variable from t, as we now have in var)
  if(!is.numeric(t) && is.character(x)) {
    # TODO: improve this, we're assuming first dimension labels are same as all
    #       other dimensions. May not be true for internally generated labels
    #       from layer routines
    x <- x[x %in% names(bm$Z.bm[[1]]$labels)]
    t <- c(bm$Z.bm[[1]]$labels[x], recursive = TRUE)
  }

  # Create dummy to hold times that are not in the skeleton for which
  # we will only return the relevant layer info (if any)
  missing <- NULL

  # Do we have a time in t?
  # If so, find index and put in x,
  # otherwise, x must be an index already
  if(is.numeric(t)) {
    t <- unique(t)
    reorder.t <- order(order(t))
    # TODO: improve this? This gives correct result, but ...
    #       We've got $t duplicated in all Z.bm objects, but this
    #       makes implementation easier (especially for layers where we fill in
    #       anything missing across dims after simulation)
    missing <- t[!(t %in% bm$Z.bm[[1]]$t)]
    x <- which(bm$Z.bm[[1]]$t %in% t)
  } else {
    reorder.t <- order(order(x))
    t <- bm$Z.bm[[1]]$t[x]
  }

  # Deal with the case where none of the times are skeleton points
  if(length(x) == 0) {
    if(is.na(var)) {
      if(nrow(bm$Z.bm[[1]]$layers) == 0) {
        # TODO: Murray says he'd love to spend some time allowing
        #       layers in only certain dimensions and for the
        #       path to be unconstrained in others.
        #       Louis places face in palms
        # TODO: Louis thinks it would be lovely to have a rainbow of
        #       colours for different time segments and overlapping transparent
        #       polygons and infilled paths.
        #       Murray places face in palms
        lyrs <- tibble(
          t.l = numeric(),
          t.u = numeric(),
          outer.L = matrix(numeric(), nrow = 0, ncol = bm$dim),
          outer.U = matrix(numeric(), nrow = 0, ncol = bm$dim),
          inner.cube = list()
        )
      } else {
        lyrs <- dplyr::bind_rows(lapply(c(bm$Z.bm[[1]]$t[x], missing), function(t) { transformlayersNd_(bm, t) }))
      }
      return(list(t = t,
                  W_t = rep(NA, length(t)),
                  W_tm = rep(NA, length(t)),
                  layers = lyrs,
                  labels = list()))
    } else if(var == "t") {
      # Get just t
      return(unname(t))
    } else if(var == "W_t") {
      # Get just W_t
      return(matrix(NA, nrow = length(t), ncol = bm$dim))
    } else if(var == "W_tm") {
      # Get just W_tm
      return(matrix(NA, nrow = length(t), ncol = bm$dim))
    } else if(var == "layers") {
      # Get just layers
      if(nrow(bm$Z.bm[[1]]$layers) == 0) {
        lyrs <- tibble(
          t.l = numeric(),
          t.u = numeric(),
          outer.L = matrix(numeric(), nrow = 0, ncol = bm$dim),
          outer.U = matrix(numeric(), nrow = 0, ncol = bm$dim),
          inner.cube = list()
        )
      } else {
        lyrs <- dplyr::bind_rows(lapply(c(bm$Z.bm[[1]]$t[x], missing), function(t) { transformlayersNd_(bm, t) }))
      }
      return(lyrs)
    } else if(var == "labels") {
      return(list())
    }
  }

  # index in x, maybe var in t
  if(length(x)>0 && any(as.integer(x)!=x)) {
    stop("must be an integer index")
  }
  if(length(x)>0 && any(x < 1 | x > length(bm$Z.bm[[1]]$t))) {
    stop("index out of bounds")
  }

  if(is.na(var) || var == "") {
    # Get everything

    if(nrow(bm$Z.bm[[1]]$layers) == 0) {
      lyrs <- tibble(
        t.l = numeric(),
        t.u = numeric(),
        outer.L = matrix(numeric(), nrow = 0, ncol = bm$dim),
        outer.U = matrix(numeric(), nrow = 0, ncol = bm$dim),
        inner.cube = list()
      )
    } else {
      lyrs <- dplyr::bind_rows(lapply(c(bm$Z.bm[[1]]$t[x], missing), function(t) { transformlayersNd_(bm, t) }))
    }

    # res <- matrix(c(unname(bm$t[x]), missing,
    #                 unname(bm$W_t[x]), rep(NA, length(missing)),
    #                 unname(bm$W_tm[x]), rep(NA, length(missing))), ncol = 3)
    # res <- res[order(res[,1]),,drop=FALSE][reorder.t,,drop=FALSE]

    res.t <- c(unname(bm$Z.bm[[1]]$t[x]), missing)
    res.W_t <- rbind(getWtNd_(bm, x), matrix(NA, nrow = length(missing), ncol = bm$dim))
    res.W_tm <- rbind(getWtmNd_(bm, x), matrix(NA, nrow = length(missing), ncol = bm$dim))
    res.labels <- lapply(bm$Z.bm[[1]]$labels, function(ts) { ts[ts %in% res.t] })

    o <- order(res.t)
    return(list(t = res.t[o][reorder.t],
                W_t = res.W_t[o,,drop=FALSE][reorder.t,,drop=FALSE],
                W_tm = res.W_tm[o,,drop=FALSE][reorder.t,,drop=FALSE],
                layers = lyrs,
                labels = res.labels[lapply(res.labels, length)>0]))
  } else if(var == "t") {
    # Get just t
    return(unname(t))
  } else if(var == "W_t") {
    # Get just W_t
    res.t <- c(unname(bm$Z.bm[[1]]$t[x]), missing)
    res.W_t <- rbind(getWtNd_(bm, x), matrix(NA, nrow = length(missing), ncol = bm$dim))

    o <- order(res.t)
    return(res.W_t[o,,drop=FALSE][reorder.t,,drop=FALSE])
  } else if(var == "W_tm") {
    # Get just W_tm
    res.t <- c(unname(bm$Z.bm[[1]]$t[x]), missing)
    res.W_tm <- rbind(getWtmNd_(bm, x), matrix(NA, nrow = length(missing), ncol = bm$dim))

    o <- order(res.t)
    return(res.W_tm[o,,drop=FALSE][reorder.t,,drop=FALSE])
  } else if(var == "layers") {
    # Get just layers
    if(nrow(bm$Z.bm[[1]]$layers) == 0) {
      lyrs <- tibble(
        t.l = numeric(),
        t.u = numeric(),
        outer.L = matrix(numeric(), nrow = 0, ncol = bm$dim),
        outer.U = matrix(numeric(), nrow = 0, ncol = bm$dim),
        inner.cube = list()
      )
    } else {
      lyrs <- dplyr::bind_rows(lapply(c(bm$Z.bm[[1]]$t[x], missing), function(t) { transformlayersNd_(bm, t) }))
    }
    return(lyrs)
  } else if(var == "labels") {
    res.t <- c(unname(bm$Z.bm[[1]]$t[x]), missing)
    res.labels <- lapply(bm$Z.bm[[1]]$labels, function(ts) { ts[ts %in% res.t] })

    return(res.labels[lapply(res.labels, length)>0])
  }
}



#' @rdname indexing
#' @export
`[.BrownianMotionNdZ` <- function(bm, i, var, t) {
  # var is read from the call below; i and t are optional
  x <- if(missing(i)) NULL else i
  if(missing(t)) t <- NULL
  bm <- bm[[1]]

  # Determine what variable has been requested (eg bm[1,"W_t"])
  # or is NA if not specified (eg bm[1])
  var <- index_var_(sys.call(), parent.frame(), c("t", "Z_t", "Z_tm", "layers", "labels"))

  if(is.null(x)) {
    x <- 1:length(bm$Z.bm[[1]]$t)
  }

  if(!is.null(t)) {
    t2 <- substitute(t)
    env <- parent.frame()
    if(as.character(t2)[1] == ":") {
      t <- c(eval(t2[[2]], envir = env),
             bm$Z.bm[[1]]$t[bm$Z.bm[[1]]$t > eval(t2[[2]], envir = env) &
                              bm$Z.bm[[1]]$t < eval(t2[[3]], envir = env)],
             eval(t2[[3]], envir = env))
    } else {
      t <- eval(t, envir = env)
    }
  }

  # Do we have label in x, and maybe variable to retrieve in t?
  # If so, get time for each label (and dump variable from t, as we now have in var)
  if(!is.numeric(t) && is.character(x)) {
    # TODO: improve this, we're assuming first dimension labels are same as all
    #       other dimensions. May not be true for internally generated labels
    #       from layer routines
    x <- x[x %in% names(bm$Z.bm[[1]]$labels)]
    t <- c(bm$Z.bm[[1]]$labels[x], recursive = TRUE)
  }

  # Create dummy to hold times that are not in the skeleton for which
  # we will only return the relevant layer info (if any)
  missing <- NULL

  # Do we have a time in t?
  # If so, find index and put in x,
  # otherwise, x must be an index already
  if(is.numeric(t)) {
    t <- unique(t)
    reorder.t <- order(order(t))
    # TODO: improve this? This gives correct result, but ...
    #       We've got $t duplicated in all Z.bm objects, but this
    #       makes implementation easier (especially for layers where we fill in
    #       anything missing across dims after simulation)
    missing <- t[!(t %in% bm$Z.bm[[1]]$t)]
    x <- which(bm$Z.bm[[1]]$t %in% t)
  } else {
    reorder.t <- order(order(x))
    t <- bm$Z.bm[[1]]$t[x]
  }

  # Deal with the case where none of the times are skeleton points
  if(length(x) == 0) {
    if(is.na(var)) {
      if(nrow(bm$Z.bm[[1]]$layers) == 0) {
        # TODO: Murray says he'd love to spend some time allowing
        #       layers in only certain dimensions and for the
        #       path to be unconstrained in others.
        #       Louis places face in palms
        # TODO: Louis thinks it would be lovely to have a rainbow of
        #       colours for different time segments and overlapping transparent
        #       polygons and infilled paths.
        #       Murray places face in palms
        lyrs <- tibble(
          t.l = numeric(),
          t.u = numeric(),
          outer.L = matrix(numeric(), nrow = 0, ncol = bm$dim),
          outer.U = matrix(numeric(), nrow = 0, ncol = bm$dim),
          inner.cube = list()
        )
      } else {
        lyrs <- dplyr::bind_rows(lapply(c(bm$Z.bm[[1]]$t[x], missing), function(t) { transformlayersNd_(bm, t, do.transform = FALSE) }))
      }
      return(list(t = t,
                  Z_t = rep(NA, length(t)),
                  Z_tm = rep(NA, length(t)),
                  layers = lyrs,
                  labels = list()))
    } else if(var == "t") {
      # Get just t
      return(unname(t))
    } else if(var == "Z_t") {
      # Get just Z_t
      return(matrix(NA, nrow = length(t), ncol = bm$dim))
    } else if(var == "Z_tm") {
      # Get just Z_tm
      return(matrix(NA, nrow = length(t), ncol = bm$dim))
    } else if(var == "layers") {
      # Get just layers
      if(nrow(bm$Z.bm[[1]]$layers) == 0) {
        lyrs <- tibble(
          t.l = numeric(),
          t.u = numeric(),
          outer.L = matrix(numeric(), nrow = 0, ncol = bm$dim),
          outer.U = matrix(numeric(), nrow = 0, ncol = bm$dim),
          inner.cube = list()
        )
      } else {
        lyrs <- dplyr::bind_rows(lapply(c(bm$Z.bm[[1]]$t[x], missing), function(t) { transformlayersNd_(bm, t, do.transform = FALSE) }))
      }
      return(lyrs)
    } else if(var == "labels") {
      return(list())
    }
  }

  # index in x, maybe var in t
  if(length(x)>0 && any(as.integer(x)!=x)) {
    stop("must be an integer index")
  }
  if(length(x)>0 && any(x < 1 | x > length(bm$Z.bm[[1]]$t))) {
    stop("index out of bounds")
  }

  if(is.na(var) || var == "") {
    # Get everything

    if(nrow(bm$Z.bm[[1]]$layers) == 0) {
      lyrs <- tibble(
        t.l = numeric(),
        t.u = numeric(),
        outer.L = matrix(numeric(), nrow = 0, ncol = bm$dim),
        outer.U = matrix(numeric(), nrow = 0, ncol = bm$dim),
        inner.cube = list()
      )
    } else {
      lyrs <- dplyr::bind_rows(lapply(c(bm$Z.bm[[1]]$t[x], missing), function(t) { transformlayersNd_(bm, t, do.transform = FALSE) }))
    }

    res.t <- c(unname(bm$Z.bm[[1]]$t[x]), missing)
    res.Z_t <- rbind(getZtNd_(bm, x), matrix(NA, nrow = length(missing), ncol = bm$dim))
    res.Z_tm <- rbind(getZtmNd_(bm, x), matrix(NA, nrow = length(missing), ncol = bm$dim))
    res.labels <- lapply(bm$Z.bm[[1]]$labels, function(ts) { ts[ts %in% res.t] })

    o <- order(res.t)
    return(list(t = res.t[o][reorder.t],
                Z_t = res.Z_t[o,,drop=FALSE][reorder.t,,drop=FALSE],
                Z_tm = res.Z_tm[o,,drop=FALSE][reorder.t,,drop=FALSE],
                layers = lyrs,
                labels = res.labels[lapply(res.labels, length)>0]))
  } else if(var == "t") {
    # Get just t
    return(unname(t))
  } else if(var == "Z_t") {
    # Get just Z_t
    res.t <- c(unname(bm$Z.bm[[1]]$t[x]), missing)
    res.Z_t <- rbind(getZtNd_(bm, x), matrix(NA, nrow = length(missing), ncol = bm$dim))

    o <- order(res.t)
    return(res.Z_t[o,,drop=FALSE][reorder.t,,drop=FALSE])
  } else if(var == "Z_tm") {
    # Get just Z_tm
    res.t <- c(unname(bm$Z.bm[[1]]$t[x]), missing)
    res.Z_tm <- rbind(getZtmNd_(bm, x), matrix(NA, nrow = length(missing), ncol = bm$dim))

    o <- order(res.t)
    return(res.Z_tm[o,,drop=FALSE][reorder.t,,drop=FALSE])
  } else if(var == "layers") {
    # Get just layers
    if(nrow(bm$Z.bm[[1]]$layers) == 0) {
      lyrs <- tibble(
        t.l = numeric(),
        t.u = numeric(),
        outer.L = matrix(numeric(), nrow = 0, ncol = bm$dim),
        outer.U = matrix(numeric(), nrow = 0, ncol = bm$dim),
        inner.cube = list()
      )
    } else {
      lyrs <- dplyr::bind_rows(lapply(c(bm$Z.bm[[1]]$t[x], missing), function(t) { transformlayersNd_(bm, t, do.transform = FALSE) }))
    }
    return(lyrs)
  } else if(var == "labels") {
    res.t <- c(unname(bm$Z.bm[[1]]$t[x]), missing)
    res.labels <- lapply(bm$Z.bm[[1]]$labels, function(ts) { ts[ts %in% res.t] })

    return(res.labels[lapply(res.labels, length)>0])
  }
}
