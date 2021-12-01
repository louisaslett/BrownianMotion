#' Manually specify the value of a point in the skeleton
#'
#' Note: Only valid for time points that don't lie within existing layers. Cannot be an existing time within the skeleton (in this case you should delete then add the point).
#'
#' @param bm a Brownian motion object from which simulation should continue.
#'   Note the object is updated in place
#' @param t vector of fixed times of the Brownian motion.  Defaults to value of 0.
#' @param W_t matrix of size the length of t (rows) x dimension (column). If dim = 1, a vector of the length of t is acceptable. If dim > 1 and the length of t = 1, then a vector of length 1 is acceptable (corresponding to the same value at every co-ordinate), or length d is acceptable (corresponding to the d components). Defaults to 0 for every component.
#' @param label vector of length 1 or the length of t indicating the introduction of user specified labels (if any). By default, uses any names in the vector of times.
#'
#' @export
add.points <- function(bm, t, W_t, label = names(t)) {
  UseMethod("add.points")
}

#' @export
add.points.BrownianMotionNd <- function(bm, t, W_t, label = names(t)) {
  if(!("BrownianMotionNd" %in% class(bm))) {
    stop("bm argument must be a BrownianMotionNd object.")
  }

  # t checks
  if(!isTRUE(e <- checkmate::check_vector(t, strict = TRUE, any.missing = FALSE, min.len = 1, unique = TRUE, null.ok = FALSE))) {
    stop("Error with argument t: ", e)
  }

  # W_t checks
  if(bm$dim == 1L) {
    if(!isTRUE(e <- checkmate::check_vector(W_t, strict = TRUE, any.missing = FALSE, len = length(t), null.ok = FALSE))) {
      stop("Error with argument W_t: ", e)
    }
    if(!is.matrix(W_t)) {
      W_t <- matrix(W_t, ncol = 1)
    }
  } else if(bm$dim > 1L) {
    if(length(t) == 1) { # allow W_t to be scalar or vector
      if(is.matrix(W_t) && nrow(W_t) == 1) {
        W_t <- W_t[1,]
      }
      if(!isTRUE(e <- checkmate::check_vector(W_t, strict = TRUE, any.missing = FALSE, null.ok = FALSE))) {
        stop("Error with argument W_t: ", e)
      }
      if(length(W_t) == 1) {
        W_t <- rep(W_t, bm$dim)
      } else if(length(W_t) != bm$dim) {
        stop("Error with argument W_t: must be of length 1 or of length dim")
      }
      W_t <- matrix(W_t, nrow = 1)
    } else { # W_t must be matrix
      if(!isTRUE(e <- checkmate::check_matrix(W_t, any.missing = FALSE, nrows = length(t), ncols = bm$dim))) {
        stop("Error with argument W_t: ", e)
      }
    }
  }

  for(d in 1:bm$dim) {
    add.points(bm$Z.bm[[d]], t, W_t[,d])
  }

  invisible(bm)
}

#' @export
add.points.BrownianMotion <- function(bm, t, W_t, label = names(t)) {
  if(!("BrownianMotion" %in% class(bm))) {
    stop("bm argument must be a BrownianMotion object.")
  }
  if(!is.realvector(t)) {
    stop("t must be a vector of times.")
  }
  if(!is.realvector(W_t)) {
    stop("W_t must be a vector of observations.")
  }
  if(length(t) != length(W_t)) {
    stop("t and W_t must be the same length.")
  }

  assert.bmlabel(label, t)
  if(!is.null(label) && length(label) == 1) {
    label <- rep(label, length(t))
  }

  if(any(t %in% bm$t)) {
    stop("t is already a time in the skeleton. No times added.")
  }
  if(any(sapply(t, function(x) { any(x >= bm$layers$t.l | x <= bm$layers$t.u) }))) {
    stop("t is within a layer (currently not supported, but may be added in future release). No times added.")
  }

  if(length(t) == 1) { # for single value insertion sort for speed
    if(t < min(bm$t)) { # insertion at the start
      bm$t <- c(t,
                bm$t)
      bm$W_t <- c(W_t,
                  bm$W_t)
      bm$W_tm <- c(W_t,
                   bm$W_tm)
    } else if(t > max(bm$t)) {
      bm$t <- c(bm$t,
                t)
      bm$W_t <- c(bm$W_t,
                  W_t)
      bm$W_tm <- c(bm$W_tm,
                   W_t)
    } else {
      i <- max(which(bm$t < t))
      bm$t <- c(bm$t[1:i],
                t,
                bm$t[(i+1):length(bm$t)])
      bm$W_t <- c(bm$W_t[1:i],
                  W_t,
                  bm$W_t[(i+1):length(bm$W_t)])
      bm$W_tm <- c(bm$W_tm[1:i],
                   W_t,
                   bm$W_tm[(i+1):length(bm$W_tm)])
    }
  } else { # otherwise resort to append and sort
    bm$t <- c(bm$t, t)
    bm$W_t <- c(bm$W_t, W_t)
    bm$W_tm <- c(bm$W_tm, W_t)
    o <- order(bm$t)
    bm$W_t <- bm$W_t[o]
    bm$W_tm <- bm$W_tm[o]
    bm$t <- bm$t[o]
  }

  # Labels
  add.labels_(bm, "user", t)
  add.labels_(bm, "forced", t)
  add.labels_(bm, label, t)
  if(any(bm$t < bm$labels[["start"]] | bm$t > bm$labels[["end"]])) {
    bm$labels[["start"]] <- min(bm$t)
    bm$labels[["end"]] <- max(bm$t)
    bm$labels[["seg.start"]] <- c(bm$labels[["start"]], bm$t[which(bm$W_t != bm$W_tm)])
    bm$labels[["seg.end"]] <- c(bm$t[which(bm$W_t != bm$W_tm)], bm$labels[["end"]])
  }

  invisible(bm)
}
