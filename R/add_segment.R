#' Create / move a segment
#'
#' @param bm a Brownian motion object from which simulation should continue.
#'   Note the object is updated in place
#' @param l closed left end of interval to add a segment.  By default, `-Inf` which
#'     results in the introduction of a segment by shifting the entire path up to `r`
#' @param r open right end of interval to add a segment.  By default, `Inf` which
#'     results in the introduction of a segment by shifting the entire path from `l`
#' @param W_t vector of length the dimension of the Brownian motion specifying the coordinate of the Brownian motion at the time l (specified by the user). A vector of length 1 is acceptable (corresponding to the same value at every co-ordinate). Defaults to NULL, meaning this argument is ignored in favor of delta.
#' @param delta vector of length the dimension of the Brownian motion specifying the shift of the Brownian motion at the time l (specified by the user). A vector of length 1 is acceptable (corresponding to the same shift for every co-ordinate). Defaults to NULL, meaning this argument is ignored in favor of W_t.
#' @param label vector of length 1 or the length of t indicating the introduction of user specified labels (if any). By default, uses any names in the vector of times.
#'
#'
#' @export
add.segment <- function(bm, l = -Inf, r = Inf, W_t = NULL, delta = NULL, label = NULL) {
  UseMethod("add.segment")
}

#' @export
add.segment.BrownianMotionNd <- function(bm, l = -Inf, r = Inf, W_t = NULL, delta = NULL, label = NULL) {
  if(!("BrownianMotionNd" %in% class(bm))) {
    stop("bm argument must be a BrownianMotionNd object.")
  }

  if(!is.realscalar(l)) {
    stop("l must be a scalar.")
  }
  if(!is.realscalar(r)) {
    stop("r must be a scalar.")
  }
  if(l >= r) {
    stop("l must be less than r.")
  }
  if(l!=-Inf && !(l %in% bm$Z.bm[[1]]$t)) {
    stop("l must be a sampled path point")
  }
  if(r!=Inf && !(r %in% bm$Z.bm[[1]]$t)) {
    stop("r must be a sampled path point")
  }
  if(is.null(delta) && is.null(W_t)) {
    stop("One of delta or W_t must be specified")
  }
  if(!is.null(delta) && !is.null(W_t)) {
    stop("Can't specify both delta and W_t, please choose only one")
  }

  if(!is.null(W_t)) {
    # W_t checks
    if(bm$dim == 1L) {
      if(!isTRUE(e <- checkmate::check_vector(W_t, strict = TRUE, any.missing = FALSE, len = 1, null.ok = FALSE))) {
        stop("Error with argument W_t: ", e)
      }
      if(!is.matrix(W_t)) {
        W_t <- matrix(W_t, ncol = 1)
      }
    } else if(bm$dim > 1L) {
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
    }
  }
  if(!is.null(delta)) {
    # delta checks
    if(bm$dim == 1L) {
      if(!isTRUE(e <- checkmate::check_vector(delta, strict = TRUE, any.missing = FALSE, len = 1, null.ok = FALSE))) {
        stop("Error with argument delta: ", e)
      }
      if(!is.matrix(delta)) {
        delta <- matrix(delta, ncol = 1)
      }
    } else if(bm$dim > 1L) {
      if(is.matrix(delta) && nrow(delta) == 1) {
        delta <- delta[1,]
      }
      if(!isTRUE(e <- checkmate::check_vector(delta, strict = TRUE, any.missing = FALSE, null.ok = FALSE))) {
        stop("Error with argument delta: ", e)
      }
      if(length(delta) == 1) {
        delta <- rep(delta, bm$dim)
      } else if(length(delta) != bm$dim) {
        stop("Error with argument delta: must be of length 1 or of length dim")
      }
      delta <- matrix(delta, nrow = 1)
    }
  }

  for(d in 1:bm$dim) {
    add.segment(bm$Z.bm[[d]], l, r, W_t[,d], delta[,d], label)
  }

  invisible(bm)
}

#' @export
add.segment.BrownianMotion <- function(bm, l = -Inf, r = Inf, W_t = NULL, delta = NULL, label = NULL) {
  if(!("BrownianMotion" %in% class(bm))) {
    stop("bm argument must be a BrownianMotion object.")
  }
  if(!is.realscalar(l)) {
    stop("l must be a scalar.")
  }
  if(!is.realscalar(r)) {
    stop("r must be a scalar.")
  }
  if(l >= r) {
    stop("l must be less than r.")
  }
  if(l!=-Inf && !(l %in% bm$t)) {
    stop("l must be a sampled path point")
  }
  if(r!=Inf && !(r %in% bm$t)) {
    stop("r must be a sampled path point")
  }
  if(is.null(delta) && is.null(W_t)) {
    stop("One of delta or W_t must be specified")
  }
  if(!is.null(delta) && !is.null(W_t)) {
    stop("Can't specify both delta and W_t, please choose only one")
  }

  assert.bmlabel(label, c(l,r))
  if(!is.null(label) && length(label) == 1) {
    label <- rep(label, 2)
  }

  if(is.null(delta)) {
    if(l == -Inf) {
      delta <- W_t - bm$W_t[which.min(bm$t)]
    } else {
      delta <- W_t - bm$W_t[bm$t==l]
    }
  }
  # From here we can act only in terms of delta

  if(l == -Inf && r == Inf) {
    # Shift everything!
    bm$W_t <- bm$W_t + delta
    bm$W_tm <- bm$W_tm + delta
    bm$layers[,4:7] <- bm$layers[,4:7] + delta
    bm$user.layers[,3:4] <- bm$user.layers[,3:4] + delta
    bm$bb.local$W_t <- bm$bb.local$W_t + delta
    bm$bb.local$W_tm <- bm$bb.local$W_tm + delta
    bm$bb.local$layers[,4:7] <- bm$bb.local$layers[,4:7] + delta
    bm$bb.local$user.layers[,3:4] <- bm$bb.local$user.layers[,3:4] + delta
    add.labels_(bm, label, range(bm$t))
  } else {
    idx <- min(which(bm$t >= l)):max(which(bm$t <= r))
    if(r == Inf) {
      idx2 <- idx
    } else {
      idx2 <- head(idx,-1)
    }
    bm$W_t[idx2]  <- bm$W_t[idx2] + delta
    if(l == -Inf) {
      idx2 <- idx
    } else {
      idx2 <- tail(idx,-1)
    }
    bm$W_tm[idx2] <- bm$W_tm[idx2] + delta

    if(sum(bm$bb.local$t >= l & bm$bb.local$t <= r) > 0) {
      idx <- min(which(bm$bb.local$t >= l)):max(which(bm$bb.local$t <= r))
      if(r == Inf) {
        idx2 <- idx
      } else {
        idx2 <- head(idx,-1)
      }
      bm$bb.local$W_t[idx2]  <- bm$bb.local$W_t[idx2] + delta
      if(l == -Inf) {
        idx2 <- idx
      } else {
        idx2 <- tail(idx,-1)
      }
      bm$bb.local$W_tm[idx2] <- bm$bb.local$W_tm[idx2] + delta
    }

    bm$layers[bm$layers$t.l >= l & bm$layers$t.u <= r,4:7] <- bm$layers[bm$layers$t.l >= l & bm$layers$t.u <= r,4:7] + delta
    bm$user.layers[bm$user.layers$t.l >= l & bm$user.layers$t.u <= r,3:4] <- bm$user.layers[bm$user.layers$t.l >= l & bm$user.layers$t.u <= r,3:4] + delta
    bm$bb.local$layers[bm$bb.local$layers$t.l >= l & bm$bb.local$layers$t.u <= r,4:7] <- bm$bb.local$layers[bm$bb.local$layers$t.l >= l & bm$bb.local$layers$t.u <= r,4:7] + delta
    bm$bb.local$user.layers[bm$bb.local$user.layers$t.l >= l & bm$bb.local$user.layers$t.u <= r,3:4] <- bm$bb.local$user.layers[bm$bb.local$user.layers$t.l >= l & bm$bb.local$user.layers$t.u <= r,3:4] + delta

    add.labels_(bm, label, c(max(l,min(bm$t)),min(r,max(bm$t))))
  }
  bm$labels[["start"]] <- min(bm$t)
  bm$labels[["end"]] <- max(bm$t)
  bm$labels[["seg.start"]] <- c(bm$labels[["start"]], bm$t[which(bm$W_t != bm$W_tm)])
  bm$labels[["seg.end"]] <- c(bm$t[which(bm$W_t != bm$W_tm)], bm$labels[["end"]])

  invisible(bm)
}
