#' Manually specify the value of a point in the skeleton
#'
#' Only valid outside layers. Cannot be existing time (should delete then add)
#'
#' @param t time
#' @param W_t value
#'
#' @export
add.points <- function(bm, t, W_t, label = names(t)) {
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
