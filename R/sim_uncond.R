#' Simulate Brownian motion unconditionally
#'
#' Simulates Brownian motion unconditionally at times beyond the current
#' endpoint stored in the supplied object.  Providing times in \code{t} which
#' are not beyond the latest simulated endpoint will result in an error.
#'
#' @param bm a Brownian motion object from which simulation should continue.
#'   Note the object is updated in place.
#' @param t  a vector of times at which to simulate the path
#'
#' @return the Brownian motion object which was passed in argument \code{bm} is
#'   updated in place and returned, enabling chaining of commands with
#'   dplyr (and other) style pipes.
#'
#' @export
sim.uncond <- function(bm, t, label = names(t)) {
  # Arg types
  if(!("BrownianMotion" %in% class(bm))) {
    stop("bm argument must be a BrownianMotion object.")
  }
  assert.timevector(t)
  if(any(t < min(bm$t))) {
    stop(paste0("cannot simulate path at times before the path was initialised (at time ", min(bm$t), ")"))
  }
  if(any(t < max(bm$t))) {
    stop(paste0("cannot unconditionally simulate path at times before the end of the path (at time ", max(bm$t), ")"))
  }
  assert.bmlabel(label, t)

  # Transformation
  if(!is.null(label) && length(label) == 1) {
    label <- rep(label, length(t))
  }
  if(is.unsorted(t)) {
    o <- order(t)
    label <- label[o]
    t <- t[o]
  }

  # Eliminate times we know
  t <- setdiff(t, bm$t)

  invisible(sim.uncond_(bm, t, label))
}

# Internal implementation of above callable without error checking
sim.uncond_ <- function(bm, t, label) {
  t_delta <- c(t[1] - tail(bm$t, 1), diff(t))
  bm$t <- c(bm$t, unname(t))

  W_delta <- rnorm(length(t_delta), sd = sqrt(t_delta))
  W_delta[1] <- W_delta[1] + tail(bm$W_t, 1)
  bm$W_t <- c(bm$W_t, cumsum(W_delta))

  add.labels_(bm, "user", t)
  add.labels_(bm, label, t)
  bm$labels[["end"]] <- max(t)
  bm
}
