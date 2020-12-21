#' Simulate Brownian motion unconditionally
#'
#' Simulates Brownian motion unconditionally at times beyond the current
#' endpoint stored in the supplied object.  Providing times in \code{t} which
#' are not beyond the latest simulated endpoint will result in an error.
#'
#' @param bm a Brownian motion object from which simulation should continue.
#'   Note the object is updated in place.
#' @param t either a vector of times to simulate at, or a single endpoint time
#'   if t.grid is supplied.
#' @param t.grid by default this is NULL so that only the times supplied in
#'   \code{t} are simulated.  If a scalar integer is supplied then an evenly
#'   spaced grid of this many times will be simulated up to \code{t}.
#'
#' @return the Brownian motion object which was passed in argument \code{bm} is
#'   updated in place and returned, enabling chaining of commands with
#'   dplyr (and other) style pipes.
#'
#' @export
sim.uncond <- function(bm, t, t.grid = NULL) {
  # Arg types
  if(!("BrownianMotion" %in% class(bm))) {
    stop("bm argument must be a BrownianMotion object.")
  }
  if(!is.realscalar(t) && !is.realvector(t)) {
    stop("t must be either a scalar or vector of times.")
  }
  if(!is.null(t.grid) && !is.intscalar(t.grid)) {
    stop("t.grid must be either NULL or an integer number of grid points.")
  }

  # Arg combos
  if(is.realscalar(t) && !is.null(t.grid)) {
    t <- seq(tail(bm$t, 1), t, length.out = t.grid+1)[-1]
  }
  if(!is.realscalar(t) && !is.null(t.grid)) {
    stop("t.grid argument can only be used when t is a scalar.")
  }

  # Sanity
  if(any(tail(bm$t, 1) >= t)) {
    stop("Brownian motion has already passed some of the target times.  Please use a Brownian Bridge to infill.")
  }

  # Transformation
  if(is.unsorted(t)) {
    t <- sort(t)
  }

  # Eliminate times we know
  t <- setdiff(t, bm$t)

  invisible(sim.uncond_(bm, t))
}

# Internal implementation of above callable without error checking
sim.uncond_ <- function(bm, t) {
  t_delta <- c(t[1] - tail(bm$t, 1), diff(t))
  bm$t <- c(bm$t, t)

  W_delta <- rnorm(length(t_delta), sd = sqrt(t_delta))
  W_delta[1] <- W_delta[1] + tail(bm$W_t, 1)
  bm$W_t <- c(bm$W_t, cumsum(W_delta))

  bm
}
