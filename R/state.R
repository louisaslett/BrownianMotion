#' Return process state
#'
#' Returns the state of the Brownian motion at a given time.
#'
#' @param bm a Brownian motion object from which simulation should continue.
#'   Note the object is updated in place.
#' @param t a scalar time at which the state is required.  This time must have
#'   already been simulated and be present in the Brownian motion object
#'   supplied in the \code{bm} argument.  If omitted, the most recent state is
#'   returned.
#'
#' @return a vector containing the state at the specified time.
#'
#' @export
state <- function(bm, t = NULL) {
  if(!("BrownianMotion" %in% class(bm))) {
    stop("bm argument must be a BrownianMotion object.")
  }
  if(!is.null(t) && !is.realvector(t)) {
    stop("t must be either NULL or a scalar time.")
  }

  if(is.null(t)) {
    t <- max(bm$t)
  }

  invisible(state_(bm, t))
}

state_ <- function(bm, t) {
  t_idx <- match(t, bm$t)
  bm$W_t[t_idx]
}
