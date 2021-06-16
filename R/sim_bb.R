#' Simulate a Brownian bridge
#'
#' Simulates Brownian motion conditional upon two (previously simulated) endpoints, outside of any layer.
#' Exactly one of the arguments \code{q} or \code{q.grid} should be supplied.
#'
#' @param bm a Brownian motion object from which simulation should continue.
#'   Note the object is updated in place.
#' @param s a scalar time upon which to condition the left hand end of the
#'   Brownian bridge.  This time must have already been simulated and be present
#'   in the Brownian motion object supplied in the \code{bm} argument.
#' @param t a scalar time upon which to condition the right hand end of the
#'   Brownian bridge.  This time must have already been simulated and be present
#'   in the Brownian motion object supplied in the \code{bm} argument.
#' @param q a vector of times at which to simulate the Brownian bridge, or NULL
#'   (the default) if the argument \code{q.grid} is to be used.  Note that all
#'   times must be strictly contained in the interval from \code{s} to \code{t}.
#' @param q.grid how many times to simulate the Brownian bridge at, which will
#'   be equally spaced between \code{s} and \code{t}, or NULL (the default) if
#'   the argument \code{q} is to be supplied.
#'
#' @return the Brownian motion object which was passed in argument \code{bm} is
#'   updated in place and returned, enabling chaining of commands with
#'   dplyr (and other) style pipes.
#'
#' @export
sim.bb <- function(bm, s, t, q = NULL, q.grid = NULL, label = names(q)) {
  # Arg types & combos
  if(!("BrownianMotion" %in% class(bm))) {
    stop("bm argument must be a BrownianMotion object.")
  }
  if(as.integer(is.null(q)) + as.integer(is.null(q.grid)) != 1) {
    stop("Exactly one of q or q.grid must be specified")
  }
  if(!is.null(q.grid) && !is.intscalar(q.grid)) {
    stop("q.grid must be either NULL or an integer number of grid points.")
  }
  if(!is.null(q) && !is.realvector(q)) {
    stop("q must be either NULL or a vector of times.")
  }
  if(!is.realscalar(t)) {
    stop("t must be a scalar.")
  }
  if(!is.realscalar(s)) {
    stop("s must be a scalar.")
  }

  # Checks
  if(!(s %in% bm$t)) {
    stop("s not found in bm path.")
  }
  if(!(t %in% bm$t)) {
    stop("t not found in bm path.")
  }
  s_idx <- match(s, bm$t)
  t_idx <- match(t, bm$t)
  if(t_idx != s_idx+1) {
    stop("Cannot have other brownian motion observations between times s and t")
  }

  # Make grid if it wasn't supplied
  if(is.null(q)) {
    q <- head(tail(seq(s, t, length.out = q.grid+2), -1), -1)
  }
  if(any(q < s) || any(q > t)) {
    stop("All q must lie between s and t")
  }
  # Label checks now we have q
  assert.bmlabel(label, q)
  if(!is.null(label) && length(label) == 1) {
    label <- rep(label, length(q))
  }
  if(is.unsorted(q)) {
    o <- order(q)
    label <- label[o]
    q <- q[o]
  }
  # Form sim vector
  q <- c(s, q, t)

  # Final check that the q values are not inside of a layer
  if(any(colSums(matrix(apply(bm$layers, 1, function(y) { q > y[2] & q < y[3] }), ncol = 2)) > 0)) {
    stop("A standard Brownian Bridge cannot be simulated inside of a layer.")
  }

  invisible(sim.bb_(bm, s_idx, q, t_idx, label))
}

sim.bb_ <- function(bm, s_idx, q, t_idx, label) {
  n <- length(q)

  W_t <- rep(0, n)
  W_t[1] <- bm$W_t[s_idx]
  W_t[n] <- bm$W_t[t_idx]

  for(i in 2:(n-1)) {
    W_t[i] <- rnorm(1,
                    mean = W_t[i-1] + (q[i]-q[i-1])/(q[n]-q[i-1])*(W_t[n]-W_t[i-1]),
                    sd = sqrt((q[n]-q[i])*(q[i]-q[i-1])/(q[n]-q[i-1])))
  }
  bm$t <- c(bm$t[1:s_idx],
            head(tail(q, -1), -1),
            bm$t[t_idx:length(bm$t)])
  bm$W_t <- c(bm$W_t[1:s_idx],
              head(tail(W_t, -1), -1),
              bm$W_t[t_idx:length(bm$W_t)])

  add.labels_(bm, "user", q[2:(n-1)])
  add.labels_(bm, label, q[2:(n-1)])

  bm
}
