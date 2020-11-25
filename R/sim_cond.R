#' Simulate Brownian motion conditional on
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
sim.cond <- function(bm, s, t, q = NULL, q.grid = NULL) {
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
  if(!(t %in% bm$bounds$t.u)) {
    stop("t is not the time of a realised maximum/minimum.  At present this package only supports conditional bridge simulation where the right-end conditioning time is a maximum or minimum.")
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
  if(is.unsorted(q)) {
    q <- sort(q)
  }
  for(qq in q) {
    bm.res <- sim.cond_(bm, s_idx, qq, t_idx)
    s_idx <- s_idx+1
    t_idx <- t_idx+1
  }

  invisible(bm.res)
}

sim.cond_ <- function(bm, s_idx, q, t_idx) {

# TODO: tidy this up, vectorise wrt q and remove dependency on scale pkg C++ code!

# Internal implementation of above callable without error checking
# q = new time
# s = left end time
# tau = attained min/max time
# x = W_s
# y = W_tau
# minI = +1 if tau is the time of attaining min, -1 if the time of attaining max
# bdry = value of the other boundary that is not attained at tau

  s <- bm$t[s_idx]
  tau <- bm$t[t_idx]
  x <- bm$W_t[s_idx]
  y <- bm$W_t[t_idx]
  minmax_idx <- match(tau, bm$bounds$t.u)
  if(y == bm$bounds$L[minmax_idx]) {
    minI <- +1
    bdry <- bm$bounds$U[minmax_idx]
  } else if(y == bm$bounds$U[minmax_idx]) {
    minI <- -1
    bdry <- bm$bounds$L[minmax_idx]
  } else {
    stop("Error: min/max mismatch at tau")
  }

# sim.cond_ <- function(q, s, tau, x, y, minI, bdry) {


  ## Repeat until acceptance
  repeat{
    ### Simulation of Bessel proposal
    c.scaling <- c(tau-s,minI*(x-y)*(tau-q)/(tau-s)^1.5,sqrt((tau-q)*(q-s))/(tau-s))
    bb.sims   <- rnorm(3,0,sd=c.scaling[3])
    w         <- y + minI*sqrt(c.scaling[1]*((c.scaling[2]+bb.sims[1])^2+bb.sims[2]^2+bb.sims[3]^2))

    ### Determine acceptance / rejection of proposal
    u       <- runif(1,0,1) # Uniform RV to make accept / reject decision
    counter <- ceiling(sqrt((tau-s)+(bdry-y)^2)/(2*abs(bdry-y))); if((counter%%2 == 0)==1){counter <- counter + 1} # Determining the minimum threshold for computing the boundary
    repeat{
      bounds <- eadelC_(counter,s,q,x,w,y,bdry)*eadelC_(counter,q,tau,w,y,y,bdry) # Determine crossing probability
      if(u <= bounds[1]) { # accept?
        accI <- 1
        break
      }
      if(u > bounds[2]) { # reject?
        accI <- 0
        break
      }
      counter <- counter + 2
    }

    ### If accept then break loop
    if(accI == 1) {
      break
    }
  }

  W_t <- w

  bm$t <- c(bm$t[1:s_idx],
            q,
            bm$t[t_idx:length(bm$t)])
  bm$W_t <- c(bm$W_t[1:s_idx],
              W_t,
              bm$W_t[t_idx:length(bm$W_t)])

  bm
}
