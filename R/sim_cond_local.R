#' Simulate Brownian motion conditional on localised layer
#'
#' Simulates Brownian motion
#'
#'
#' @param bm a Brownian motion object from which simulation should continue.
#'   Note the object is updated in place.
#'
#' @return the Brownian motion object which was passed in argument \code{bm} is
#'   updated in place and returned, enabling chaining of commands with
#'   dplyr (and other) style pipes.
#'
#' @export
sim.condlocal <- function(bm, s, t, q = NULL, q.grid = NULL) {
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
  localised <- bm$layers[bm$layers$type == "localised",]
  if(!(t %in% localised$t.u)) {
    stop("t is not the time of a realised maximum/minimum in a localised layer.  At present this package only supports conditional bridge simulation where the right-end conditioning time is a maximum or minimum.")
  }
  # Do we then need to check that s is the start of this same localised layer?
  # Actually this is taken care of by the following which allows no intermediate obs
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
  # Eliminate times we know
  q <- setdiff(q, bm$t)

  for(qq in q) {
    bm.res <- sim.condlocal_(bm, s_idx, qq, t_idx)
    s_idx <- s_idx+1
    t_idx <- t_idx+1
  }

  bm.res$layers <- bm.res$layers[order(bm$layers$t.l),]
  invisible(bm.res)
}

sim.condlocal_ <- function(bm, s_idx, q, t_idx) {

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
  minmax_idx <- match(tau, bm$layers$t.u)
  if(y == bm$layers$Ld[minmax_idx]) {
    minI <- +1
    bdry <- bm$layers$Uu[minmax_idx]
  } else if(y == bm$layers$Uu[minmax_idx]) {
    minI <- -1
    bdry <- bm$layers$Ld[minmax_idx]
  } else {
    stop("Error: min/max mismatch at tau")
  }

# sim.condlocal_ <- function(q, s, tau, x, y, minI, bdry) {


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

  # Update layer info
  # We split the layer in two, adding left and right layers either side of the
  # newly simulated time, then remove the old layer
  cur.layer <- which(bm$layers$t.l == s & bm$layers$t.u == tau)
  bm$layers <- add_row(bm$layers,
                       type = "intersection",
                       t.l = s,
                       t.u = q,
                       Ld = bm$layers[cur.layer,"Ld",TRUE],
                       Uu = bm$layers[cur.layer,"Uu",TRUE],
                       Lu = min(x, W_t),
                       Ud = max(x, W_t),
                       Lu.hard = TRUE, #ifelse(head(c(0, x.new), -1) > x.new, TRUE, FALSE)
                       Ud.hard = TRUE)
  bm$layers <- add_row(bm$layers,
                       type = "localised",
                       t.l = q,
                       t.u = tau,
                       Ld = bm$layers[cur.layer,"Ld",TRUE],
                       Uu = bm$layers[cur.layer,"Uu",TRUE],
                       Lu = min(y, W_t),
                       Ud = max(y, W_t),
                       Lu.hard = TRUE, #ifelse(head(c(0, x.new), -1) > x.new, TRUE, FALSE)
                       Ud.hard = TRUE)
  bm$layers <- bm$layers[-cur.layer,]

  bm
}
