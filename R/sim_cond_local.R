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
    stop("t is not the end time of a localised layer.")
  }
  # NOTE: Do we then need to check that s is the start of this same localised layer? NO ...
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

  lyr.hard <- localised[localised$t.u==t, c("Lu.hard","Ud.hard")]
  for(qq in q) {
    if(any(!lyr.hard)) {
      sim.condlocalsoft_(bm, s_idx, qq, t_idx)
      lyr.hard <- tail(bm$layers)[,c("Lu.hard","Ud.hard")]
    } else {
      sim.condlocal_(bm, s_idx, qq, t_idx)
      # don't need to update as hard stays hard
    }
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
                       Ud.hard = TRUE) # Always hard ... "conditional simulation is the viagra of Brownian motion" - M Pollock, 3/2/2021
  bm$layers <- bm$layers[-cur.layer,]

  bm
}

sim.condlocalsoft_ <- function(bm, s_idx, q, t_idx) {
  s <- bm$t[s_idx]
  t <- bm$t[t_idx]
  x <- bm$W_t[s_idx]
  y <- bm$W_t[t_idx]

  lyr_idx <- match(t, bm$layers$t.u)
  Ll <- bm$layers$Ld[lyr_idx]
  Lu <- bm$layers$Lu[lyr_idx]
  Ul <- bm$layers$Ud[lyr_idx]
  Uu <- bm$layers$Uu[lyr_idx]

  if(y == bm$layers$Ld[lyr_idx]) {
    minI <- +1
  } else if(y == bm$layers$Uu[lyr_idx]) {
    minI <- -1
  } else {
    stop("Error: min/max mismatch at right end of layer")
  }

  # Minimum indicator = + 1 if minimum, -1 if maximum
  # Repeat until acceptance
  repeat{
    ### Simulation of Bessel proposal
    c.scaling     <- c(t-s,minI*(x-y)*(t-q)/(t-s)^1.5,sqrt((t-q)*(q-s))/(t-s))
    bb.sims       <- rnorm(3,0,sd=c.scaling[3])
    w             <- y + minI*sqrt(c.scaling[1]*((c.scaling[2]+bb.sims[1])^2+bb.sims[2]^2+bb.sims[3]^2))

    ### Determine acceptance (and if so which layer) or rejection of the proposal
    u.sq <- runif(1,0,1); u.qt <- runif(1,0,1) # Uniform RV's for the left and right hand intervals
    if(minI==1){bdry.l <- Ul; bdry.u <- Uu}else{bdry.l <- Lu; bdry.u <- Ll} # Determination of boundaries for interval
    counter.sq <- counter.qt <- 3 # Initialisation of the counter

    ##### Resolve interval [s,q]

    repeat{
      bounds.l.sq <- eadelC_(counter.sq,s,q,x,w,y,bdry.l)  # Determine crossing probability of bdry.l for the interval [s,q]
      bounds.u.sq <- eadelC_(counter.sq,s,q,x,w,y,bdry.u)  # Determine crossing probability of bdry.u for the interval [s,q]

      if(u.sq <= bounds.l.sq[1] | u.sq > bounds.u.sq[2]){break} # Bounds for interval [s,q] have been sufficiently resolved
      if(u.sq > bounds.l.sq[2] & u.sq <= bounds.u.sq[1]){break} # Bounds for interval [s,q] have been sufficiently resolved
      counter.sq <- counter.sq + 2 # Otherwise index counter
    }

    ##### Resolve interval [q,t]

    repeat{ # Resolve interval [q,t]
      bounds.l.qt <- eadelC_(counter.qt,q,t,w,y,y,bdry.l)  # Determine crossing probability of bdry.l for the interval [q,t]
      bounds.u.qt <- eadelC_(counter.qt,q,t,w,y,y,bdry.u)  # Determine crossing probability of bdry.u for the interval [q,t]

      if(u.qt <= bounds.l.qt[1] | u.qt > bounds.u.qt[2]){break} # Bounds for interval [q,t] have been sufficiently resolved
      if(u.qt > bounds.l.qt[2] & u.qt <= bounds.u.qt[1]){break} # Bounds for interval [q,t] have been sufficiently resolved
      counter.qt <- counter.qt + 2 # Otherwise index counter
    }

    ##### Determine acceptance / rejection and output layers of applicable

    ########## Case 1: Either u.sq and u.qt determine rejection by exceedence
    if(u.sq >= bounds.u.sq[2] | u.qt >= bounds.u.qt[2]){accI <- 0}

    ########## Case 2: Both u.sq and u.qt determine rejection by not meeting inner layer
    if(u.sq <= bounds.l.sq[1] & u.qt <= bounds.l.qt[1]){accI <- 0}

    ########## Case 3: u.sq and u.qt fall in interval
    if(u.sq > bounds.l.sq[2] & u.sq <= bounds.u.sq[1] & u.qt > bounds.l.qt[2] & u.qt <= bounds.u.qt[1]){
      accI <- 1 # Accept
      if(minI==1){ # Determine layers
        layer.sq <- c(s,q,x,w,Ll,min(x,w),max(x,w,bdry.l),bdry.u,0) # Determine [s,q] layer information
        layer.qt <- c(q,t,w,y,Ll,Lu,max(w,y,bdry.l),bdry.u,minI) # Determine [q,t] layer information
      }else{
        layer.sq <- c(s,q,x,w,bdry.u,min(x,w,bdry.l),max(x,w),Uu,0) # Determine [s,q] layer information
        layer.qt <- c(q,t,w,y,bdry.u,min(w,y,bdry.l),Ul,Uu,minI) # Determine [q,t] layer information
      }
    }

    ########## Case 4: u.sq falls in interval and u.qt does not reach or exceed interval
    if(u.sq > bounds.l.sq[2] & u.sq <= bounds.u.sq[1] & u.qt <= bounds.l.qt[1]){
      accI <- 1 # Accept
      if(minI==1){ # Determine layers
        layer.sq <- c(s,q,x,w,Ll,min(x,w),max(x,w,bdry.l),bdry.u,0) # Determine [s,q] layer information
        layer.qt <- c(q,t,w,y,Ll,Lu,max(w,y),max(w,y,bdry.l),minI) # Determine [q,t] layer information
      }else{
        layer.sq <- c(s,q,x,w,bdry.u,min(x,w,bdry.l),max(x,w),Uu,0) # Determine [s,q] layer information
        layer.qt <- c(q,t,w,y,min(w,y,bdry.l),min(w,y),Ul,Uu,minI) # Determine [q,t] layer information
      }
    }
    ########## Case 5: u.qt falls in interval and u.sq does not reach or exceed interval
    if(u.qt > bounds.l.qt[2] & u.qt <= bounds.u.qt[1] & u.sq <= bounds.l.sq [1]){
      accI <- 1 # Accept
      if(minI==1){ # Determine layers
        layer.sq <- c(s,q,x,w,Ll,min(x,w),max(x,w),max(x,w,bdry.l),0) # Determine [s,q] layer information
        layer.qt <- c(q,t,w,y,Ll,Lu,max(w,y,bdry.l),bdry.u,minI) # Determine [q,t] layer information
      }else{
        layer.sq <- c(s,q,x,w,min(x,w,bdry.l),min(x,w),max(x,w),Uu,0) # Determine [s,q] layer information
        layer.qt <- c(q,t,w,y,bdry.u,min(w,y,bdry.l),Ul,Uu,minI) # Determine [q,t] layer information
      }
    }

    ### If accepted then break
    if(accI == 1){break}
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
  cur.layer <- which(bm$layers$t.l == s & bm$layers$t.u == t)
  bm$layers <- add_row(bm$layers,
                       type = "intersection",
                       t.l = layer.sq[1],
                       t.u = layer.sq[2],
                       Ld = layer.sq[5],
                       Uu = layer.sq[8],
                       Lu = layer.sq[6],
                       Ud = layer.sq[7],
                       Lu.hard = TRUE,
                       Ud.hard = TRUE)
  bm$layers <- add_row(bm$layers,
                       type = "localised",
                       t.l = layer.qt[1],
                       t.u = layer.qt[2],
                       Ld = layer.qt[5],
                       Uu = layer.qt[8],
                       Lu = layer.qt[6],
                       Ud = layer.qt[7],
                       Lu.hard = ifelse(minI==1 | layer.qt[6]==W_t, TRUE, FALSE),
                       Ud.hard = ifelse(minI!=1 | layer.qt[7]==W_t, TRUE, FALSE))
  bm$layers <- bm$layers[-cur.layer,]

# minimum & Ud==w
# want: - Lu.hard = T, Ud.hard = T
# alg LA: Lu.hard = T, Ud.hard = T
# alg MP:
#
# minimum & Ud > w
# want: - - Lu.hard = T, Ud.hard = F
# alg LA: Lu.hard = T, Ud.hard = F
#   alg MP:
#
# maximum & Lu = w
# want: - - Lu.hard = T, Ud.hard = T
# alg LA: Lu.hard = T, Ud.hard = T
#   alg MP:
#
# maximum & Lu < w
# want: - - Lu.hard = F, Ud.hard = T
# alg LA: Lu.hard = F, Ud.hard = T
#   alg MP:
### Output
# list(w=w, layer.sq=layer.sq, layer.sq.type="Intersection", layer.qt=layer.qt, layer.qt.type=if(minI==1){if(layer.qt[7]==w){"Hard Localised"}else{"Soft Localised"}}else{if(layer.qt[6]==w){"Hard Localised"}else{"Soft Localised"}},u.sq=u.sq, sq.bds=c(bounds.l.sq,bounds.u.sq), u.qt = u.qt, qt.bds = c(bounds.l.qt,bounds.u.qt))}

  bm
}
