#' Simulate Brownian motion conditional on bessel layer
#'
#' This currently converts the Bessel layer to an intersection layer and then
#' uses the conditional intersection layer sampling.
#'
#' @export
sim.condbessel <- function(bm, s, t, q = NULL, q.grid = NULL) {
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
  bessel <- bm$layers[bm$layers$type == "bessel",]
  if(!(t %in% bessel$t.u)) {
    stop("t is not the time of a realised maximum/minimum in a bessel layer.")
  }
  # Do we then need to check that s is the start of this same bessel layer?
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
    bm.res <- sim.condbessel_(bm, s_idx, qq, t_idx)
    s_idx <- s_idx+1
    t_idx <- t_idx+1
  }

  bm.res$layers <- bm.res$layers[order(bm.res$layers$t.l),]
  invisible(bm.res)
}

# Old version converts to intersection and then simulates
sim.condbesselOLD_ <- function(bm, s_idx, q, t_idx) {
  s <- bm$t[s_idx]
  t <- bm$t[t_idx]

  # Update layer info
  # We split the layer in two, adding left and right layers either side of the
  # newly simulated time, then remove the old layer
  cur.layer <- which(bm$layers$t.l == s & bm$layers$t.u == t)

  if(bm$layers$type[cur.layer] == "bessel") {
    bm <- intersection.layers.IL_(bm, s, t, 1)
    bm <- sim.condintersection_(bm, s_idx, q, t_idx)
  } else if(bm$layers$type[cur.layer] == "intersection") {
    bm <- sim.condintersection_(bm, s_idx, q, t_idx)
  } else {
    stop("Error: must have bessel or intersection layer at this point.")
  }

  bm
}

sim.condbessel_ <- function(bm, s_idx, q, t_idx) {

  s <- bm$t[s_idx]
  t <- bm$t[t_idx]
  x <- bm$W_t[s_idx]
  y <- bm$W_t[t_idx]
  cur.layer <- which(bm$layers$t.l == s & bm$layers$t.u == t)
  Ll <- bm$layers$Ld[cur.layer]
  Lu <- bm$layers$Lu[cur.layer]
  Ul <- bm$layers$Ud[cur.layer]
  Uu <- bm$layers$Uu[cur.layer]

  res <- sim.condbessel2_(q, s, t, x, y, Ll, Lu, Ul, Uu)

  bm$t <- c(bm$t[1:s_idx],
            q,
            bm$t[t_idx:length(bm$t)])
  bm$W_t <- c(bm$W_t[1:s_idx],
              res$w,
              bm$W_t[t_idx:length(bm$W_t)])

  # Update layer info
  # We split the layer in two, adding left and right layers either side of the
  # newly simulated time, then remove the old layer
  bm$layers <- add_row(bm$layers,
                       type = "bessel",
                       t.l = s,
                       t.u = q,
                       Ld = res$layer.sq[5],
                       Uu = res$layer.sq[8],
                       Lu = res$layer.sq[6],
                       Ud = res$layer.sq[7],
                       Lu.hard = ifelse(res$layer.sq[6] == min(x,res$w), TRUE, FALSE),
                       Ud.hard = ifelse(res$layer.sq[7] == max(x,res$w), TRUE, FALSE))
  bm$layers <- add_row(bm$layers,
                       type = "bessel",
                       t.l = q,
                       t.u = t,
                       Ld = res$layer.qt[5],
                       Uu = res$layer.qt[8],
                       Lu = res$layer.qt[6],
                       Ud = res$layer.qt[7],
                       Lu.hard = ifelse(res$layer.qt[6] == min(res$w,y), TRUE, FALSE),
                       Ud.hard = ifelse(res$layer.qt[7] == max(res$w,y), TRUE, FALSE))
  bm$layers <- bm$layers[-cur.layer,]

  bm
}

sim.condbessel2_ <- function(q, s, t, x, y, Ll, Lu, Ul, Uu) {
  # Repeat until acceptance
  repeat{
    # Simulate min/max for Bessel process
    mnmx <- eaminmax_(s,t,x,y,Ll,Lu,Ul,Uu); m<-mnmx$m; tau<-mnmx$tau; minI<-mnmx$minI
    if(q < tau){int.type <- 2}else{int.type <- 1} # Type 1 corresponds to tau occurring prior to q, Type 2 corresponds to tau occurring after q

    # Make draw at q conditional on m
    dr <- eabesmid_(q,s,tau,t,x,m,y,minI)$w

    # Realisation
    tp.mat <- matsort_(cbind(c(s,t,tau,q),c(x,y,m,dr)),1)

    ### Determine acceptance (and if so which layer) or rejection of the proposal
    type.vec <- numeric(3) # Storage of the numeric types
    u.ints <- runif(3,0,1) # Uniform RV's for intervals
    if(minI==1){bdry.l <- Ul; bdry.u <- Uu}else{bdry.l <- Lu; bdry.u <- Ll} # Determination of boundaries for interval
    m.counter <- rep(3,3) # Initialisation of the counter

    ##### Resolve interval 1

    repeat{
      bounds.l.one <- eadelC_(m.counter[1],tp.mat[1,1],tp.mat[2,1],tp.mat[1,2],tp.mat[2,2],m,bdry.l)  # Determine crossing probability of bdry.l for the interval one
      bounds.u.one <- eadelC_(m.counter[1],tp.mat[1,1],tp.mat[2,1],tp.mat[1,2],tp.mat[2,2],m,bdry.u)  # Determine crossing probability of bdry.u for the interval one

      if(u.ints[1] <= bounds.l.one[1]){type.vec[1] <- 1; break} # Bounds for interval one have been sufficiently resolved
      if(u.ints[1] >  bounds.u.one[2]){type.vec[1] <- 3; break} # Bounds for interval one have been sufficiently resolved
      if(u.ints[1] > bounds.l.one[2] & u.ints[1] <= bounds.u.one[1]){type.vec[1] <- 2; break} # Bounds for interval one have been sufficiently resolved
      m.counter[1] <- m.counter[1] + 2 # Otherwise index counter
    }

    ##### Resolve interval 2

    repeat{
      bounds.l.two <- eadelC_(m.counter[2],tp.mat[2,1],tp.mat[3,1],tp.mat[2,2],tp.mat[3,2],m,bdry.l)  # Determine crossing probability of bdry.l for the interval two
      bounds.u.two <- eadelC_(m.counter[2],tp.mat[2,1],tp.mat[3,1],tp.mat[2,2],tp.mat[3,2],m,bdry.u)  # Determine crossing probability of bdry.u for the interval two

      if(u.ints[2] <= bounds.l.two[1]){type.vec[2] <- 1; break} # Bounds for interval two have been sufficiently resolved
      if(u.ints[2] >  bounds.u.two[2]){type.vec[2] <- 3; break} # Bounds for interval two have been sufficiently resolved
      if(u.ints[2] > bounds.l.two[2] & u.ints[2] <= bounds.u.two[1]){type.vec[2] <- 2; break} # Bounds for interval two have been sufficiently resolved
      m.counter[2] <- m.counter[2] + 2 # Otherwise index counter
    }

    ##### Resolve interval 3

    repeat{
      bounds.l.three <- eadelC_(m.counter[3],tp.mat[3,1],tp.mat[4,1],tp.mat[3,2],tp.mat[4,2],m,bdry.l)  # Determine crossing probability of bdry.l for the interval three
      bounds.u.three <- eadelC_(m.counter[3],tp.mat[3,1],tp.mat[4,1],tp.mat[3,2],tp.mat[4,2],m,bdry.u)  # Determine crossing probability of bdry.u for the interval three

      if(u.ints[3] <= bounds.l.three[1]){type.vec[3] <- 1; break} # Bounds for interval three have been sufficiently resolved
      if(u.ints[3] >  bounds.u.three[2]){type.vec[3] <- 3; break} # Bounds for interval three have been sufficiently resolved
      if(u.ints[3] > bounds.l.three[2] & u.ints[3] <= bounds.u.three[1]){type.vec[3] <- 2; break} # Bounds for interval three have been sufficiently resolved
      m.counter[3] <- m.counter[3] + 2 # Otherwise index counter
    }

    #### Determine acceptance / rejection
    if(max(type.vec)==1){accI <- 1}else{if(max(type.vec)==2){if(rbinom(1,1,0.5)==1){accI <- 1}else{accI <- 0}}else{accI <- 0}}

    #### If accepted then break
    if(accI == 1){break}
  }

  ### Based on the appearance of the minimum / maximum determine laters for sq and qt

  if(int.type == 1){
    int.sq.type <- max(type.vec[1:2]); int.qt.type <- type.vec[3]
  }else{
    int.sq.type <- type.vec[1];int.qt.type <- max(type.vec[2:3])
  }

  #### Determine Layers
  if(minI == 1){ # Case where we have a minimum
    if(int.sq.type == 1){ # If interval sq is type 1
      layer.sq <- c(s,q,x,dr,Ll,min(Lu,x,dr),max(x,dr),bdry.l)
    }else{
      layer.sq <- c(s,q,x,dr,Ll,min(Lu,x,dr),max(x,dr,bdry.l),bdry.u)
    }
    if(int.qt.type == 1){ # If interval sq is type 1
      layer.qt <- c(q,t,dr,y,Ll,min(Lu,dr,y),max(dr,y),bdry.l)
    }else{
      layer.qt <- c(q,t,dr,y,Ll,min(Lu,dr,y),max(dr,y,bdry.l),bdry.u)
    }
  }else{
    if(int.sq.type == 1){ # If interval sq is type 1
      layer.sq <- c(s,q,x,dr,bdry.l,min(x,dr),max(Ul,x,dr),Uu)
    }else{
      layer.sq <- c(s,q,x,dr,bdry.u,min(x,dr,bdry.l),max(Ul,x,dr),Uu)
    }
    if(int.qt.type == 1){ # If interval sq is type 1
      layer.qt <- c(q,t,dr,y,bdry.l,min(dr,y),max(Ul,dr,y),Uu)
    }else{
      layer.qt <- c(q,t,dr,y,bdry.u,min(dr,y,bdry.l),max(Ul,dr,y),Uu)
    }
  }

  ### Output
  list(w=dr, layer.sq=layer.sq, layer.sq.type="Bessel", layer.qt=layer.qt, layer.qt.type="Bessel", u.ints = u.ints, l.bds = rbind(bounds.l.one,bounds.l.two,bounds.l.three), u.bds = rbind(bounds.u.one,bounds.u.two,bounds.u.three),mnmx=mnmx,tp.mat=tp.mat)
}
