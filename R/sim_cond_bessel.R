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
  # Repeat until acceptance for the simulation of the mid-point
  repeat{
    # Simulate min/max for Bessel process
    mnmx <- eaminmax_(s,t,x,y,Ll,Lu,Ul,Uu); m<-mnmx$m; tau<-mnmx$tau; minI<-mnmx$minI
    # Determine whether auxiliary minimum / maximum and corresponding boundaries
    if(minI==1){B1<-Ul;B2<-Uu}else{B1<-Lu;B2<-Ll}
    # Make draw at q conditional on m
    dr <- eabesmid_(q,s,tau,t,x,m,y,minI)$w
    # Determine accept reject
    ar.prob <- eabesex_(c(s,q),c(q,t),c(x,dr),c(dr,y),m,B1,B2,minI)
    if(ar.prob$accI==1){break}
  }
  # Modification of Ul / Lu fiven new intermediate point
  nLuu.left <- min(x,dr)
  nLuu.right <- min(dr,y)
  nUll.left <- max(x,dr)
  nUll.right <- max(dr,y)
  nLu.left <- min(x,dr,Lu)
  nLu.right <- min(dr,y,Lu)
  nUl.left <- max(x,dr,Ul)
  nUl.right <- max(dr,y,Ul)

  # Determine left and right layers
  u.bisect <- runif(1,0,1) # Uniform for embedded rejection sampler
  m.counter <- 3 # Counter to detemine resolution of retrospective inversion sampler

  # Determine intervals: interval type 1
  repeat{
    p.int1 <- eabetaC_(m.counter,s,q,x,dr,Ll,nLu.left,nUl.left,Uu)*eabetaC_(m.counter,q,t,dr,y,Ll,nLu.right,nUl.right,Uu) # Case 1 probability
    p.int2 <- eabetaC_(m.counter,s,q,x,dr,nLu.left,nLuu.left,nUll.left,nUl.left)*eabetaC_(m.counter,q,t,dr,y,Ll,nLu.right,nUl.right,Uu) # Case 2 probability
    p.int3 <- eabetaC_(m.counter,s,q,x,dr,Ll,nLu.left,nUl.left,Uu)*eabetaC_(m.counter,q,t,dr,y,nLu.right,nLuu.right,nUll.right,nUl.right) # Case 3 probability
    denom <- p.int1+p.int2+p.int3 # Dominating probabilities
    if(u.bisect <= p.int1[1]/denom[1] | u.bisect > p.int1[2]/denom[1]){break} # If resolved sufficiently then break
    m.counter <- m.counter + 2 # Else index counter
  }

  # Determine intervals: interval type 2 / 3
  repeat{
    p.int1 <- eabetaC_(m.counter,s,q,x,dr,Ll,nLu.left,nUl.left,Uu)*eabetaC_(m.counter,q,t,dr,y,Ll,nLu.right,nUl.right,Uu) # Case 1 probability
    p.int2 <- eabetaC_(m.counter,s,q,x,dr,nLu.left,nLuu.left,nUll.left,nUl.left)*eabetaC_(m.counter,q,t,dr,y,Ll,nLu.right,nUl.right,Uu) # Case 2 probability
    p.int3 <- eabetaC_(m.counter,s,q,x,dr,Ll,nLu.left,nUl.left,Uu)*eabetaC_(m.counter,q,t,dr,y,nLu.right,nLuu.right,nUll.right,nUl.right) # Case 3 probability
    denom <- p.int1+p.int2+p.int3 # Dominating probabilities
    if(u.bisect <= (p.int1[1]+p.int2[1])/denom[1] | u.bisect > (p.int1[2]+p.int2[2])/denom[1]){break} # If resolved sufficiently then break
    m.counter <- m.counter + 2 # Else index counter
  }

  # Cases and outputting layers
  ## Case 1
  if(u.bisect <= p.int1[1]/denom[1]){
    layer.sq <- c(s,q,x,dr,Ll,nLu.left,nUl.left,Uu)
    layer.qt <- c(q,t,dr,y,Ll,nLu.right,nUl.right,Uu)
  }
  ## Case 2
  if(u.bisect > p.int1[1]/denom[1] & u.bisect <= (p.int1[2]+p.int2[2])/denom[1]){
    layer.sq <- c(s,q,x,dr,nLu.left,nLuu.left,nUll.left,nUl.left)
    layer.qt <- c(q,t,dr,y,Ll,nLu.right,nUl.right,Uu)
  }
  ## Case 3
  if(u.bisect > (p.int1[2]+p.int2[2])/denom[1]){
    layer.sq <- c(s,q,x,dr,Ll,nLu.left,nUl.left,Uu)
    layer.qt <- c(q,t,dr,y,nLu.right,nLuu.right,nUll.right,nUl.right)
  }

  ### Output
  list(w=dr, layer.sq=layer.sq, layer.sq.type="Bessel", layer.qt=layer.qt, layer.qt.type="Bessel", u.bisect = u.bisect)
}
