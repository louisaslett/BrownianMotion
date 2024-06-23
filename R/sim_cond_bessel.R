# Old version converts to intersection and then simulates
# sim.condbesselOLD_ <- function(bm, s_idx, q, t_idx) {
#   s <- bm$t[s_idx]
#   t <- bm$t[t_idx]
#
#   # Update layer info
#   # We split the layer in two, adding left and right layers either side of the
#   # newly simulated time, then remove the old layer
#   cur.layer <- which(bm$layers$t.l == s & bm$layers$t.u == t)
#
#   if(bm$layers$type[cur.layer] == "bessel") {
#     bm <- intersection.layers.IL_(bm, s, t, 1)
#     bm <- sim.condintersection_(bm, s_idx, q, t_idx)
#   } else if(bm$layers$type[cur.layer] == "intersection") {
#     bm <- sim.condintersection_(bm, s_idx, q, t_idx)
#   } else {
#     stop("Error: must have bessel or intersection layer at this point.")
#   }
#
#   bm
# }

sim.condbessel_ <- function(bm, s_idx, q, t_idx, label, dr = NA) {

  s <- bm$t[s_idx]
  t <- bm$t[t_idx]
  x <- bm$W_t[s_idx]
  y <- bm$W_tm[t_idx]
  cur.layer <- which(bm$layers$t.l == s & bm$layers$t.u == t)
  Ll <- bm$layers$Ld[cur.layer]
  Lu <- bm$layers$Lu[cur.layer]
  Ul <- bm$layers$Ud[cur.layer]
  Uu <- bm$layers$Uu[cur.layer]

  if(is.na(dr)) {
    res <- sim.condbessel2_(q, s, t, x, y, Ll, Lu, Ul, Uu)
  } else {
    res <- sim.condbessel3_(dr, q, s, t, x, y, Ll, Lu, Ul, Uu)
  }

  bm$t <- c(bm$t[1:s_idx],
            q,
            bm$t[t_idx:length(bm$t)])
  bm$W_t <- c(bm$W_t[1:s_idx],
              res$w,
              bm$W_t[t_idx:length(bm$W_t)])
  bm$W_tm <- c(bm$W_tm[1:s_idx],
               res$w,
               bm$W_tm[t_idx:length(bm$W_tm)])

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

  add.labels_(bm, "user", q)
  add.labels_(bm, label, q)

  bm
}

sim.condbessel2_ <- function(q, s, t, x, y, Ll, Lu, Ul, Uu) {
  # Repeat until acceptance for the simulation of the mid-point
  cnt <- 0
  cnt2 <- 0
  cnt3 <- c()
  repeat{
    # Simulate min/max for Bessel process
    mnmx <- eaminmax_(s,t,x,y,Ll,Lu,Ul,Uu); m<-mnmx$m; tau<-mnmx$tau; minI<-mnmx$minI
    # Determine whether auxiliary minimum / maximum and corresponding boundaries
    if(minI==1){B1<-Ul;B2<-Uu}else{B1<-Lu;B2<-Ll}
    # Make draw at q conditional on m
    dr <- eabesmid_(q,s,tau,t,x,m,y,minI)$w
    # Determine accept reject
    cnt <- cnt+1
    if((m<dr & dr>B2) | (m>dr & dr<B2))
      cnt2 <- cnt2+1
    if(cnt > 1000000)
      browser()
    if(!((m<dr & dr>B2) | (m>dr & dr<B2))) {
      ar.prob <- eabesex_(c(s,q),c(q,t),c(x,dr),c(dr,y),m,B1,B2,minI)
      cnt3 <- c(cnt3, prod(ar.prob$em[,8]))
      if(ar.prob$accI==1){break}
    }
  }

  # HACK
  # Modification of Ul / Lu given new intermediate point
  nLuu.left <- min(x,dr)
  nLuu.right <- min(dr,y)
  nUll.left <- max(x,dr)
  nUll.right <- max(dr,y)
  nLu.left <- min(x,dr,Lu)
  nLu.right <- min(dr,y,Lu)
  nUl.left <- max(x,dr,Ul)
  nUl.right <- max(dr,y,Ul)

  if(tau < q){
    layer.sq <- c(s,q,x,dr,Ll,nLu.left,nUl.left,Uu)
    layer.qt <- c(q,t,dr,y,Ll,nLuu.right,nUll.right,Uu)
  } else {
    layer.sq <- c(s,q,x,dr,Ll,nLuu.left,nUll.left,Uu)
    layer.qt <- c(q,t,dr,y,Ll,nLu.right,nUl.right,Uu)
  }

  list(w=dr, layer.sq=layer.sq, layer.sq.type="Bessel", layer.qt=layer.qt, layer.qt.type="Bessel")
  # sim.condbessel3_(dr, q, s, t, x, y, Ll, Lu, Ul, Uu)
}

sim.condbessel3_ <- function(dr, q, s, t, x, y, Ll, Lu, Ul, Uu) {
  # Modification of Ul / Lu given new intermediate point
  nLuu.left <- min(x,dr)
  nLuu.right <- min(dr,y)
  nUll.left <- max(x,dr)
  nUll.right <- max(dr,y)
  nLu.left <- min(x,dr,Lu)
  nLu.right <- min(dr,y,Lu)
  nUl.left <- max(x,dr,Ul)
  nUl.right <- max(dr,y,Ul)

  # Check whether the layers co-incide with the end points

  left.coinc <- (x==nLu.left  | x==nUl.left  | dr==nLu.left  | dr==nUl.left)
  right.coinc <- (y==nLu.right | y==nUl.right | dr==nLu.right | dr==nUl.right)

  if(left.coinc && right.coinc) { # Case where either end- or mid-points co-incide with the layer
    # Case 1
    layer.sq <- c(s,q,x,dr,Ll,nLu.left,nUl.left,Uu)
    layer.qt <- c(q,t,dr,y,Ll,nLu.right,nUl.right,Uu)
  } else if(left.coinc && !right.coinc) { # Case where left end- or mid-points co-incide with the layer
    # Bisection variables
    u.bisect <- runif(1,0,1) # Uniform for embedded rejection sampler
    m.counter <- 3 # Counter to determine resolution of retrospective inversion sample

    # Determine intervals
    repeat{
      p.int1 <- eabetaC_(m.counter,q,t,dr,y,Ll,nLu.right,nUl.right,Uu)[1:2] # Case 1 probability
      p.int3 <- eabetaC_(m.counter,q,t,dr,y,nLu.right,nLuu.right,nUll.right,nUl.right)[1:2] # Case 3 probability
      denom <- p.int1+p.int3 # Dominating probabilities
      if(denom[1] > 0){ if(u.bisect <= p.int1[1]/denom[2] | u.bisect > p.int1[2]/denom[1]){break} } # If resolved sufficiently then break
      m.counter <- m.counter + 2 # Else index counter
    }
    # Cases and outputting layers
    ## Case 1
    if(u.bisect <= p.int1[1]/denom[2]){
      layer.sq <- c(s,q,x,dr,Ll,nLu.left,nUl.left,Uu)
      layer.qt <- c(q,t,dr,y,Ll,nLu.right,nUl.right,Uu)
    }
    ## Case 3
    if(u.bisect > p.int1[1]/denom[2]){
      layer.sq <- c(s,q,x,dr,Ll,nLu.left,nUl.left,Uu)
      layer.qt <- c(q,t,dr,y,nLu.right,nLuu.right,nUll.right,nUl.right)
    }
  } else if(!left.coinc && right.coinc) { # Case where right end- or mid-points co-incide with the layer
    # Bisection variables
    u.bisect <- runif(1,0,1) # Uniform for embedded rejection sampler
    m.counter <- 3 # Counter to detemine resolution of retrospective inversion sample

    # Determine intervals
    repeat{
      p.int1 <- eabetaC_(m.counter,s,q,x,dr,Ll,nLu.left,nUl.left,Uu)[1:2] # Case 1 probability
      p.int2 <- eabetaC_(m.counter,s,q,x,dr,nLu.left,nLuu.left,nUll.left,nUl.left)[1:2] # Case 2 probability
      denom <- p.int1+p.int2 # Dominating probabilities
      if(denom[1] > 0){ if(u.bisect <= p.int1[1]/denom[2] | u.bisect > p.int1[2]/denom[1]){break} } # If resolved sufficiently then break
      m.counter <- m.counter + 2 # Else index counter
    }
    # Cases and outputting layers
    ## Case 1
    if(u.bisect <= p.int1[1]/denom[2]){
      layer.sq <- c(s,q,x,dr,Ll,nLu.left,nUl.left,Uu)
      layer.qt <- c(q,t,dr,y,Ll,nLu.right,nUl.right,Uu)
    }
    ## Case 2
    if(u.bisect > p.int1[1]/denom[2]){
      layer.sq <- c(s,q,x,dr,nLu.left,nLuu.left,nUll.left,nUl.left)
      layer.qt <- c(q,t,dr,y,Ll,nLu.right,nUl.right,Uu)
    }
  } else if(!left.coinc && !right.coinc) { # Case where end- or mid-points co-incide with the layer
    # Determine left and right layers
    u.bisect <- runif(1,0,1) # Uniform for embedded rejection sampler
    m.counter <- 3 # Counter to detemine resolution of retrospective inversion sampler

    # Determine intervals: interval type 1
    repeat{
      p.int1 <- eabetaC_(m.counter,s,q,x,dr,Ll,nLu.left,nUl.left,Uu)[1:2]*eabetaC_(m.counter,q,t,dr,y,Ll,nLu.right,nUl.right,Uu)[1:2] # Case 1 probability
      p.int2 <- eabetaC_(m.counter,s,q,x,dr,nLu.left,nLuu.left,nUll.left,nUl.left)[1:2]*eabetaC_(m.counter,q,t,dr,y,Ll,nLu.right,nUl.right,Uu)[1:2] # Case 2 probability
      p.int3 <- eabetaC_(m.counter,s,q,x,dr,Ll,nLu.left,nUl.left,Uu)[1:2]*eabetaC_(m.counter,q,t,dr,y,nLu.right,nLuu.right,nUll.right,nUl.right)[1:2] # Case 3 probability
      denom <- p.int1+p.int2+p.int3 # Dominating probabilities
      if(denom[1] > 0){ if(u.bisect <= p.int1[1]/denom[2] | u.bisect > p.int1[2]/denom[1]){break} } # If resolved sufficiently then break
      m.counter <- m.counter + 2 # Else index counter
    }

    # Determine intervals: interval type 2 / 3
    repeat{
      p.int1 <- eabetaC_(m.counter,s,q,x,dr,Ll,nLu.left,nUl.left,Uu)[1:2]*eabetaC_(m.counter,q,t,dr,y,Ll,nLu.right,nUl.right,Uu)[1:2] # Case 1 probability
      p.int2 <- eabetaC_(m.counter,s,q,x,dr,nLu.left,nLuu.left,nUll.left,nUl.left)[1:2]*eabetaC_(m.counter,q,t,dr,y,Ll,nLu.right,nUl.right,Uu)[1:2] # Case 2 probability
      p.int3 <- eabetaC_(m.counter,s,q,x,dr,Ll,nLu.left,nUl.left,Uu)[1:2]*eabetaC_(m.counter,q,t,dr,y,nLu.right,nLuu.right,nUll.right,nUl.right)[1:2] # Case 3 probability
      denom <- p.int1+p.int2+p.int3 # Dominating probabilities
      if(denom[1] > 0){ if(u.bisect <= (p.int1[1]+p.int2[1])/denom[2] | u.bisect > (p.int1[2]+p.int2[2])/denom[1]){break} } # If resolved sufficiently then break
      m.counter <- m.counter + 2 # Else index counter
    }

    # Cases and outputting layers
    ## Case 1
    if(u.bisect <= p.int1[1]/denom[2]){
      layer.sq <- c(s,q,x,dr,Ll,nLu.left,nUl.left,Uu)
      layer.qt <- c(q,t,dr,y,Ll,nLu.right,nUl.right,Uu)
    }
    ## Case 2
    if(u.bisect > p.int1[1]/denom[2] & u.bisect <= (p.int1[1]+p.int2[1])/denom[2]){
      layer.sq <- c(s,q,x,dr,nLu.left,nLuu.left,nUll.left,nUl.left)
      layer.qt <- c(q,t,dr,y,Ll,nLu.right,nUl.right,Uu)
    }
    ## Case 3
    if(u.bisect > (p.int1[1]+p.int2[1])/denom[2]){
      layer.sq <- c(s,q,x,dr,Ll,nLu.left,nUl.left,Uu)
      layer.qt <- c(q,t,dr,y,nLu.right,nLuu.right,nUll.right,nUl.right)
    }
 }
  ### Output
  list(w=dr, layer.sq=layer.sq, layer.sq.type="Bessel", layer.qt=layer.qt, layer.qt.type="Bessel")
}
