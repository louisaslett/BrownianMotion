#' Refine an intersection layer
#'
#'
#'
#' @export
refine <- function(bm, s, t, mult = bm$mult) {
  UseMethod("refine")
}

#' @export
refine.BrownianMotionNd <- function(bm, ...) {
  if(!("BrownianMotionNd" %in% class(bm))) {
    stop("bm argument must be a BrownianMotionNd object.")
  }

  for(d in 1:bm$dim) {
    refine(bm$Z.bm[[d]], ...)
  }

  invisible(bm)
}

#' @export
refine.BrownianMotion <- function(bm, s, t, mult = bm$mult) {
  # TODO: tidy up and add error handling, this is a quick prototype

  lyrs <- which(bm$layers$t.u > s & bm$layers$t.l<t)
  for(l in lyrs) {
    if(bm$layers$type[l] == "intersection") {
      refine.intersection_(bm, l, mult)
    } else if(bm$layers$type[l] == "localised") {
      refine.local_(bm, l, mult)
    } else if(bm$layers$type[l] == "bessel") {
      refine.bessel_(bm, l, mult)
    } else if(bm$layers$type[l] == "intersection-bb") {
      refine.bbintersection_(bm, l, mult)
    } else if(bm$layers$type[l] == "localised-bb") {
      refine.bblocal_(bm, l, mult)
    } else if(bm$layers$type[l] == "bessel-bb") {
      refine.bbbessel_(bm, l, mult)
    } else {
      warning("Cannot refine layer on interval [{bm$layers$t.l[l]},{bm$layers$t.u[l]}] of type {bm$layers$type[l]}.\n")
    }
  }

  invisible(bm)
}

refine_ <- function(bm, lyr.idx, mult) {
  if(bm$layers$type[lyr.idx] == "intersection") {
    refine.intersection_(bm, lyr.idx, mult)
  } else if(bm$layers$type[lyr.idx] == "localised") {
    refine.local_(bm, lyr.idx, mult)
  } else if(bm$layers$type[lyr.idx] == "bessel") {
    refine.bessel_(bm, lyr.idx, mult)
  } else if(bm$layers$type[lyr.idx] == "intersection-bb") {
    refine.bbintersection_(bm, lyr.idx, mult)
  } else if(bm$layers$type[lyr.idx] == "localised-bb") {
    refine.bblocal_(bm, lyr.idx, mult)
  } else if(bm$layers$type[lyr.idx] == "bessel-bb") {
    refine.bbbessel_(bm, lyr.idx, mult)
  } else {
    warning("Cannot refine layer on interval [{bm$layers$t.l[lyr.idx]},{bm$layers$t.u[lyr.idx]}] of type {bm$layers$type[lyr.idx]}.\n")
  }

  invisible(bm)
}

refine.intersection_ <- function(bm, lyr.idx, mult) {
  s <- bm$layers[lyr.idx,"t.l",drop=TRUE]
  t <- bm$layers[lyr.idx,"t.u",drop=TRUE]
  x <- bm$W_t[match(s, bm$t)]
  y <- bm$W_tm[match(t, bm$t)]
  Ll <- bm$layers[lyr.idx,"Ld",drop=TRUE]
  Lu <- bm$layers[lyr.idx,"Lu",drop=TRUE]
  Ul <- bm$layers[lyr.idx,"Ud",drop=TRUE]
  Uu <- bm$layers[lyr.idx,"Uu",drop=TRUE]

  if(max((Uu-Ul), (Lu-Ll)) > mult*sqrt(t-s)) {
    while(max((Uu-Ul), (Lu-Ll)) > mult*sqrt(t-s)) {
      Ls<-Ll+(Lu-Ll)/2
      Us<-Uu-(Uu-Ul)/2
      mat <- matrix(c(Ll,Ls,Us,Uu,Ls,Lu,Us,Uu,Ll,Ls,Ul,Us,Ls,Lu,Ul,Us), 4, 4, byrow=TRUE)
      bbind <- deind <- 0
      bbct <- 1
      m1 <- 5
      u1 <- runif(1, 0, 1)
      while(bbind == 0) {
        while(deind == 0) {
          debd <- eabetaC_(m1,s,t,x,y,Ll,Lu,Ul,Uu)[2:1]
          if(debd[2] <= 0) {
            m1 <- m1+2
          } else {
            deind<-1
          }
        }
        be <- matrix(0, bbct, 2)
        for(i in 1:bbct) {
          be[i,] <- eabetaC_(m1,s,t,x,y,mat[i,1],mat[i,2],mat[i,3],mat[i,4])[1:2]
        }
        bd <- c(sum(be[,1])/debd[1],sum(be[,2])/debd[2])
        if(u1 <= bd[1]) {
          bbind <- 1
        } else {
          if(u1 >= bd[2]) {
            bbct <- bbct+1
            deind <- 0
          } else {
            m1 <- m1+2
            deind <- 0
          }
        }
      }
      Ll <- mat[bbct,1]
      Lu <- mat[bbct,2]
      Ul <- mat[bbct,3]
      Uu <- mat[bbct,4]
    }
  }

  bm$layers[lyr.idx,"Ld"] <- Ll
  bm$layers[lyr.idx,"Lu"] <- Lu
  bm$layers[lyr.idx,"Ud"] <- Ul
  bm$layers[lyr.idx,"Uu"] <- Uu

  bm
}

refine.bbintersection_ <- function(bm, lyr.idx, mult) {
  bblyr.idx <- which(bm$bb.local$layers$t.l == bm$layers$t.l[lyr.idx] &
                       bm$bb.local$layers$t.u == bm$layers$t.u[lyr.idx])

  if(bm$bb.local$layers$type[bblyr.idx] != "intersection")
    stop("Auxilliary and main path layer mismatch detected.")

  refine.intersection_(bm$bb.local, bblyr.idx, mult)

  s_idx <- match(bm$layers$t.l[lyr.idx], bm$t)
  t_idx <- match(bm$layers$t.u[lyr.idx], bm$t)
  bb.s_idx <- match(bm$layers$t.l[lyr.idx], bm$bb.local$t)
  bb.t_idx <- match(bm$layers$t.u[lyr.idx], bm$bb.local$t)

  Bs <- bm$W_t[s_idx]
  Bt <- bm$W_tm[t_idx]
  Ws <- bm$bb.local$W_t[bb.s_idx]
  Wt <- bm$bb.local$W_t[bb.t_idx]

  # bm$layers[lyr.idx,"Ld"] <- Bs + bm$bb.local$layers[bblyr.idx,"Ld"] - Ws
  # bm$layers[lyr.idx,"Lu"] <- Bs + bm$bb.local$layers[bblyr.idx,"Ld"] - Ws + ((Bt-Bs)-(Wt-Ws))
  # bm$layers[lyr.idx,"Ud"] <- Bs + bm$bb.local$layers[bblyr.idx,"Uu"] - Ws
  # bm$layers[lyr.idx,"Uu"] <- Bs + bm$bb.local$layers[bblyr.idx,"Uu"] - Ws + ((Bt-Bs)-(Wt-Ws))
  bm$layers[lyr.idx,"Ld"] <- Bs - Ws
  bm$layers[lyr.idx,"Lu"] <- NA
  bm$layers[lyr.idx,"Ud"] <- NA
  bm$layers[lyr.idx,"Uu"] <- Bs - Ws + ((Bt-Bs)-(Wt-Ws))
  bm$layers[lyr.idx,"Lu.hard"] <- bm$bb.local$layers[bblyr.idx,"Lu.hard"]
  bm$layers[lyr.idx,"Ud.hard"] <- bm$bb.local$layers[bblyr.idx,"Ud.hard"]

  bm
}

refine.local_ <- function(bm, lyr.idx, mult) {
  s <- bm$layers[lyr.idx,"t.l",drop=TRUE]
  t <- bm$layers[lyr.idx,"t.u",drop=TRUE]
  x <- bm$W_t[match(s, bm$t)]
  y <- bm$W_tm[match(t, bm$t)]
  Ll <- bm$layers[lyr.idx,"Ld",drop=TRUE]
  Lu <- bm$layers[lyr.idx,"Lu",drop=TRUE]
  Ul <- bm$layers[lyr.idx,"Ud",drop=TRUE]
  Uu <- bm$layers[lyr.idx,"Uu",drop=TRUE]
  if(y == Ll) {
    minI <- +1
  } else if(y == Uu) {
    minI <- 0
    flip <- -c(x,y,Ll,Lu,Ul,Uu)
    x <- flip[1]
    y <- flip[2]
    Ll <- flip[6]
    Lu <- flip[5]
    Ul <- flip[4]
    Uu <- flip[3]
  } else {
    stop("Error: min/max mismatch at tau")
  }

  if(Uu-Ul > mult*2*(t-s)^(0.5)) {
    while(Uu-Ul > mult*2*(t-s)^(0.5)) { # Condition to determine whether refinement is necessary
      m <- 3 # Base for the truncation of alternating series
      Um <- (Ul+Uu)/2 # Mid point of upper intersection
      u <- runif(1,0,1) # Uniform draw to determine layer
      evalind <- 0
      while(evalind == 0) {
        LlUuprob <- eadelC_(m,s,t,x,y,Ll,Uu) # Probability associated with Ll->UU
        LlUmprob <- eadelC_(m,s,t,x,y,Ll,Um) # Probability associated with Ll->Um
        LlUlprob <- eadelC_(m,s,t,x,y,Ll,Ul) # Probability associated with Ll->Ul
        if(min(LlUuprob[1],LlUmprob[1],LlUlprob[1]) < 0) {
          m <- m + 2
        } else {
          UlUuprob <- LlUuprob - LlUlprob
          UlUmprob <- LlUmprob - LlUlprob
          if(u <= UlUmprob[1]/UlUuprob[2]) {
            evalind <- 1 # We have now refined
            Ul <- Ul
            Uu <- Um
          } # Update the layers
          if(u > UlUmprob[2]/UlUuprob[1]) {
            evalind <- 1 # We have now refined
            Ul <- Um
            Uu <- Uu
          } # Update the layers
          if(evalind == 0) { # Not refined enough so index
            m <- m + 2
          }
        }
      }
    }
  }

  if(Ul == x) {
    hard <- TRUE
  } else {
    hard <- FALSE
  }
  if(minI == 0) {
    flip <- -c(x,y,Ll,Lu,Ul,Uu)
    x <- flip[1]
    y <- flip[2]
    Ll <- flip[6]
    Lu <- flip[5]
    Ul <- flip[4]
    Uu <- flip[3]
    if(!hard) {
      bm$layers[lyr.idx,"Lu.hard"] <- FALSE
    } else {
      if(!bm$layers[lyr.idx,"Lu.hard"]) { # "Once soft you can't become hard again." - M Pollock, 6/1/2021
        stop("Error: Soft layer has become hard")
      }
      bm$layers[lyr.idx,"Lu.hard"] <- TRUE
    }
  } else {
    if(!hard) {
      bm$layers[lyr.idx,"Ud.hard"] <- FALSE
    } else {
      if(!bm$layers[lyr.idx,"Ud.hard"]) { # "Once soft you can't become hard again." - M Pollock, 6/1/2021
        stop("Error: Soft layer has become hard")
      }
      bm$layers[lyr.idx,"Ud.hard"] <- TRUE
    }
  }

  bm$layers[lyr.idx,"Ld"] <- Ll
  bm$layers[lyr.idx,"Lu"] <- Lu
  bm$layers[lyr.idx,"Ud"] <- Ul
  bm$layers[lyr.idx,"Uu"] <- Uu

  bm
}

refine.bblocal_ <- function(bm, lyr.idx, mult) {
  bblyr.idx <- which(bm$bb.local$layers$t.l == bm$layers$t.l[lyr.idx] &
                       bm$bb.local$layers$t.u == bm$layers$t.u[lyr.idx])

  if(bm$bb.local$layers$type[bblyr.idx] != "localised")
    stop("Auxilliary and main path layer mismatch detected.")

  refine.local_(bm$bb.local, bblyr.idx, mult)

  s_idx <- match(bm$layers$t.l[lyr.idx], bm$t)
  t_idx <- match(bm$layers$t.u[lyr.idx], bm$t)
  bb.s_idx <- match(bm$layers$t.l[lyr.idx], bm$bb.local$t)
  bb.t_idx <- match(bm$layers$t.u[lyr.idx], bm$bb.local$t)

  Bs <- bm$W_t[s_idx]
  Bt <- bm$W_tm[t_idx]
  Ws <- bm$bb.local$W_t[bb.s_idx]
  Wt <- bm$bb.local$W_t[bb.t_idx]

  # bm$layers[lyr.idx,"Ld"] <- Bs + bm$bb.local$layers[bblyr.idx,"Ld"] - Ws
  # bm$layers[lyr.idx,"Lu"] <- Bs + bm$bb.local$layers[bblyr.idx,"Ld"] - Ws + ((Bt-Bs)-(Wt-Ws))
  # bm$layers[lyr.idx,"Ud"] <- Bs + bm$bb.local$layers[bblyr.idx,"Uu"] - Ws
  # bm$layers[lyr.idx,"Uu"] <- Bs + bm$bb.local$layers[bblyr.idx,"Uu"] - Ws + ((Bt-Bs)-(Wt-Ws))
  bm$layers[lyr.idx,"Ld"] <- Bs - Ws
  bm$layers[lyr.idx,"Lu"] <- NA
  bm$layers[lyr.idx,"Ud"] <- NA
  bm$layers[lyr.idx,"Uu"] <- Bs - Ws + ((Bt-Bs)-(Wt-Ws))
  bm$layers[lyr.idx,"Lu.hard"] <- bm$bb.local$layers[bblyr.idx,"Lu.hard"]
  bm$layers[lyr.idx,"Ud.hard"] <- bm$bb.local$layers[bblyr.idx,"Ud.hard"]

  bm
}

refine.bessel_ <- function(bm, lyr.idx, mult) { # s, t, x, y, Ll, Lu, Ul, Uu, mult = 1) {
  s <- bm$layers[lyr.idx,"t.l",drop=TRUE]
  t <- bm$layers[lyr.idx,"t.u",drop=TRUE]
  x <- bm$W_t[match(s, bm$t)]
  y <- bm$W_tm[match(t, bm$t)]
  Ll <- bm$layers[lyr.idx,"Ld",drop=TRUE]
  Lu <- bm$layers[lyr.idx,"Lu",drop=TRUE]
  Ul <- bm$layers[lyr.idx,"Ud",drop=TRUE]
  Uu <- bm$layers[lyr.idx,"Uu",drop=TRUE]

  if(max((Uu-Ul),(Lu-Ll))>mult*((t-s)/2)^(0.5)){
    m.counter <- 3 # Initialise counter
    prob.min <- eagammaC_(m.counter,s,t,x,y,Lu,Ul) # Determine initial (outer) probabilities (min)
    prob.max <- eagammaC_(m.counter,s,t,x,y,Ll,Uu) # Determine initial (outer) probabilities (max)

    # Simulate the uniform random variable to determine layer
    repeat{
      u.refine <- runif(1,prob.min[1],prob.max[2]) # Simulate random variable to determine successor layer
      if(prob.min[1] == prob.max[2]){break} # Outwith tolerance for resolution - collapse onto interior ### WARNING MESSAGE NEEDED

      # Resolve outer minimum
      repeat{
        prob.min <- eagammaC_(m.counter,s,t,x,y,Lu,Ul) # Compute minimum outer probability
        if(u.refine <= prob.min[1] | u.refine > prob.min[2]){break} # If the uniform is sufficiently resolved then break
        m.counter <- m.counter + 2 # Else index counter
      }

      # Resolve outer maximum
      repeat{
        prob.max <- eagammaC_(m.counter,s,t,x,y,Ll,Uu) # Compute maximum outer probability
        if(u.refine <= prob.max[1] | u.refine > prob.max[2]){break} # If the uniform is sufficiently resolved then break
        m.counter <- m.counter + 2 # Else index counter
      }

      # Determine whether to resimulate uniform
      if(u.refine > prob.min[2] & u.refine <= prob.max[1]){break}
    }

    while(max((Uu-Ul),(Lu-Ll))>mult*((t-s)/2)^(0.5)){ # Criteria for refinement
      Ls<-(Ll+Lu)/2; Us<-(Ul+Uu)/2 # Compute the mid-points of the layers
      # Determine which layer refinement we fall into
      repeat{
        prob.mid <- eagammaC_(m.counter,s,t,x,y,Ls,Us) # Determine mid probability
        if(u.refine <= prob.mid[1]){Ll <- Ls; Lu <- Lu; Ul <- Ul; Uu <- Us; break} # If the uniform is sufficiently and we are in the interior redefine layers and then break
        if(u.refine > prob.mid[2]){Ll <- Ll; Lu <- Ls; Ul <- Us; Uu <- Uu; break} # If the uniform is sufficiently and we are in the outer then redefine layers and then break
        m.counter <- m.counter + 2 # Else index counter
      }
    }
  }

  ### Output
  #list(layer = c(Ll,Lu,Ul,Uu))

  bm$layers[lyr.idx,"Ld"] <- Ll
  bm$layers[lyr.idx,"Lu"] <- Lu
  bm$layers[lyr.idx,"Ud"] <- Ul
  bm$layers[lyr.idx,"Uu"] <- Uu
  bm$layers[lyr.idx,"Lu.hard"] <- ifelse(Lu == min(x,y), TRUE, FALSE)
  bm$layers[lyr.idx,"Ud.hard"] <- ifelse(Ul == max(x,y), TRUE, FALSE)

  bm
}

refine.bbbessel_ <- function(bm, lyr.idx, mult) {
  bblyr.idx <- which(bm$bb.local$layers$t.l == bm$layers$t.l[lyr.idx] &
                       bm$bb.local$layers$t.u == bm$layers$t.u[lyr.idx])

  if(bm$bb.local$layers$type[bblyr.idx] != "bessel")
    stop("Auxilliary and main path layer mismatch detected.")

  refine.bessel_(bm$bb.local, bblyr.idx, mult)

  s_idx <- match(bm$layers$t.l[lyr.idx], bm$t)
  t_idx <- match(bm$layers$t.u[lyr.idx], bm$t)
  bb.s_idx <- match(bm$layers$t.l[lyr.idx], bm$bb.local$t)
  bb.t_idx <- match(bm$layers$t.u[lyr.idx], bm$bb.local$t)

  Bs <- bm$W_t[s_idx]
  Bt <- bm$W_tm[t_idx]
  Ws <- bm$bb.local$W_t[bb.s_idx]
  Wt <- bm$bb.local$W_t[bb.t_idx]

  # bm$layers[lyr.idx,"Ld"] <- Bs + bm$bb.local$layers[bblyr.idx,"Ld"] - Ws
  # bm$layers[lyr.idx,"Lu"] <- Bs + bm$bb.local$layers[bblyr.idx,"Ld"] - Ws + ((Bt-Bs)-(Wt-Ws))
  # bm$layers[lyr.idx,"Ud"] <- Bs + bm$bb.local$layers[bblyr.idx,"Uu"] - Ws
  # bm$layers[lyr.idx,"Uu"] <- Bs + bm$bb.local$layers[bblyr.idx,"Uu"] - Ws + ((Bt-Bs)-(Wt-Ws))
  bm$layers[lyr.idx,"Ld"] <- Bs - Ws
  bm$layers[lyr.idx,"Lu"] <- NA
  bm$layers[lyr.idx,"Ud"] <- NA
  bm$layers[lyr.idx,"Uu"] <- Bs - Ws + ((Bt-Bs)-(Wt-Ws))
  bm$layers[lyr.idx,"Lu.hard"] <- bm$bb.local$layers[bblyr.idx,"Lu.hard"]
  bm$layers[lyr.idx,"Ud.hard"] <- bm$bb.local$layers[bblyr.idx,"Ud.hard"]

  bm
}
