#' Construct Bessel layers
#'
#' Performs unbiased simulation of the smallest Bessel layers which contain the
#' Brownian motion between two known sample points by retrospective Bernoulli
#' sampling and inversion sampling.
#'
#' @param bm a Brownian motion object from which simulation should continue.
#'   Note the object is updated in place
#' @param s left hand time point
#' @param t right hand time point
#' @param mult the default layer size is \code{sqrt(t-s)}.  You can scale
#'   this by specifying the \code{mult} argument which will result in layer
#'   sizes of \code{mult*sqrt(t-s)}.
#'
#' @return
#' The Brownian motion object which was passed in argument \code{bm} is
#' updated in place and returned, enabling chaining of commands with
#' dplyr (and other) style pipes.
#'
#' @export
bessel.layers <- function(bm, s, t, refine = bm$refine, mult = bm$mult, prefer = bm$prefer) {
  ## NOTE TO LOUIS: This shares a lot of setup with Bessel Layers -- pull out into utility code (except the find Bessel layers part)
  if(!("BrownianMotion" %in% class(bm))) {
    stop("bm argument must be a BrownianMotion object.")
  }
  if(!is.realscalar(s) || !is.realscalar(t) || !is.realscalar(mult)) {
    stop("s, t and mult must be real scalars.")
  }
  if(s >= t) {
    stop("s must be less than t.")
  }
  if(mult <= 0) {
    stop("mult must be strictly positive.")
  }

  # Are these endpoints part of the skeleton?
  # If not, and the point is beyond the end, simulate it.  Otherwise, error
  if(!(s %in% bm$t)) {
    sim(bm, s, refine, mult, prefer)
  }
  if(!(t %in% bm$t)) {
    sim(bm, t, refine, mult, prefer)
  }

  # Make sure no part of [s,t] interval lies inside any other layer
  # if(any(bm$bounds$t.l>s & bm$bounds$t.u<t) || any(bm$bessel.layers$t.l>s & bm$bessel.layers$t.u<t)) {
  #   stop("Cannot create bessel layer within another layer specification.")
  # }

  # Find all intervals between s and t not currently inside any other layer or inside an Intersection layer
  s.i <- which(s == bm$t)
  t.i <- which(t == bm$t)
  all.pairs.noL <- matrix(c(bm$t[s.i:(t.i-1)], bm$t[(s.i+1):t.i]), byrow = FALSE, ncol = 2)
  all.pairs.IL  <- matrix(c(bm$t[s.i:(t.i-1)], bm$t[(s.i+1):t.i]), byrow = FALSE, ncol = 2)
  # No layer
  mid <- (all.pairs.noL[,2]+all.pairs.noL[,1])/2
  incl <- rep(TRUE, nrow(all.pairs.noL))
  for(i in 1:length(mid)) {
    if(any(mid[i] >= bm$layers$t.l & mid[i] <= bm$layers$t.u) ||
       any(mid[i] >= bm$user.layers$t.l & mid[i] <= bm$user.layers$t.u)) {
      incl[i] <- FALSE
    }
  }
  all.pairs.noL <- all.pairs.noL[incl,,drop = FALSE]
  # cat("No layer:"); print(all.pairs.noL)
  # Intersection layer
  mid <- (all.pairs.IL[,2]+all.pairs.IL[,1])/2
  incl <- rep(FALSE, nrow(all.pairs.IL))
  all.pairs.IL.lyr_idx <- c()
  for(i in 1:length(mid)) {
    if(any(mid[i] >= bm$layers$t.l & mid[i] <= bm$layers$t.u & bm$layers$type == "intersection")) {
      incl[i] <- TRUE
      all.pairs.IL.lyr_idx <- c(all.pairs.IL.lyr_idx, which(mid[i] >= bm$layers$t.l & mid[i] <= bm$layers$t.u & bm$layers$type == "intersection"))
    }
  }
  all.pairs.IL <- all.pairs.IL[incl,,drop = FALSE]
  # cat("Intersection layer:"); print(all.pairs.IL)

  if(nrow(all.pairs.noL) > 0) {
    for(i in 1:nrow(all.pairs.noL)) {
      bessel.layers_(bm, all.pairs.noL[i,1], all.pairs.noL[i,2], mult)
    }
  }
  if(length(all.pairs.IL.lyr_idx) > 0) {
    intersection.to.bessel_(bm, all.pairs.IL.lyr_idx)
  }

  bm$layers <- bm$layers[order(bm$layers$t.l),]
  invisible(bm)
}

bessel.layers_ <- function(bm, s, t, mult) {
  res <- bessel.layers.sim_(bm, s, t, bm$W_t[match(s, bm$t)], bm$W_t[match(t, bm$t)], mult)

  bm$layers <- add_row(bm$layers,
                       type = "bessel",
                       t.l = s,
                       t.u = t,
                       Ld = res$xb - res$au,
                       Uu = res$yb + res$au,
                       Lu = res$xb - res$al,
                       Ud = res$yb + res$al,
                       Lu.hard = ifelse(res$act == 1, TRUE, FALSE),
                       Ud.hard = ifelse(res$act == 1, TRUE, FALSE))

  bm
}

bessel.layers.sim_ <- function(bm, s, t, x, y, mult) {
  xb <- min(x,y)
  yb <- max(x,y)
  adf <- sqrt(t-s) * mult
  act <- 1
  gind <- 0
  u1 <- runif(1,0,1)
  m1 <- 3
  il <- eagammaC_(m1,s,t,x,y,xb-adf*act,yb+adf*act)
  while(gind == 0) {
    if(u1 >= il[2]) {
      act <- act+1
      m1 <- 3
      il <- eagammaC_(m1,s,t,x,y,xb-adf*act,yb+adf*act)
    } else {
      if(u1 <= il[1]) {
        gind <- 1
      } else {
        m1 <- m1+2
        il <- eagammaC_(m1,s,t,x,y,xb-adf*act,yb+adf*act)
      }
    }
  }
  list(xb = xb, yb = yb, adf = adf, act = act, al = (act-1)*adf, au = act*adf)
}

