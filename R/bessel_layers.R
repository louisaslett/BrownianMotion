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
#' @param mult the default layer size is \code{sqrt((t-s)/2)}.  You can scale
#'   this by specifying the \code{mult} argument which will result in layer
#'   sizes of \code{mult*sqrt((t-s)/2)}.
#'
#' @return
#' The Brownian motion object which was passed in argument \code{bm} is
#' updated in place and returned, enabling chaining of commands with
#' dplyr (and other) style pipes.
#'
#' @export
bessel.layers <- function(bm, s, t, mult = 1) {
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
  if(!(s %in% bm$t)) {
    stop("s is not a known time point in the supplied Brownian motion object.")
  }
  if(!(t %in% bm$t)) {
    stop("t is not a known time point in the supplied Brownian motion object.")
  }

  # Make sure no part of [s,t] interval lies inside any other layer
  # if(any(bm$bounds$t.l>s & bm$bounds$t.u<t) || any(bm$bessel.layers$t.l>s & bm$bessel.layers$t.u<t)) {
  #   stop("Cannot create bessel layer within another layer specification.")
  # }

  # Find all intervals between s and t not currently inside any other layer
  s.i <- which(s == bm$t)
  t.i <- which(t == bm$t)
  all.pairs <- matrix(c(bm$t[s.i:(t.i-1)], bm$t[(s.i+1):t.i]), byrow = FALSE, ncol = 2)
  mid <- (all.pairs[,2]+all.pairs[,1])/2
  incl <- rep(TRUE, nrow(all.pairs))
  for(i in 1:length(mid)) {
    if(any(mid[i] >= bm$bounds$t.l & mid[i] <= bm$bounds$t.u) ||
       any(mid[i] >= bm$bessel.layers$t.l & mid[i] <= bm$bessel.layers$t.u)) {
      incl[i] <- FALSE
    }
  }
  all.pairs <- all.pairs[incl,,drop = FALSE]

  if(nrow(all.pairs) == 0) {
    stop("No intervals in the skeleton between s and t found that do not have a layer specification already.")
  }

  for(i in 1:nrow(all.pairs)) {
    bessel.layers_(bm, all.pairs[i,1], all.pairs[i,2], mult)
  }

  invisible(bm)
}

bessel.layers_ <- function(bm, s, t, mult) {
  res <- bessel.layers.sim_(bm, s, t, bm$W_t[match(s, bm$t)], bm$W_t[match(t, bm$t)])

  bm$bessel.layers <- add_row(bm$bessel.layers,
                              t.l = s,
                              t.u = t,
                              l = res$xb - res$au,
                              u = res$yb + res$au,
                              L = res$xb - res$al,
                              U = res$yb + res$al)

  bm
}

bessel.layers.sim_ <- function(bm, s, t, x, y, mult = 1) {
  xb <- min(x,y)
  yb <- max(x,y)
  adf <- sqrt((t-s)/2) * mult
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

