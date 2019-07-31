#' Construct Bessel layers
#'
#' Unbiasedly simulates the smallest Bessel layers which contains the Brownian
#' motion between two known sample points by retrospective Bernoulli sampling
#' and inversion sampling.
#'
#' @param bm a Brownian motion object from which simulation should continue.
#'   Note the object is updated in place.
#' @param s left hand time point
#' @param t right hand time point
#'
#' @return
#' The Brownian motion object which was passed in argument \code{bm} is
#' updated in place and returned, enabling chaining of commands with
#' dplyr (and other) style pipes.
#'
#' @export
bessel.layers <- function(bm, s, t) {
  if(!("BrownianMotion" %in% class(bm))) {
    stop("bm argument must be a BrownianMotion object.")
  }
  if(!is.realscalar(s) || !is.realscalar(t)) {
    stop("s and t must be real scalars.")
  }
  if(!(s %in% bm$t)) {
    stop("s is not a known time point in the supplied Brownian motion object.")
  }
  if(!(t %in% bm$t)) {
    stop("t is not a known time point in the supplied Brownian motion object.")
  }

  invisible(bessel.layers_(bm, s, t))
}

bessel.layers_ <- function(bm, s, t) {
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

bessel.layers.sim_ <- function(bm, s, t, x, y) {
  xb <- min(x,y)
  yb <- max(x,y)
  adf <- sqrt((t-s)/2)
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

