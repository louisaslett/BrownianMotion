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
bessel.layers <- function(bm, s, t, refine = bm$refine, mult = bm$mult, prefer = bm$prefer, label = c(names(s), names(t))) {
  new.layers(bm, s, t, "bessel", refine, mult, prefer, label)
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

