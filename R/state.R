#' Return process state
#'
#' Returns the state of the Brownian motion at a given time.
#'
#' @param bm a Brownian motion object from which simulation should continue.
#'   Note the object is updated in place.
#' @param t a scalar time at which the state is required.  This time must have
#'   already been simulated and be present in the Brownian motion object
#'   supplied in the \code{bm} argument.  If omitted, the most recent state is
#'   returned.
#'
#' @return a vector containing the state at the specified time.
#'
#' @export
state <- function(bm, t = NULL) {
  if(!("BrownianMotion" %in% class(bm))) {
    stop("bm argument must be a BrownianMotion object.")
  }
  if(!is.null(t) && !is.realvector(t)) {
    stop("t must be either NULL or a scalar time.")
  }

  if(is.null(t)) {
    t <- max(bm$t)
  }

  invisible(state_(bm, t))
}

state_ <- function(bm, t) {
  t_idx <- match(t, bm$t)
  return(bm$W_t[t_idx])

  if(any(t >= bm$layers$t.l & t <= bm$layers$t.u)) {
    l <- which(t >= bm$layers$t.l & t <= bm$layers)
    if(bm$layers$type[l] == "localised-bb") {
      ss <- bm$layers$t.l[l]
      tt <- bm$layers$t.u[l]
      x <- bm$bb.local$x[l]
      y <- bm$bb.local$y[l]
      l <- which(t >= bm$bb.local$t.l & t <= bm$bb.local)
      list(t = t,
           W_t = bm$W_t[t_idx],
           in.layer = "localised",
           t.l = bm$bb.local$t.l[l],
           t.u = bm$bb.local$t.u[l],
           Ld = bm$bb.local$Ld[l] + new.bm$t[i]/(t-s)*((y-x)-W_T) + x,
           Uu = bm$bb.local$Uu[l],
           Lu = bm$bb.local$Lu[l],
           Ud = bm$bb.local$Ud[l],
           Lu.hard = bm$bb.local$Lu.hard[l],
           Ud.hard = bm$bb.local$Ud.hard[l])
    } else {
      list(t = t,
           W_t = bm$W_t[t_idx],
           in.layer = bm$layers$type[l],
           t.l = bm$layers$t.l[l],
           t.u = bm$layers$t.u[l],
           Ld = bm$layers$Ld[l],
           Uu = bm$layers$Uu[l],
           Lu = bm$layers$Lu[l],
           Ud = bm$layers$Ud[l],
           Lu.hard = bm$layers$Lu.hard[l],
           Ud.hard = bm$layers$Ud.hard[l])
    }
  } else {
    list(t = t,
         W_t = bm$W_t[t_idx])
  }
}
