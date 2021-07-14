#' Find a first passage time and state
#'
#' Finds the first passage time of a Brownian motion object to a specified lower,
#' upper, or lower and upper limit.
#'
#' This function updates a Brownian motion object adding the first passage time to
#' the specified lower (\code{l}), upper (\code{u}), or lower and upper limits.
#' The limits can also be specfied via \code{delta} which places the limits as the
#' current endpoint state +/- \code{delta}.
#' Note that simulation of the first passage implicitly imposes constraints on the
#' underlying Brownian motion between the endpoint and the first passage
#' time.
#'
#' In particular, in the case that the lower and upper limits are asymmetric around
#' the endpoint, there may be a recursive simulation of nested symmetric first passages
#' between the endpoint and the ultimate first passage, resulting in implicit imposition
#' of nested constraints.  When plotting a Brownian motion these are visible as dashed
#' red lines nested within the solid red line constraint imposed by the overall first
#' passage simulation.
#'
#' @param bm a Brownian motion object from which simulation should continue.
#'   Note the object is updated in place.
#' @param l scalar giving lower bound for first passage.
#' @param u scalar giving upper bound for first passage.
#' @param delta scalar providing the bounds by way of the positive/negative
#'   change in the Brownian motion.
#'
#' @return
#' The Brownian motion object which was passed in argument \code{bm} is
#' updated in place and returned, enabling chaining of commands with
#' dplyr (and other) style pipes.
#'
#' @export
first.passage <- function(bm, l = NULL, u = NULL, delta.l = NULL, delta.u = NULL, delta = NULL, label = NULL) {
  if(!("BrownianMotion" %in% class(bm))) {
    stop("bm argument must be a BrownianMotion object.")
  }
  # Check for a valid combination of arguments
  if(!is.null(delta) && (!is.null(l) || !is.null(u) || !is.null(delta.l) || !is.null(delta.u))) {
    stop("if delta is specified, the other limit arguments must not be specified.")
  }
  if((!is.null(delta.l) || !is.null(delta.u)) && (!is.null(l) || !is.null(u) || !is.null(delta))) {
    stop("if delta.l or delta.u are specified, the other limit arguments must not be specified.")
  }
  if((!is.null(l) || !is.null(u)) && (!is.null(delta.l) || !is.null(delta.u) || !is.null(delta))) {
    stop("if delta.l or delta.u are specified, the other limit arguments must not be specified.")
  }
  if(is.null(l) && is.null(u) && is.null(delta.l) && is.null(delta.u) && is.null(delta)) {
    # default to 1 if no limits specified
    delta <- 1
  }
  assert.bmlabel(label, 1)

  # Convert every combination to an l&u pair
  x <- state(bm)
  if(!is.null(delta)) {
    if(!is.realscalar(delta) || delta <= 0) {
      stop("if delta is specified, it must be a positive real scalar.")
    }
    l <- x - delta
    u <- x + delta
  }
  if(!is.null(delta.l)) {
    if(!is.realscalar(delta.l) || delta.l >= 0) {
      stop("if delta.l is specified, it must be a negative real scalar.")
    }
    l <- x + delta.l
  }
  if(!is.null(delta.u)) {
    if(!is.realscalar(delta.u) || delta.u <= 0) {
      stop("if delta.u is specified, it must be a positive real scalar.")
    }
    u <- x + delta.u
  }
  if(!is.realscalar(l)) {
    stop("l must specify a scalar lower first passage bound.")
  }
  if(!is.realscalar(u)) {
    stop("l must specify a scalar upper first passage bound.")
  }
  if(l >= u) {
    stop("l cannot be greater than or equal to u.")
  }
  if(l >= x) {
    stop("l must be less than the current endpoint state of the Brownian motion.")
  }
  if(u <= x) {
    stop("u must be greater than the current endpoint state of the Brownian motion.")
  }

  invisible(first.passage_(bm, l, u, label))
}

first.passage_ <- function(bm, l, u, label) {
  # Get current time/position
  t <- max(bm$t)
  x <- bm$W_t[match(t, bm$t)]

  # Storage for extension of path
  # NB x's will be zero centred and adjusted before appending to bm$W_t at the end
  t.new <- c()
  x.new <- c()

  t.l.new <- c()
  t.u.new <- c()
  L.new <- c()
  U.new <- c()

  # Transform limits to zero centred
  l <- l - x
  u <- u - x

  if(isTRUE(all.equal(abs(l), u))) {
    # Symmetric interval, so go straight to simulate passage time + location and done
    t.new <- t + stdbm.first.passage_()*u^2
    x.new <- sample(c(l, u), 1)

    t.l.new <- t
    t.u.new <- t.new
    L.new <- l
    U.new <- u
  } else {
    # Asymmetry in the limits so need to do successive sims
    delta <- min(abs(l), u)
# cat(glue::glue("Asymmetry detected, choosing a delta of {delta}\n\n"))
    # Preempt lower bound and new path for next iteration
    t.l.new <- t
    x.new <- 0

    # Loop until BM reaches either upper or lower boundary
    while(!isTRUE(all.equal(l, tail(x.new, 1), check.attributes = FALSE)) && !isTRUE(all.equal(u, tail(x.new, 1), check.attributes = FALSE))) {
      t.new2 <- tail(t.l.new, 1) + stdbm.first.passage_()*delta^2
      x.new2 <- sample(c(+delta, -delta), 1) + tail(x.new, 1)
# cat(glue::glue("Hitting time {t.new2}, at {x.new2}.\n\n"))

      t.u.new <- c(t.u.new, t.new2)
      L.new <- c(L.new, tail(x.new, 1) - delta)
      U.new <- c(U.new, tail(x.new, 1) + delta)
      # Preempt lower bound for next iteration
      t.l.new <- c(t.l.new, t.new2)

      # NB in >1D make sure to do this more efficiently since we are going to always add at least d
      t.new <- c(t.new, t.new2)
      x.new <- c(x.new, x.new2)

      delta <- min(x.new2-l, u-x.new2)
# cat(glue::glue("Next delta {delta}.\n\n"))
# if(is.infinite(delta)) stop()
    }
    # Once we're out, we'll have one too many t.l.new and x.new due to prempting next iteration, so
    # remove the trailing/heading one
    t.l.new <- head(t.l.new, -1)
    x.new <- tail(x.new, -1)
  }

  # Now append the new sims to the BM object
  bm$t <- c(bm$t, unname(t.new))
  bm$W_t <- c(bm$W_t, x.new+x)
  bm$W_tm <- c(bm$W_tm, x.new+x)
  bm$layers <- add_row(bm$layers,
                       type = "localised",
                       t.l = t.l.new,
                       t.u = t.u.new,
                       Ld = L.new+x,
                       Uu = U.new+x,
                       Lu = pmin(head(c(x, x.new+x), -1), x.new+x),
                       Ud = pmax(head(c(x, x.new+x), -1), x.new+x),
                       Lu.hard = TRUE, #ifelse(head(c(0, x.new), -1) > x.new, TRUE, FALSE)
                       Ud.hard = TRUE)
  bm$user.layers <- add_row(bm$user.layers,
                            t.l = min(t.l.new),
                            t.u = max(t.u.new),
                            L = l+x,
                            U = u+x)

  add.labels_(bm, "fpt", t.new)
  add.labels_(bm, "internal", t.new)
  add.labels_(bm, label, t.new)
  bm$labels[["end"]] <- max(t.new)
  bm$labels[["seg.end"]][length(bm$labels[["seg.end"]])] <- bm$labels[["end"]]

  bm
}

# p <- 4/pi*exp(-pi^2*0.64/8)
# p <- ceiling(p*1e6)/1e6
# q <- 4*pnorm(1/sqrt(0.64), lower = FALSE)
# q <- ceiling(q*1e6)/1e6
# print(p, digits = 22)
# print(q, digits = 22)
# print(p+q, digits = 22)
# print(p/(p+q), digits = 22)
# print(1-q/(p+q), digits = 22)
# identical(p/(p+q), 1-q/(p+q))

# Adapted from Scale.R (git:mpoll/scale:93fcd7f) dev.pr:159
stdbm.first.passage.proposal_ <- function() {
  p.over.pplusq <- 0.5776968790939969178311

  U <- runif(1)
  if(U < p.over.pplusq) {
    return(0.64 + 8*rexp(1)/pi^2)
  } else {
    repeat {
      E <- rexp(2)
      if(E[1]^2 <= 2*E[2]/0.64) {
        return(0.64/(1+0.64*E[1])^2)
      }
    }
  }
}

# Adapted from Scale.R (git:mpoll/scale:93fcd7f) dev.ref:167
stdbm.first.passage_ <- function() {
  accept <- FALSE
  while(!accept) {
    t <- stdbm.first.passage.proposal_()
    if(t <= 0.64) { # Series sampler 1, used around critical point 0.64 (t*<= 0.64, see Appendix C ScaLE)
      S <- exp(-1/(2*t))/2
      n <- -1
      Y <- runif(1)*S
      repeat {
        n <- n+2 # Indexation of sequence for subsequent odd / even terms
        S <- S - (n+0.5)*exp(-2*(n+0.5)^2/t)
        if(Y <= S) {
          accept <- TRUE
          break
        }
        S <- S + (n+1.5)*exp(-2*(n+1.5)^2/t)
        if(Y >= S) {
          break
        }
      }
    } else { # Series sampler 2, used around critical point 0.64 (t*> 0.64, see Appendix C  ScaLE)
      S <- exp(-pi^2*t/8)/2
      n <- -1
      Y <- runif(1,0,1)*S
      repeat {
        n <- n+2 # Indexation of sequence for subsequent odd / even terms
        S <- S - (n+0.5)*exp(-(n+0.5)^2*pi^2*t/2)
        if(Y <= S) {
          accept <- TRUE
          break
        }
        S <- S + (n+1.5)*exp(-(n+1.5)^2*pi^2*t/2)
        if(Y >= S) {
          break
        }
      }
    }
  }
  t
}
