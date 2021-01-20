#' Simulate Brownian motion at a target time
#'
#' Simulates Brownian motion at chosen times, automatically applying the correct
#' algorithm to condition on the current known (ie previously simulated) states
#' of the trajectory.
#'
#' @param bm a Brownian motion object from which simulation should continue.
#'   Note the object is updated in place.
#' @param t a vector of times to simulate at.
#' @param refine whether to automatically refine layers where the simulated
#'   time results in bisection of an existing layer.  Defaults to \code{FALSE}.
#'
#' @return the Brownian motion object which was passed in argument \code{bm} is
#'   updated in place and returned, enabling chaining of commands with
#'   dplyr (and other) style pipes.
#'
#' @export
sim <- function(bm, t, refine = FALSE) {
  if(!("BrownianMotion" %in% class(bm))) {
    stop("bm argument must be a BrownianMotion object.")
  }
  if(!is.realvector(t)) {
    stop("t must be a vector of times.")
  }
  if(is.unsorted(t)) {
    t <- sort(t)
  }
  if(!isFALSE(refine) && refine <= 0) {
    stop("refine must be FALSE (to disable refinement), TRUE (to do automatic refinement), or a strictly positive value (to supply the mult parameter during refinement).")
  }

  # Eliminate times we know
  t <- setdiff(t, bm$t)

  # Find times which are standard forward simulation at the end
  t.bmfwd <- t[which(t>tail(bm$t, 1))]
  # Find Brownian bridge candidate times
  # Reverse sort them so that index in bm$t is correct over many iterations
  t.bmbb <- sort(t[which(t<tail(bm$t, 1))], decreasing = TRUE)
  t.local <- c()
  t.intersection <- c()
  t.bessel <- c()

  # Do forward sims
  if(length(t.bmfwd) > 0) {
    sim.uncond_(bm, t.bmfwd)
  }
  # Do BB sims
  if(length(t.bmbb) > 0) {
    # Check if any are in bounded regions
    if(nrow(bm$layers) > 0) {
      t.local <- t.bmbb[colSums(matrix(apply(bm$layers[bm$layers$type == "localised",2:3], 1, function(y) { t.bmbb > y[1] & t.bmbb < y[2] }), ncol = length(t.bmbb), byrow = TRUE)) > 0]
      # Keep rest as standard Brownian bridges
      t.bmbb <- setdiff(t.bmbb, t.local)

      t.intersection <- t.bmbb[colSums(matrix(apply(bm$layers[bm$layers$type == "intersection",2:3], 1, function(y) { t.bmbb > y[1] & t.bmbb < y[2] }), ncol = length(t.bmbb), byrow = TRUE)) > 0]
      # Keep rest as standard Brownian bridges
      t.bmbb <- setdiff(t.bmbb, t.intersection)

      t.bessel <- t.bmbb[colSums(matrix(apply(bm$layers[bm$layers$type == "bessel",2:3], 1, function(y) { t.bmbb > y[1] & t.bmbb < y[2] }), ncol = length(t.bmbb), byrow = TRUE)) > 0]
      # Keep rest as standard Brownian bridges
      t.bmbb <- setdiff(t.bmbb, t.bessel)
    }
    if(length(t.local) > 0) {
      # Do simulation conditional on localised layers
      t.local <- sapply(t.local,
                        function(t) {
                          c(l = tail(which(bm$t<t), 1), r = head(which(bm$t>t), 1), q = t)
                        })
      for(l in unique(t.local[1,])) {
        i <- which(t.local[1,]==l)
        r <- t.local[2,i[1]]
        if(any(!bm$layers[match(bm$t[r], bm$layers$t.u),c("Lu.hard","Ud.hard")])) {
          warning(glue("Cannot currently simulate conditional on soft localised layer (at time {t.local[3,i[1]]}), skipping ...\n"))
          next
        }
        for(qq in rev(t.local[3,i])) {
          bm.res <- sim.condlocal_(bm, l, qq, r)
          l <- l+1
          r <- r+1
          if(!isFALSE(refine)) {
            refine.intersection_(bm, match(qq, bm$layers$t.u), refine)
            refine.local_(bm, match(qq, bm$layers$t.l), refine)
          }
        }
      }
    }
    if(length(t.intersection) > 0) {
      # Do simulation conditional on intersection layers
      t.intersection <- sapply(t.intersection,
                               function(t) {
                                 c(l = tail(which(bm$t<t), 1), r = head(which(bm$t>t), 1), q = t)
                               })
      for(l in unique(t.intersection[1,])) {
        i <- which(t.intersection[1,]==l)
        r <- t.intersection[2,i[1]]
        for(qq in rev(t.intersection[3,i])) {
          bm.res <- sim.condintersection_(bm, l, qq, r)
          l <- l+1
          r <- r+1
          if(!isFALSE(refine)) {
            refine.intersection_(bm, match(qq, bm$layers$t.u), refine)
            refine.intersection_(bm, match(qq, bm$layers$t.l), refine)
          }
        }
      }
    }
    if(length(t.bessel) > 0) {
      # Do simulation conditional on intersection layers
      t.bessel <- sapply(t.bessel,
                         function(t) {
                           c(l = tail(which(bm$t<t), 1), r = head(which(bm$t>t), 1), q = t)
                         })
      for(l in unique(t.bessel[1,])) {
        i <- which(t.bessel[1,]==l)
        r <- t.bessel[2,i[1]]
        for(qq in rev(t.bessel[3,i])) {
          bm.res <- sim.condbessel_(bm, l, qq, r)
          l <- l+1
          r <- r+1
          if(!isFALSE(refine)) {
            refine.intersection_(bm, match(qq, bm$layers$t.u), refine)
            refine.intersection_(bm, match(qq, bm$layers$t.l), refine)
          }
        }
      }
    }
    if(length(t.bmbb) > 0) {
      # Do standard Brownian bridges
      t.bmbb <- sapply(t.bmbb,
                       function(t) {
                         c(l = tail(which(bm$t<t), 1), r = head(which(bm$t>t), 1), q = t)
                       })
      for(l in unique(t.bmbb[1,])) {
        i <- which(t.bmbb[1,]==l)
        r <- t.bmbb[2,i[1]]
        sim.bb_(bm, l, c(bm$t[l], rev(t.bmbb[3,i]), bm$t[r]), r)
      }
    }
  }

  invisible(bm)
}
