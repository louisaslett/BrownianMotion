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
#'   time results in bisection of an existing layer.  Defaults to the option
#'   specified at the time of creation of \code{bm}.
#' @param mult the mult parameter to be passed through to any layer operations
#'   including for the level of refinement if \code{refine=TRUE}.   Defaults to the option
#'   specified at the time of creation of \code{bm}.
#' @param prefer indicate the preferred layer type which should arise after
#'   conditional simulation bisects an existing layer (this will be respected
#'   where possible, but it not always achievable).    Defaults to the option
#'   specified at the time of creation of \code{bm}.
#'
#' @return the Brownian motion object which was passed in argument \code{bm} is
#'   updated in place and returned, enabling chaining of commands with
#'   dplyr (and other) style pipes.
#'
#' @export
sim <- function(bm, t, refine = bm$refine, mult = bm$mult, prefer = bm$prefer) {
  if(!("BrownianMotion" %in% class(bm))) {
    stop("bm argument must be a BrownianMotion object.")
  }
  if(!is.realvector(t)) {
    stop("t must be a vector of times.")
  }
  if(is.unsorted(t)) {
    t <- sort(t)
  }
  if(length(refine) != 1 || !is.logical(refine)) {
    stop("refine must be a scalar logical value")
  }
  if(length(mult) != 1 || !is.numeric(mult) || mult <= 0) {
    stop("mult must be a strictly positive scalar")
  }
  if(length(prefer) != 1 || !(prefer %in% c("bessel", "intersection"))) {
    stop("prefer must be one of 'bessel' or 'intersection'")
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
  t.bblocal <- c()
  t.bbintersection <- c()
  t.bbbessel <- c()

  # Do forward sims
  if(length(t.bmfwd) > 0) {
    sim.uncond_(bm, t.bmfwd)
  }
  # Do BB sims
  if(length(t.bmbb) > 0) {
    # Check if any are in bounded regions
    if(nrow(bm$layers) > 0) {
      lyr.idxs <- which(apply(bm$layers[,2:3], 1, function(y) { t.bmbb > y[1] & t.bmbb < y[2] }))
      if(prefer == "bessel")
        intersection.to.bessel_(bm, lyr.idxs)
      else
        bessel.to.intersection_(bm, lyr.idxs, mult)

      t.local <- t.bmbb[colSums(matrix(apply(bm$layers[bm$layers$type == "localised",2:3], 1, function(y) { t.bmbb > y[1] & t.bmbb < y[2] }), ncol = length(t.bmbb), byrow = TRUE)) > 0]
      # Keep rest as standard Brownian bridges
      t.bmbb <- setdiff(t.bmbb, t.local)

      t.intersection <- t.bmbb[colSums(matrix(apply(bm$layers[bm$layers$type == "intersection",2:3], 1, function(y) { t.bmbb > y[1] & t.bmbb < y[2] }), ncol = length(t.bmbb), byrow = TRUE)) > 0]
      # Keep rest as standard Brownian bridges
      t.bmbb <- setdiff(t.bmbb, t.intersection)

      t.bessel <- t.bmbb[colSums(matrix(apply(bm$layers[bm$layers$type == "bessel",2:3], 1, function(y) { t.bmbb > y[1] & t.bmbb < y[2] }), ncol = length(t.bmbb), byrow = TRUE)) > 0]
      # Keep rest as standard Brownian bridges
      t.bmbb <- setdiff(t.bmbb, t.bessel)

      t.bblocal <- t.bmbb[colSums(matrix(apply(bm$layers[bm$layers$type == "localised-bb",2:3], 1, function(y) { t.bmbb > y[1] & t.bmbb < y[2] }), ncol = length(t.bmbb), byrow = TRUE)) > 0]
      # Keep rest as standard Brownian bridges
      t.bmbb <- setdiff(t.bmbb, t.bblocal)

      t.bbintersection <- t.bmbb[colSums(matrix(apply(bm$layers[bm$layers$type == "intersection-bb",2:3], 1, function(y) { t.bmbb > y[1] & t.bmbb < y[2] }), ncol = length(t.bmbb), byrow = TRUE)) > 0]
      # Keep rest as standard Brownian bridges
      t.bmbb <- setdiff(t.bmbb, t.bbintersection)

      t.bbbessel <- t.bmbb[colSums(matrix(apply(bm$layers[bm$layers$type == "bessel-bb",2:3], 1, function(y) { t.bmbb > y[1] & t.bmbb < y[2] }), ncol = length(t.bmbb), byrow = TRUE)) > 0]
      # Keep rest as standard Brownian bridges
      t.bmbb <- setdiff(t.bmbb, t.bbbessel)
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
          soft <- TRUE
        } else {
          soft <- FALSE
        }
        for(qq in rev(t.local[3,i])) {
          if(soft) {
            sim.condlocalsoft_(bm, l, qq, r)
          } else {
            sim.condlocal_(bm, l, qq, r)
          }
          if(prefer == "bessel") {
            intersection.to.bessel_(bm, match(qq, bm$layers$t.u))
          }
          l <- l+1
          r <- r+1
          if(!isFALSE(refine)) {
            refine_(bm, match(qq, bm$layers$t.u), mult)
            refine_(bm, match(qq, bm$layers$t.l), mult)
          }
          soft <- any(!tail(bm$layers)[,c("Lu.hard","Ud.hard")])
        }
      }
    }
    while(length(t.intersection) > 0) {
      # Do simulation conditional on intersection layers
      lrq <- sapply(t.intersection,
                    function(t) {
                      c(l = tail(which(bm$t<t), 1), r = head(which(bm$t>t), 1), q = t)
                    })
      i <- which.min(abs(abs(lrq[3,]-bm$t[lrq[1,]])/(bm$t[lrq[2,]]-bm$t[lrq[1,]])-0.5))
      l.idx <- lrq[1,i]
      r.idx <- lrq[2,i]
      q <- lrq[3,i]

      # print(glue("intersection sim at {q} in interval ({bm$t[l.idx]},{bm$t[r.idx]})"))

      sim.condintersection_(bm, l.idx, q, r.idx)
      if(prefer == "bessel") {
        intersection.to.bessel_(bm, c(match(q, bm$layers$t.u), match(q, bm$layers$t.l)))
      }
      if(!isFALSE(refine)) {
        refine_(bm, match(q, bm$layers$t.u), mult)
        refine_(bm, match(q, bm$layers$t.l), mult)
      }

      t.intersection <- t.intersection[-i]
    }
    while(length(t.bessel) > 0) {
      # Do simulation conditional on bessel layers
      lrq <- sapply(t.bessel,
                    function(t) {
                      c(l = tail(which(bm$t<t), 1), r = head(which(bm$t>t), 1), q = t)
                    })
      i <- which.min(abs(abs(lrq[3,]-bm$t[lrq[1,]])/(bm$t[lrq[2,]]-bm$t[lrq[1,]])-0.5))
      l.idx <- lrq[1,i]
      r.idx <- lrq[2,i]
      q <- lrq[3,i]

      # print(glue("bessel sim at {q} in interval ({bm$t[l.idx]},{bm$t[r.idx]})"))

      sim.condbessel_(bm, l.idx, q, r.idx)
      if(prefer == "intersection") {
        bessel.to.intersection_(bm, c(match(q, bm$layers$t.u), match(q, bm$layers$t.l)))
      }
      if(!isFALSE(refine)) {
        refine_(bm, match(q, bm$layers$t.u), mult)
        refine_(bm, match(q, bm$layers$t.l), mult)
      }

      t.bessel <- t.bessel[-i]
    }
    if(length(t.bblocal) > 0) {
      # Do simulation conditional on bb-localised layers
      t.bblocal <- sapply(t.bblocal,
                          function(t) {
                            c(l = tail(which(bm$t<t), 1), r = head(which(bm$t>t), 1), q = t)
                          })
      for(l in unique(t.bblocal[1,])) {
        i <- which(t.bblocal[1,]==l)
        r <- t.bblocal[2,i[1]]
        # Following not needed since internally it dispatches based on hard/soft to aux path
        # if(any(!bm$layers[match(bm$t[r], bm$layers$t.u),c("Lu.hard","Ud.hard")])) {
        #   soft <- TRUE
        # } else {
        #   soft <- FALSE
        # }
        for(qq in rev(t.bblocal[3,i])) {
          sim.condbblocal_(bm, l, qq, r)
          if(prefer == "bessel") {
            intersection.to.bessel_(bm, match(qq, bm$layers$t.u))
          }
          l <- l+1
          r <- r+1
          if(!isFALSE(refine)) {
            refine_(bm, match(qq, bm$layers$t.u), mult)
            refine_(bm, match(qq, bm$layers$t.l), mult)
          }
        }
      }
    }
    while(length(t.bbintersection) > 0) {
      # Do simulation conditional on BB intersection layers
      lrq <- sapply(t.bbintersection,
                    function(t) {
                      c(l = tail(which(bm$t<t), 1), r = head(which(bm$t>t), 1), q = t)
                    })
      i <- which.min(abs(abs(lrq[3,]-bm$t[lrq[1,]])/(bm$t[lrq[2,]]-bm$t[lrq[1,]])-0.5))
      l.idx <- lrq[1,i]
      r.idx <- lrq[2,i]
      q <- lrq[3,i]

      # print(glue("bb-intersection sim at {q} in interval ({bm$t[l.idx]},{bm$t[r.idx]})"))

      sim.condbbintersection_(bm, l.idx, q, r.idx)
      if(prefer == "bessel") {
        intersection.to.bessel_(bm, c(match(q, bm$layers$t.u), match(q, bm$layers$t.l)))
      }
      if(!isFALSE(refine)) {
        refine_(bm, match(q, bm$layers$t.u), mult)
        refine_(bm, match(q, bm$layers$t.l), mult)
      }

      t.bbintersection <- t.bbintersection[-i]
    }
    while(length(t.bbbessel) > 0) {
      # Do simulation conditional on BB bessel layers
      lrq <- sapply(t.bbbessel,
                    function(t) {
                      c(l = tail(which(bm$t<t), 1), r = head(which(bm$t>t), 1), q = t)
                    })
      i <- which.min(abs(abs(lrq[3,]-bm$t[lrq[1,]])/(bm$t[lrq[2,]]-bm$t[lrq[1,]])-0.5))
      l.idx <- lrq[1,i]
      r.idx <- lrq[2,i]
      q <- lrq[3,i]

      # print(glue("bb-bessel sim at {q} in interval ({bm$t[l.idx]},{bm$t[r.idx]})"))

      sim.condbbbessel_(bm, l.idx, q, r.idx)
      if(prefer == "intersection") {
        bessel.to.intersection_(bm, c(match(q, bm$layers$t.u), match(q, bm$layers$t.l)))
      }
      if(!isFALSE(refine)) {
        refine_(bm, match(q, bm$layers$t.u), mult)
        refine_(bm, match(q, bm$layers$t.l), mult)
      }

      t.bbbessel <- t.bbbessel[-i]
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
