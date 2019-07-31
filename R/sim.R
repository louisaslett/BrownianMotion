#' Simulate Brownian motion at a target time
#'
#' Simulates Brownian motion at chosen times, automatically applying the correct
#' algorithm to condition on the current known (ie previously simulated) states
#' of the trajectory.
#'
#' @param bm a Brownian motion object from which simulation should continue.
#'   Note the object is updated in place.
#' @param t a vector of times to simulate at.
#'
#' @return the Brownian motion object which was passed in argument \code{bm} is
#'   updated in place and returned, enabling chaining of commands with
#'   dplyr (and other) style pipes.
#'
#' @export
sim <- function(bm, t) {
  if(!("BrownianMotion" %in% class(bm))) {
    stop("bm argument must be a BrownianMotion object.")
  }
  if(!is.realvector(t)) {
    stop("t must be a vector of times.")
  }
  if(is.unsorted(t)) {
    t <- sort(t)
  }

  # Eliminate times we know
  t <- setdiff(t, bm$t)

  # Find times which are standard forward simulation at the end
  t.bmfwd <- t[which(t>tail(bm$t, 1))]
  # Find Brownian bridge candidate times
  # Reverse sort them so that index in bm$t is correct over many iterations
  t.bmbb <- sort(t[which(t<tail(bm$t, 1))], decreasing = TRUE)

  # Do forward sims
  if(length(t.bmfwd) > 0) {
    sim.uncond_(bm, t.bmfwd)
  }
  # Do BB sims
  if(length(t.bmbb) > 0) {
    # Check if any are in bounded regions
    if(nrow(bm$bounds) > 0) {
      t.bdd <- t.bmbb[rowSums(matrix(apply(bm$bounds, 1, function(y) { t.bmbb > y[1] & t.bmbb < y[2] }), ncol = 2)) > 0]
      warning(glue("Currently unable to simulate constrained Brownian motion at timepoint {t.bdd}.\n\n"))
      # Keep rest as standard Brownian bridges
      t.bmbb <- setdiff(t.bmbb, t.bdd)
    }
    if(nrow(bm$bessel.layers) > 0) {
      t.lyr <- t.bmbb[rowSums(matrix(apply(bm$bessel.layers, 1, function(y) { t.bmbb > y[1] & t.bmbb < y[2] }), ncol = 2)) > 0]
      warning(glue("Currently unable to simulate Brownian motion within Bessel layers at timepoint {t.lyr}.\n\n"))
      # Keep rest as standard Brownian bridges
      t.bmbb <- setdiff(t.bmbb, t.lyr)
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
