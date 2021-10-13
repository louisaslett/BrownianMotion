#' Create a Brownian motion object
#'
#' Creates an R environment to contain the trajectories and layer information of
#' a Brownian motion.
#'
#' @param t vector of fixed times of the Brownian motion.  Defaults to 0.
#' @param W_t vector of corresponding locations of the Brownian motion at the times t.  Defaults to 0.
#' @param refine refine is
#' @param mult mult is
#' @param prefer prefer is
#'
#' @export
create.bm <- function(t = 0, W_t = 0, refine = TRUE, mult = 1, prefer = "bessel") {
  if(!is.numeric(t)) {
    stop("Argument t must be a numeric vector")
  }
  if(!is.numeric(W_t)) {
    stop("Argument W_t must be a numeric vector")
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

  create.bm_(t, W_t, refine, mult, prefer)
}

create.bm_ <- function(t, W_t, refine, mult, prefer, nested = FALSE) {
  bm <- new.env(parent = emptyenv())
  bm$W_t <- W_t[order(t)]
  bm$W_tm <- W_t[order(t)]
  bm$t <- sort(t)
  bm$labels <- list()
  bm$labels[["start"]] <- min(t)
  bm$labels[["end"]] <- max(t)
  bm$labels[["seg.start"]] <- min(t)
  bm$labels[["seg.end"]] <- max(t)
  bm$labels[["forced"]] <- sort(t)
  bm$refine <- refine
  bm$mult <- mult
  bm$prefer <- prefer
  bm$layers <- tibble(type = factor(levels = c("localised", "localised-bb", "bessel", "bessel-bb", "intersection", "intersection-bb")),
                      t.l = numeric(),
                      t.u = numeric(),
                      Ld = numeric(),
                      Uu = numeric(),
                      Lu = numeric(),
                      Ud = numeric(),
                      Lu.hard = logical(),
                      Ud.hard = logical())
  bm$user.layers <- tibble(t.l = numeric(),
                           t.u = numeric(),
                           L = numeric(),
                           U = numeric()) # What was user specified vs what was actually computed to achieve the simulation
  if(!nested) {
    bm$bb.local <- create.bm_(t, W_t, refine, mult, prefer, nested = TRUE) # Auxilliary info for Brownian bridge localised layers
  }
  class(bm) <- "BrownianMotion"
  bm
}
