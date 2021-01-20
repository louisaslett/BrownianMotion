#' Create a Brownian motion object
#'
#' Creates an R environment to contain the trajectories and layer information of
#' a Brownian motion.
#'
#' @param W_0 the starting location of the Brownian motion at time 0.
#'
#' @export
create.bm <- function(W_0 = 0) {
  if(!is.realscalar(W_0)) {
    stop("Argument W_0 must be a scalar")
  }

  create.bm_(W_0)
}

create.bm_ <- function(W_0 = 0, nested = FALSE) {
  bm <- new.env(parent = emptyenv())
  bm$t <- 0
  bm$W_t <- W_0
  bm$layers <- tibble(type = factor(levels = c("localised", "localised-bb", "bessel", "intersection")),
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
    bm$bb.local <- create.bm_(nested = TRUE) # Auxilliary info for Brownian bridge localised layers
  }
  class(bm) <- "BrownianMotion"
  bm
}
