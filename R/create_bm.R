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

create.bm_ <- function(W_0 = 0) {
  bm <- new.env(parent = emptyenv())
  bm$t <- 0
  bm$W_t <- W_0
  bm$bounds <- tibble(t.l = numeric(),
                      t.u = numeric(),
                      l = numeric(),
                      u = numeric(),
                      L = numeric(),
                      U = numeric())
  bm$bessel.layers <- tibble(t.l = numeric(),
                             t.u = numeric(),
                             l = numeric(),
                             u = numeric(),
                             L = numeric(),
                             U = numeric())
  class(bm) <- "BrownianMotion"
  bm
}
