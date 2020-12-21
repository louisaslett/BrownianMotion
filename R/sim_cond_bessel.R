#' Simulate Brownian motion conditional on bessel layer
#'
#' This currently converts the Bessel layer to an intersection layer and then
#' uses the conditional intersection layer sampling.
#'
#' @export
sim.condbessel <- function(bm, s, t, q = NULL, q.grid = NULL) {
  # Arg types & combos
  if(!("BrownianMotion" %in% class(bm))) {
    stop("bm argument must be a BrownianMotion object.")
  }
  if(as.integer(is.null(q)) + as.integer(is.null(q.grid)) != 1) {
    stop("Exactly one of q or q.grid must be specified")
  }
  if(!is.null(q.grid) && !is.intscalar(q.grid)) {
    stop("q.grid must be either NULL or an integer number of grid points.")
  }
  if(!is.null(q) && !is.realvector(q)) {
    stop("q must be either NULL or a vector of times.")
  }
  if(!is.realscalar(t)) {
    stop("t must be a scalar.")
  }
  if(!is.realscalar(s)) {
    stop("s must be a scalar.")
  }

  # Checks
  if(!(s %in% bm$t)) {
    stop("s not found in bm path.")
  }
  if(!(t %in% bm$t)) {
    stop("t not found in bm path.")
  }
  bessel <- bm$layers[bm$layers$type == "bessel",]
  if(!(t %in% bessel$t.u)) {
    stop("t is not the time of a realised maximum/minimum in a bessel layer.")
  }
  # Do we then need to check that s is the start of this same bessel layer?
  # Actually this is taken care of by the following which allows no intermediate obs
  s_idx <- match(s, bm$t)
  t_idx <- match(t, bm$t)
  if(t_idx != s_idx+1) {
    stop("Cannot have other brownian motion observations between times s and t")
  }

  # Make grid if it wasn't supplied
  if(is.null(q)) {
    q <- head(tail(seq(s, t, length.out = q.grid+2), -1), -1)
  }
  if(any(q < s) || any(q > t)) {
    stop("All q must lie between s and t")
  }
  if(is.unsorted(q)) {
    q <- sort(q)
  }
  # Eliminate times we know
  q <- setdiff(q, bm$t)

  for(qq in q) {
    bm.res <- sim.condbessel_(bm, s_idx, qq, t_idx)
    s_idx <- s_idx+1
    t_idx <- t_idx+1
  }

  bm.res$layers <- bm.res$layers[order(bm.res$layers$t.l),]
  invisible(bm.res)
}

sim.condbessel_ <- function(bm, s_idx, q, t_idx) {
  s <- bm$t[s_idx]
  t <- bm$t[t_idx]

  # Update layer info
  # We split the layer in two, adding left and right layers either side of the
  # newly simulated time, then remove the old layer
  cur.layer <- which(bm$layers$t.l == s & bm$layers$t.u == t)

  if(bm$layers$type[cur.layer] == "bessel") {
    bm <- intersection.layers.IL_(bm, s, t, 1)
    bm <- sim.condintersection_(bm, s_idx, q, t_idx)
  } else if(bm$layers$type[cur.layer] == "intersection") {
    bm <- sim.condintersection_(bm, s_idx, q, t_idx)
  } else {
    stop("Error: must have bessel or intersection layer at this point.")
  }

  bm
}
