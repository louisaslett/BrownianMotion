#' @export
localise.layers <- function(bm, s, t, refine = bm$refine, mult = bm$mult, prefer = bm$prefer) {
  ## NOTE TO LOUIS: This shares a lot of setup with Bessel Layers -- pull out into utility code (except the find Bessel layers part)
  if(!("BrownianMotion" %in% class(bm))) {
    stop("bm argument must be a BrownianMotion object.")
  }
  if(!is.realscalar(s) || !is.realscalar(t) || !is.realscalar(mult)) {
    stop("s, t and mult must be real scalars.")
  }
  if(s >= t) {
    stop("s must be less than t.")
  }
  if(mult <= 0) {
    stop("mult must be strictly positive.")
  }

  # Are these endpoints part of the skeleton?
  # If not, and the point is beyond the end, simulate it.  Otherwise, error
  if(!(s %in% bm$t)) {
    sim(bm, s)
  }
  if(!(t %in% bm$t)) {
    # In the localised case we only sim t directly if we are not past the
    # end of the skeleton, otherwise we do repeated first passages
    if(t < max(bm$t)) {
      sim(bm, t, refine, mult, prefer)
    } else {
      s.tmp <- max(bm$t)
      theta <- mult*sqrt(t-s.tmp)
      while(max(bm$t) < t) {
        first.passage(bm, delta = theta/2)
      }
      # Remove green user layers as these first passages were not a result of a user first passage request
      bm$user.layers <- bm$user.layers[!(bm$user.layers$t.l >= s.tmp & bm$user.layers$t.u <= max(bm$t)),]
      if(!(t %in% bm$t)) {
        # t is not the end point, so conditionally simulate at t and then chop off the end
        sim(bm, t, refine, mult, prefer)
        bm$t <- head(bm$t, -1)
        bm$W_t <- head(bm$W_t, -1)
        # Last layer is now starting at t, so remove it
        rm.lyr_(bm, get.lyr_(bm, t)$idx)
      }
    }
  }

  # Find all intervals between s and t not currently inside any other layer or inside a Bessel layer
  s.i <- which(s == bm$t)
  t.i <- which(t == bm$t)
  all.pairs.noL <- matrix(c(bm$t[s.i:(t.i-1)], bm$t[(s.i+1):t.i]), byrow = FALSE, ncol = 2)
  # No layer
  mid <- (all.pairs.noL[,2]+all.pairs.noL[,1])/2
  incl <- rep(TRUE, nrow(all.pairs.noL))
  for(i in 1:length(mid)) {
    if(any(mid[i] >= bm$layers$t.l & mid[i] <= bm$layers$t.u) ||
       any(mid[i] >= bm$user.layers$t.l & mid[i] <= bm$user.layers$t.u)) {
      incl[i] <- FALSE
    }
  }
  all.pairs.noL <- all.pairs.noL[incl,,drop = FALSE]
  # cat("No layer:"); print(all.pairs.noL)

  if(nrow(all.pairs.noL) == 0) {
    stop("No intervals in the skeleton between s and t found that have no layer specification.")
  }

  if(nrow(all.pairs.noL) > 0) {
    for(i in 1:nrow(all.pairs.noL)) {
      bb.localise_(bm, all.pairs.noL[i,1], all.pairs.noL[i,2], refine, mult, prefer)
    }
  }

  bm$layers <- bm$layers[order(bm$layers$t.l),]
  invisible(bm)
}
