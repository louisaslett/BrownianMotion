#' Construct new layers
#'
#' TODO UPDATE ME
#' Performs unbiased simulation of the smallest specified layers which contain the
#' Brownian motion between two known sample points by retrospective Bernoulli
#' sampling and inversion sampling.
#'
#' @param bm a Brownian motion object from which simulation should continue.
#'   Note the object is updated in place
#' @param s left hand time point
#' @param t right hand time point
#' @param mult the default layer size is \code{sqrt(t-s)}.  You can scale
#'   this by specifying the \code{mult} argument which will result in layer
#'   sizes of \code{mult*sqrt(t-s)}.
#'
#' @return
#' The Brownian motion object which was passed in argument \code{bm} is
#' updated in place and returned, enabling chaining of commands with
#' dplyr (and other) style pipes.
#'
#' @export
layers <- function(bm, s, t, type = bm$prefer, refine = bm$refine, mult = bm$mult, prefer = bm$prefer, label = c(names(s), names(t))) {
  UseMethod("layers")
}

#' @export
layers.BrownianMotionNd <- function(bm, ...) {
  if(!("BrownianMotionNd" %in% class(bm))) {
    stop("bm argument must be a BrownianMotionNd object.")
  }

  old.t <- bm$Z.bm[[1]]$t # Extract a vector of the current times simulated in the BrownianMotionNd object
  new.t <- numeric(0) # Storage vector for new times to be simulated to nsure all dimensions have the same temporal resolution

  for(d in 1:bm$dim) { # Simulate layers for each dimension
    layers(bm$Z.bm[[d]], ...)
  }

  for(d in 1:bm$dim) { # Extract the newly simulated times
    new.t <- unique(c(new.t, setdiff(bm$Z.bm[[d]]$t, old.t)))
  }

  for(d in 1:bm$dim) { # Cross-populate times on each dimension
    sim(bm$Z.bm[[d]], t = new.t)
  }

  invisible(bm)
}

#' @export
layers.BrownianMotion <- function(bm, s, t, type = bm$prefer, refine = bm$refine, mult = bm$mult, prefer = bm$prefer, label = c(names(s), names(t))) {
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
  if(!(type %in% c("bessel", "intersection", "localised"))) {
    stop("layer type must be one of 'bessel', 'intersection', or 'localised'")
  }
  if(!(prefer %in% c("bessel", "intersection"))) {
    stop("layer type must be one of 'bessel' or 'intersection'")
  }
  assert.bmlabel(label, 1:2)
  if(!is.null(label) && length(label) == 1) {
    names(s) <- names(t) <- if(is.na(label)) { NULL } else { label }
  } else if(!is.null(label)) {
    names(s) <- if(is.na(label[1])) { NULL } else { label[1] }
    names(t) <- if(is.na(label[2])) { NULL } else { label[2] }
  }

  # Are these endpoints part of the skeleton?
  # If not, and the point is beyond the end, simulate it.  Otherwise, error
  if(!(s %in% bm$t)) {
    sim(bm, s, refine, mult, prefer)
  }

  if(type != "localised") {
    if(!(t %in% bm$t)) {
      sim(bm, t, refine, mult, prefer)
    }
  } else {
    # In the localised case we only sim t directly if we are not past the
    # end of the skeleton, otherwise we do repeated first passages
    if(!(t %in% bm$t)) {
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
          delete.skeleton(bm, l = t)
        }
      }
    }
  }

  # Extracting intervals over the s&t defined, and determining existing layer type (if any)
  bm$layers <- bm$layers[order(bm$layers$t.l), ] # **** TODO: remove this once layer table is kept ordered
  s.i <- which(s == bm$t) # Index of the left hand most time
  t.i <- which(t == bm$t) # Index of the right hand most time
  all.pairs <- matrix(c(bm$t[s.i:(t.i-1)], bm$t[(s.i+1):t.i]), byrow = FALSE, ncol = 2) # All pairs of consecutive times between s + t, which define intervals
  is.in.layer <- which(all.pairs[,1] %in% bm$layers$t.l) # Extracting intervals which have existing intervals
  all.types <- rep("none", nrow(all.pairs))
  all.types[is.in.layer] <- bm$layers$type[which(bm$layers$t.l >= s & bm$layers$t.u <= t)] # Setting intervals to "none" layer by default and adding in existing layers

  # Eliminate any layer pairs which are not of a type we're going to convert (keeping also none for new layer sims)
  if(type == "localised") {
    keep <- which(all.types == "none")
  } else {
    keep <- which(!(all.types %in% c("localised", "localised-bb", type, paste0(type,"-bb"))))
  }
  all.pairs <- all.pairs[keep,,drop = FALSE]
  all.types <- all.types[keep]

  if(nrow(all.pairs) > 0) {
    for(i in 1:nrow(all.pairs)) {
      if(all.types[i] == "none") {
        if(type == "bessel") {
          bessel.layers_(bm, all.pairs[i,1], all.pairs[i,2], mult)
        } else if(type == "intersection") {
          intersection.layers.BLIL_(bm, all.pairs[i,1], all.pairs[i,2], mult)
        } else if(type == "localised") {
          bb.localise_(bm, all.pairs[i,1], all.pairs[i,2], refine, mult, prefer)
        } else {
          stop("Error new.layers() #1")
        }
      } else if(all.types[i] == "intersection" || all.types[i] == "intersection-bb") {
        intersection.to.bessel_(bm, which(bm$layers$t.l == all.pairs[i,1], bm$layers$t.u == all.pairs[i,2]))
      } else if(all.types[i] == "bessel" || all.types[i] == "bessel-bb") {
        # intersection.layers.IL_(bm, all.pairs[i,1], all.pairs[i,2], mult) # Old way of doing
        bessel.to.intersection_(bm, which(bm$layers$t.l == all.pairs[i,1], bm$layers$t.u == all.pairs[i,2]))
      } else {
        stop("Error new.layers() #2")
      }
    }
  }

  bm$layers <- bm$layers[order(bm$layers$t.l),]
  invisible(bm)
}
