#' Copy / Concatenate Brownian motion objects
#'
#' @param ... multiple Brownian motion objects which will be concatenated together in the order provided into a new BrownianMotion object.
#'   Note the object is updated in place
#' @param t0 scalar representing the time the new Brownian motion should be initialised.  Defaults to value of initial time of the first Brownian motion object to be concatenated.
#'
#' @export
concat.bm <- function(..., t0 = NULL) {
  bms <- list(...)

  if(length(bms) == 1 && is.list(bms[[1]]) && length(bms[[1]]) > 1) {
    return(do.call(concat.bm, c(bms[[1]], t0 = t0)))
  } else {
    UseMethod("concat.bm")
  }
}

#' @export
concat.bm.BrownianMotionNd <- function(..., t0 = NULL) {
  bms <- list(...)

  if(length(bms) == 1 && is.list(bms[[1]]) && length(bms[[1]]) > 1) {
    return(do.call(concat.bm, c(bms[[1]], t0 = t0)))
  }

  if(length(bms) > 1) {
    if(!all(sapply(bms, function(bm) { "BrownianMotionNd" %in% class(bm) }))) {
      stop("All arguments must be BrownianMotionNd objects")
    }
  } else {
    if(!("BrownianMotionNd" %in% class(bms[[1]]))) {
      stop("Argument must be BrownianMotionNd objects")
    }
  }

  res <- new.env(parent = emptyenv())
  class(res) <- "BrownianMotionNd"
  res$dim <- bms[[1]]$dim
  res$chol <- bms[[1]]$chol

  if(!all(sapply(bms, function(bm) { bm$dim == res$dim }))) {
    stop("All arguments must be the same dimension")
  }
  if(!all(sapply(bms, function(bm) { identical(bm$chol, res$chol) }))) {
    stop("All arguments must have the same covariance")
  }

  res$Z.bm <- list()
  for(d in 1:res$dim) {
    if(length(bms) == 1) {
      res$Z.bm[[d]] <- concat.bm(bms[[1]]$Z.bm[[d]], t0 = t0)
    } else {
      res$Z.bm[[d]] <- concat.bm(lapply(bms, function(bm) { bm$Z.bm[[d]] }), t0 = t0)
    }
  }

  res
}

#' @export
concat.bm.BrownianMotion <- function(..., t0 = NULL) {
  bms <- list(...)

  if(length(bms) == 1 && is.list(bms[[1]]) && length(bms[[1]]) > 1) {
    return(do.call(concat.bm, c(bms[[1]], t0 = t0)))
  }

  if(length(bms) > 1) {
    if(!all(sapply(bms, function(bm) { "BrownianMotion" %in% class(bm) }))) {
      stop("All arguments must be BrownianMotion objects")
    }
    if(!all(sapply(bms, function(bm) { bm$refine }) == bms[[1]]$refine)) {
      warning("Not all BrownianMotion objects have the same refine parameter ... using setting of first supplied object")
    }
    if(!all(sapply(bms, function(bm) { bm$mult }) == bms[[1]]$mult)) {
      warning("Not all BrownianMotion objects have the same mult parameter ... using setting of first supplied object")
    }
    if(!all(sapply(bms, function(bm) { bm$prefer }) == bms[[1]]$prefer)) {
      warning("Not all BrownianMotion objects have the same prefer parameter ... using setting of first supplied object")
    }
  } else {
    if(!("BrownianMotion" %in% class(bms[[1]]))) {
      stop("Argument must be BrownianMotion objects")
    }
  }

  # Copy first bm object into res
  res <- new.env(parent = emptyenv())
  class(res) <- "BrownianMotion"
  for(nm in ls(bms[[1]], all.names=TRUE)) assign(nm, get(nm, bms[[1]]), res)
  res$bb.local <- new.env(parent = emptyenv())
  class(res$bb.local) <- "BrownianMotion"
  for(nm in ls(bms[[1]]$bb.local, all.names=TRUE)) assign(nm, get(nm, bms[[1]]$bb.local), res$bb.local)

  if(!is.null(t0)) {
    del.t <- t0 - min(res$t)
    res$t <- res$t + del.t
    res$layers[,2:3] <- res$layers[,2:3] + del.t
    res$user.layers[,1:2] <- res$user.layers[,1:2] + del.t
    res$bb.local$t <- res$bb.local$t + del.t
    res$bb.local$layers[,2:3] <- res$bb.local$layers[,2:3] + del.t
    res$bb.local$user.layers[,1:2] <- res$bb.local$user.layers[,1:2] + del.t
    for(lab in names(res$labels)) {
      res$labels[[lab]] <- res$labels[[lab]] + del.t
    }
  }
  if(length(bms) > 1) {
    for(i in 2:length(bms)) {
      del.t <- max(res$t) - min(bms[[i]]$t)

      res$t <- c(res$t, bms[[i]]$t[-1] + del.t)
      res$W_t <- c(head(res$W_t, -1), bms[[i]]$W_t)
      res$W_tm <- c(res$W_tm, bms[[i]]$W_tm[-1])

      res$bb.local$t <- c(res$bb.local$t, bms[[i]]$bb.local$t[-1] + del.t)
      res$bb.local$W_t <- c(head(res$bb.local$W_t, -1), bms[[i]]$bb.local$W_t)
      res$bb.local$W_tm <- c(res$bb.local$W_tm, bms[[i]]$bb.local$W_tm[-1])

      tmp <- bms[[i]]$layers
      tmp[,2:3] <- tmp[,2:3] + del.t
      res$layers <- rbind(res$layers, tmp)
      tmp <- bms[[i]]$user.layers
      tmp[,1:2] <- tmp[,1:2] + del.t
      res$user.layers <- rbind(res$user.layers, tmp)

      tmp <- bms[[i]]$bb.local$layers
      tmp[,2:3] <- tmp[,2:3] + del.t
      res$bb.local$layers <- rbind(res$bb.local$layers, tmp)
      tmp <- bms[[i]]$bb.local$user.layers
      tmp[,1:2] <- tmp[,1:2] + del.t
      res$bb.local$user.layers <- rbind(res$bb.local$user.layers, tmp)

      labs <- unique(c(names(res$labels), names(bms[[i]]$labels)))
      for(lab in labs) {
        res$labels[[lab]] <- c(res$labels[[lab]], bms[[i]]$labels[[lab]] + del.t)
      }
    }
  }

  res$labels[["start"]] <- min(res$t)
  res$labels[["end"]] <- max(res$t)
  res$labels[["seg.start"]] <- c(res$labels[["start"]], res$t[which(res$W_t != res$W_tm)])
  res$labels[["seg.end"]] <- c(res$t[which(res$W_t != res$W_tm)], res$labels[["end"]])

  res
}
