#' Manually add labels after simulation
#'
#' Description of adding labels
#'
#' @param bm a Brownian motion object from which simulation should continue.
#'   Note the object is updated in place
#' @param t vector of fixed times of the Brownian motion.
#' @param label vector of length 1 or the length of t indicating the introduction of user specified labels (if any). By default, uses any names in the vector of times.
#'
#' @export
add.labels <- function(bm, t, label = names(t)) {
  UseMethod("add.labels")
}

#' @export
add.labels.BrownianMotionNd <- function(bm, t, label = names(t)) {
  if(!("BrownianMotionNd" %in% class(bm))) {
    stop("bm argument must be a BrownianMotionNd object.")
  }

  for(d in 1:bm$dim) {
    add.labels(bm$Z.bm[[d]], t, label)
  }

  invisible(bm)
}

#' @export
add.labels.BrownianMotion <- function(bm, t, label = names(t)) {
  if(!("BrownianMotion" %in% class(bm))) {
    stop("bm argument must be a BrownianMotion object.")
  }

  # Check times
  if(any(!(t %in% bm$t))) {
    stop("Can only add labels to existing observed time points.")
  }

  assert.bmlabel(label, t)

  add.labels_(bm, label, t)

  invisible(bm)
}



#' Manually delete labels after simulation
#'
#' Description of deleting labels
#'
#' @param bm a Brownian motion object from which simulation should continue.
#'   Note the object is updated in place
#' @param t vector of fixed times of the Brownian motion.
#' @param label vector of length 1 or the length of t indicating the introduction of user specified labels (if any). By default, uses any names in the vector of times.
#'
#' @export
delete.labels <- function(bm, t = NA, label = NA) {
  UseMethod("delete.labels")
}

#' @export
delete.labels.BrownianMotionNd <- function(bm, t = NA, label = NA) {
  if(!("BrownianMotionNd" %in% class(bm))) {
    stop("bm argument must be a BrownianMotionNd object.")
  }

  for(d in 1:bm$dim) {
    delete.labels(bm$Z.bm[[d]], t, label)
  }

  invisible(bm)
}

#' @export
delete.labels.BrownianMotion <- function(bm, t = NA, label = NA) {
  if(!("BrownianMotion" %in% class(bm))) {
    stop("bm argument must be a BrownianMotion object.")
  }

  if(any(label %in% c("start", "end", "user", "seg.start", "seg.end", "fpt", "internal", "forced"))) {
    stop("one of the specified labels clashes with an internal label naming scheme: label not deleted.")
  }
  if(length(label) == 1 && is.na(label)) {
    if(length(t) == 1 && is.na(t)) {
      stop("must specify at least one of t and label arguments")
    }
    label <- setdiff(names(bm$labels), c("start", "end", "user", "seg.start", "seg.end", "fpt", "internal", "forced"))
  }

  if(length(t) == 1 && is.na(t)) {
    # Delete all of a label ignoring time

    for(l in label) {
      bm$labels[[l]] <- NULL
    }
  } else {
    # Time matters
    # Check times
    if(any(!(t %in% bm$t))) {
      stop("Can only delete labels for existing observed time points.")
    }

    for(tt in t) {
      for(lbl in label) {
        lbl_idxs <- which(bm$labels[[lbl]]==tt)
        if(length(lbl_idxs)>0) {
          bm$labels[[lbl]] <- bm$labels[[lbl]][-lbl_idxs]
          if(length(bm$labels[[lbl]])==0)
            bm$labels[[lbl]] <- NULL
        }
      }
    }
  }

  invisible(bm)
}
