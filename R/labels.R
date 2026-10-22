#' Manually add labels after simulation
#'
#' Attaches user labels to times already in the skeleton. Labelled times can
#' then be used for indexing, e.g. \code{bm["mylabel", "W_t"]}. A time can
#' carry several labels, and the names used internally by the package (such as
#' \code{"start"}, \code{"end"} and \code{"user"}) cannot be used.
#'
#' @param bm a Brownian motion object.
#'   Note the object is updated in place
#' @param t vector of times in the skeleton to label.
#' @param label vector of length 1 or the length of t indicating the introduction of user specified labels (if any). By default, uses any names in the vector of times.
#'
#' @return The Brownian motion object \code{bm}, updated in place and returned
#'   invisibly so that calls can be chained with pipes.
#'
#' @examples
#' set.seed(1)
#' bm <- create_bm()
#' bm |> sim(1:4)
#'
#' bm |> add_labels(t = c(1, 3), label = "odd")
#' bm |> add_labels(t = c(peak = 2))  # names of t give the labels
#' bm["odd", "t"]
#' bm[, "labels"]
#' @seealso \code{\link{delete_labels}}
#' @export
add_labels <- function(bm, t, label = names(t)) {
  UseMethod("add_labels")
}

#' @export
add_labels.BrownianMotionNd <- function(bm, t, label = names(t)) {
  if(!("BrownianMotionNd" %in% class(bm))) {
    stop("bm argument must be a BrownianMotionNd object.")
  }

  for(d in 1:bm$dim) {
    add_labels(bm$Z.bm[[d]], t, label)
  }

  invisible(bm)
}

#' @export
add_labels.BrownianMotion <- function(bm, t, label = names(t)) {
  if(!("BrownianMotion" %in% class(bm))) {
    stop("bm argument must be a BrownianMotion object.")
  }

  # Check times
  if(any(!(t %in% bm$t))) {
    stop("Can only add labels to existing observed time points.")
  }

  assert.bmlabel(label, t)

  add_labels_(bm, label, t)

  invisible(bm)
}



#' Manually delete labels after simulation
#'
#' Removes user labels. Give \code{label} alone to remove those labels from all
#' times, \code{t} alone to remove all user labels from those times, or both to
#' remove those labels from those times only. Labels used internally by the
#' package cannot be deleted.
#'
#' @param bm a Brownian motion object.
#'   Note the object is updated in place
#' @param t vector of times in the skeleton from which to remove labels.
#' @param label vector of user labels to remove.
#'
#' @return The Brownian motion object \code{bm}, updated in place and returned
#'   invisibly so that calls can be chained with pipes.
#'
#' @examples
#' set.seed(1)
#' bm <- create_bm()
#' bm |> sim(c(a = 1, b = 2, a = 3))
#'
#' bm |> delete_labels(t = 3, label = "a")  # "a" remains at t = 1
#' bm |> delete_labels(label = "b")         # "b" removed everywhere
#' bm[, "labels"]
#' @seealso \code{\link{add_labels}}
#' @export
delete_labels <- function(bm, t = NA, label = NA) {
  UseMethod("delete_labels")
}

#' @export
delete_labels.BrownianMotionNd <- function(bm, t = NA, label = NA) {
  if(!("BrownianMotionNd" %in% class(bm))) {
    stop("bm argument must be a BrownianMotionNd object.")
  }

  for(d in 1:bm$dim) {
    delete_labels(bm$Z.bm[[d]], t, label)
  }

  invisible(bm)
}

#' @export
delete_labels.BrownianMotion <- function(bm, t = NA, label = NA) {
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
