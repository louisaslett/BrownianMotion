#' Eliminate parts of path, layer or both from the skeleton
#'
#' @param bm a Brownian motion object from which simulation should continue.
#'   Note the object is updated in place
#' @param l closed left end of interval to delete.  By default, `-Inf` which
#'     results in left truncation of the path up to `r`
#' @param r open right end of interval to delete.  By default, `Inf` which
#'     results in right truncation of the path from `l`
#' @param type one of `"all"` or `"layer"` to specify whether to delete both path
#'     observations and layers, or just layers in the open interval `(l,r)`
#'
#'
#' @export
delete.skeleton <- function(bm, l = -Inf, r = Inf, type = "all") {
  UseMethod("delete.skeleton")
}

#' @export
delete.skeleton.BrownianMotionNd <- function(bm, ...) {
  if(!("BrownianMotionNd" %in% class(bm))) {
    stop("bm argument must be a BrownianMotionNd object.")
  }

  for(d in 1:bm$dim) {
    delete.skeleton(bm$Z.bm[[d]], ...)
  }

  invisible(bm)
}

#' @export
delete.skeleton.BrownianMotion <- function(bm, l = -Inf, r = Inf, type = "all") {
  if(!("BrownianMotion" %in% class(bm))) {
    stop("bm argument must be a BrownianMotion object.")
  }
  if(length(type) != 1 || !(type %in% c("all", "layer"))) {
    stop('type must be either "all" or "layer"')
  }
  if(!is.realscalar(l)) {
    stop("l must be a scalar.")
  }
  if(!is.realscalar(r)) {
    stop("r must be a scalar.")
  }
  if(l >= r) {
    stop("l must be less than r.")
  }
  if(l == -Inf && r == Inf) {
    stop("at least one of l or r must be specified")
  }
  if(l!=-Inf && !(l %in% bm$t)) {
    stop("l must be a sampled path point")
  }
  if(r!=Inf && !(r %in% bm$t)) {
    stop("r must be a sampled path point")
  }

  ### User layers (complex case) first ###
  # i) layers fully contained are just removed
  usr_idxs <- which(bm$user.layers$t.l>=l & bm$user.layers$t.u<=r)
  if(length(usr_idxs) > 0) {
    bm$user.layers <- bm$user.layers[-usr_idxs,]
  }
  # ii) layers straddling the left and right
  usr_idxs <- which(bm$user.layers$t.l<l & bm$user.layers$t.u>r)
  if(length(usr_idxs) > 1) stop("impossible case ii for user layer matching")
  if(length(usr_idxs) == 1) {
    new.lyr <- bm$user.layers[usr_idxs,,drop=FALSE]
    new.lyr$t.l <- r
    bm$user.layers[usr_idxs,]$t.u <- l
    bm$user.layers <- rbind(bm$user.layers, new.lyr)
  }
  # iii) layers straddling the left only
  usr_idxs <- which(bm$user.layers$t.l<l & bm$user.layers$t.u>l & bm$user.layers$t.u<=r)
  if(length(usr_idxs) > 1) stop("impossible case iii for user layer matching")
  if(length(usr_idxs) == 1) {
    bm$user.layers[usr_idxs,]$t.u <- l
  }
  # iv) layers straddling the right only
  usr_idxs <- which(bm$user.layers$t.l<r & bm$user.layers$t.u>r & bm$user.layers$t.l>=l)
  if(length(usr_idxs) > 1) stop("impossible case iv for user layer matching")
  if(length(usr_idxs) == 1) {
    bm$user.layers[usr_idxs,]$t.l <- r
  }

  ### Standard layers ###
  lyr_idxs <- which(bm$layers$t.l>=l & bm$layers$t.u<=r)
  if(length(lyr_idxs) > 0) {
    bm$layers <- bm$layers[-lyr_idxs,]
  }

  ### Aux layers ###
  aux_lyr_idxs <- which(bm$bb.local$layers$t.l>=l & bm$bb.local$layers$t.u<=r)
  if(length(aux_lyr_idxs) > 0) {
    bm$bb.local$layers <- bm$bb.local$layers[-aux_lyr_idxs,]
  }

  ### Path observations ###
  if(type == "all") {
    obs_idxs <- which(bm$t>l & bm$t<r)
    if(length(obs_idxs)>0) {
      bm$t <- bm$t[-obs_idxs]
      bm$W_t <- bm$W_t[-obs_idxs]
      bm$W_tm <- bm$W_tm[-obs_idxs]
    }
    aux_obs_idxs <- which(bm$bb.local$t>l & bm$bb.local$t<r)
    if(length(aux_obs_idxs)>0) {
      bm$bb.local$t <- bm$bb.local$t[-aux_obs_idxs]
      bm$bb.local$W_t <- bm$bb.local$W_t[-aux_obs_idxs]
      bm$bb.local$W_tm <- bm$bb.local$W_tm[-aux_obs_idxs]
    }
    for(lbl in names(bm$labels)) {
      lbl_idxs <- which(bm$labels[[lbl]]>l & bm$labels[[lbl]]<r)
      if(length(lbl_idxs)>0) {
        bm$labels[[lbl]] <- bm$labels[[lbl]][-lbl_idxs]
        if(length(bm$labels[[lbl]])==0)
          bm$labels[[lbl]] <- NULL
      }
    }
    ### Reconcile W_t and W_tm changes at the deletion times ###
    if(l!=-Inf) {
      l_idx <- which(bm$t == l)
      bm$W_t[l_idx] <- bm$W_tm[l_idx]
    }
    if(r!=Inf) {
      r_idx <- which(bm$t == r)
      bm$W_tm[r_idx] <- bm$W_t[r_idx]
    }
    ### Update labels ###
    bm$labels[["start"]] <- min(bm$t)
    bm$labels[["end"]] <- max(bm$t)
    bm$labels[["seg.start"]] <- c(bm$labels[["start"]], bm$t[which(bm$W_t != bm$W_tm)])
    bm$labels[["seg.end"]] <- c(bm$t[which(bm$W_t != bm$W_tm)], bm$labels[["end"]])
  }

  invisible(bm)
}
