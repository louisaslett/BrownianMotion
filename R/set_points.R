#' Manually specify the value of a point in the skeleton
#'
#' @param t time
#' @param W_t value
#' @param add.segment TRUE/FALSE.  If true, then this creates a jump
#'     discontinuity at the specified time by updating only W_t and not W_tm.
#'     If false, then this updates both W_t and W_tm.
#'     Ignored if t is a new rather than existing time in the skeleton.
#'
#' @export
set.points <- function(bm, t, W_t, add.segment = FALSE, label = names(t)) {
  if(!("BrownianMotion" %in% class(bm))) {
    stop("bm argument must be a BrownianMotion object.")
  }
  if(!is.realvector(t)) {
    stop("t must be a vector of times.")
  }
  if(!is.realvector(W_t)) {
    stop("W_t must be a vector of observations.")
  }
  if(length(t) != length(W_t)) {
    stop("t and W_t must be the same length.")
  }

  assert.bmlabel(label, t)
  if(!is.null(label) && length(label) == 1) {
    label <- rep(label, length(t))
  }

  new.t <- c()
  new.W_t <- c()
  new.W_tm <- c()
  mod.t_idx <- c()
  mod.W_t <- c()
  mod.W_tm <- c()
  for(i in 1:length(t)) {
    if(t[i] > max(bm$t)) { # new point beyond end of path, so put straight in
      new.t <- c(new.t, t[i])
      new.W_t <- c(new.W_t, W_t[i])
      new.W_tm <- c(new.W_tm, W_t[i])
    } else {
      lyr_idx <- which(bm$layers$t.l <= t[i] & bm$layers$t.u >= t[i])
      if(length(lyr_idx) != 0) { # no layer constraint, add in
        stop("Currently points can only be manually specified for times not lying within a layer.")
      }
      if(t[i] %in% bm$t) { # existing point
        t_idx <- which(bm$t == t[i])
        mod.t_idx <- c(mod.t_idx, t_idx)
        mod.W_t <- c(mod.W_t, W_t[i])
        if(!add.segment) {
          mod.W_tm <- c(mod.W_tm, W_t[i])
        } else {
          mod.W_tm <- bm$W_tm[t_idx]
        }
      } else { # new point before end
        new.t <- c(new.t, t[i])
        new.W_t <- c(new.W_t, W_t[i])
        new.W_tm <- c(new.W_tm, W_t[i])

        # This is first pass at trying to allow adding points manually within a layer, but come back to when we have knickers(!)
        # if(W_t > bm$layers$Uu[lyr_idx] || W_t < bm$layers$Ld[lyr_idx]) {
        #   stop(paste0("new observation at time ",t[i],"violates existing layer constraint."))
        # }
        # if(bm$layers$type[lyr_idx] == "intersection") {
        #   intersection.to.bessel_(bm, lyr_idx)
        # }
        # if(bm$layers$type[lyr_idx] == "bessel") {
        #   s_idx <- which(bm$t == bm$layers$t.l[lyr_idx])
        #   t_idx <- which(bm$t == bm$layers$t.u[lyr_idx])
        #   sim.condbessel_(bm, s_idx, t[i], t_idx, label, W_t[i])
        # } else {
        #
        # }
      }
    }
  }
  bm$t <- c(bm$t, new.t)
  bm$W_t <- c(bm$W_t, new.W_t)
  bm$W_tm <- c(bm$W_tm, new.W_tm)
  bm$W_t[mod.t_idx] <- mod.W_t
  bm$W_tm[mod.t_idx] <- mod.W_tm




  invisible(bm)
}
