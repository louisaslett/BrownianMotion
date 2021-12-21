getlayersNd_ <- function(bm, t) {
  if(length(t) != 1) {
    stop("Only one time point supported")
  }

  # Check if there is no layer at the specified time and return early
  if(sum(bm$Z.bm[[1]]$layers$t.l <= t & bm$Z.bm[[1]]$layers$t.u > t) == 0) {
    return(NULL) # or should this be tibble with NULL entries?
  }

  # Pull out relevant row of layers for the given time
  lyrs <- sapply(bm$Z.bm, function(Z.bm) {
    Z.bm$layers[Z.bm$layers$t.l <= t & Z.bm$layers$t.u > t,c("t.l", "t.u", "Ld", "Uu")]
  })

  # Check time intervals in all dimensions are in agreement
  if(!all(c(lyrs[1,], recursive = TRUE)==lyrs[1,1]) || !all(c(lyrs[2,], recursive = TRUE)==lyrs[2,1])) {
    stop("Inconsistent layer times in differing dimensions")
  }

  # Make a single row tibble with the transformed layer
  tibble(
    t.l = lyrs[1,1][[1]],
    t.u = lyrs[2,1][[1]],
    L = matrix(c(lyrs[3,], recursive = TRUE), nrow = 1),
    U = matrix(c(lyrs[4,], recursive = TRUE), nrow = 1)
  )
}

transformlayersNd_ <- function(bm, t) {
  lyr <- getlayersNd_(bm, t)
  if(is.null(lyr))
    return(tibble(
      t.l = numeric(),
      t.u = numeric(),
      outer.L = matrix(numeric(), nrow = 0, ncol = bm$dim),
      outer.U = matrix(numeric(), nrow = 0, ncol = bm$dim),
      inner.cube = list()
    ))

  LU <- rbind(lyr$L,lyr$U)
  Z.vertices <- as.matrix(do.call("expand.grid", lapply(1:bm$dim, function(d) { LU[,d] })))
  inner.cube <- tcrossprod(Z.vertices, bm$chol)
  outer.L <- apply(inner.cube, 2, min)
  outer.U <- apply(inner.cube, 2, max)
  # W.outer <- expand.grid(lapply(1:bm$dim, function(d) { c(outer.L[d], outer.U[d]) }))

  tibble(
    t.l = lyr[1,1][[1]],
    t.u = lyr[1,2][[1]],
    outer.L = matrix(outer.L, nrow = 1),
    outer.U = matrix(outer.U, nrow = 1),
    inner.cube = list(inner.cube)
  )
}

getZtNd_ <- function(bm, i) {
  do.call("cbind", lapply(bm$Z.bm, function(Z.bm) { unname(Z.bm$W_t[i]) }))
}

getZtmNd_ <- function(bm, i) {
  do.call("cbind", lapply(bm$Z.bm, function(Z.bm) { unname(Z.bm$W_tm[i]) }))
}

getWtNd_ <- function(bm, i) {
  tcrossprod(getZtNd_(bm, i), bm$chol)
}

getWtmNd_ <- function(bm, i) {
  tcrossprod(getZtmNd_(bm, i), bm$chol)
}

