#' Create a Brownian motion object
#'
#' Creates an R environment to contain the trajectories and layer information of
#' a Brownian motion.
#'
#' @param t vector of fixed times of the Brownian motion.  Defaults to value of 0.
#' @param W_t matrix of size the length of t (rows) x dimension (column). If dim = 1, a vector of the length of t is acceptable. If dim > 1 and the length of t = 1, then a vector of length 1 is acceptable (corresponding to the same value at every co-ordinate), or length d is acceptable (corresponding to the d components). Defaults to 0 for every component.
#' @param dim scalar corresponding to the dimension of the Brownian motion.
#' @param cov covariance matrix of the Brownian motion. If scalar corresponds to a multiple of the identity matrix. If matrix should be of size dim x dim.
#' @param refine indicates whether refinement should be by default carried out. If of length 1, then this specifies refinement for every dimension. If of length dim, then this specifies refinement for each dimension (in order).
#' @param mult indicates in the case of refinement the threshold layer size before refinement is carried out. If of length 1, then the default layer size is \code{sqrt(t-s)}, and you can scale this by specifying the \code{mult} argument which will result in layer sizes of \code{mult*sqrt(t-s)}, and this is common for every dimension. If of length dim, then this specifies the \code{mult} for each dimension (in order).
#' @param prefer indicates whether there is a preference for "bessel" or "intersection" layers where possible. If of length 1, then this specifies the preference for every dimension. If of length dim, then this specifies preference for each dimension (in order).
#'
#' @export
create.bm <- function(t = 0, W_t = 0, dim = 1, cov = 1, refine = TRUE, mult = 1, prefer = "bessel") {
  # t checks
  if(!isTRUE(e <- checkmate::check_vector(t, strict = TRUE, any.missing = FALSE, min.len = 1, unique = TRUE, null.ok = FALSE))) {
    stop("Error with argument t: ", e)
  }

  # dim checks
  if(!isTRUE(e <- checkmate::check_int(dim, lower = 1))) {
    stop("Error with argument dim: ", e)
  }
  dim <- as.integer(dim)

  # W_t checks
  if(dim == 1L) {
    if(!isTRUE(e <- checkmate::check_vector(W_t, strict = TRUE, any.missing = FALSE, len = length(t), null.ok = FALSE))) {
      stop("Error with argument W_t: ", e)
    }
    if(!is.matrix(W_t)) {
      W_t <- matrix(W_t, ncol = 1)
    }
  } else if(dim > 1L) {
    if(length(t) == 1) { # allow W_t to be scalar or vector
      if(is.matrix(W_t) && nrow(W_t) == 1) {
        W_t <- W_t[1,]
      }
      if(!isTRUE(e <- checkmate::check_vector(W_t, strict = TRUE, any.missing = FALSE, null.ok = FALSE))) {
        stop("Error with argument W_t: ", e)
      }
      if(length(W_t) == 1) {
        W_t <- rep(W_t, dim)
      } else if(length(W_t) != dim) {
        stop("Error with argument W_t: must be of length 1 or of length dim")
      }
      W_t <- matrix(W_t, nrow = 1)
    } else { # W_t must be matrix
      if(!isTRUE(e <- checkmate::check_matrix(W_t, any.missing = FALSE, nrows = length(t), ncols = dim))) {
        stop("Error with argument W_t: ", e)
      }
    }
  }

  # cov checks
  if(is.matrix(cov)) {
    if(!isTRUE(e <- checkmate::check_matrix(cov, any.missing = FALSE, nrows = dim, ncols = dim))) {
      stop("Error with argument cov: ", e)
    }
    if(!isSymmetric(cov)) {
      stop("Error with argument cov: the covariance matrix is not symmetric")
    }
    e <- eigen(cov, symmetric = TRUE)
    if(any(e$values <= 0)) {
      stop("Error with argument cov: the covariance matrix is not positive definite")
    }
  } else {
    if(!isTRUE(e <- checkmate::check_scalar(cov))) {
      stop("Error with argument cov: must be either a scalar or matrix (dim x dim)")
    }
    if(cov <= 0) {
      stop("Error with argument cov: scalar cov must be strictly positive")
    }
    cov <- diag(cov, nrow = dim, ncol = dim)
  }

  # refine, mult, prefer checks
  if(!isTRUE(e <- checkmate::check_logical(refine, any.missing = FALSE))) {
    stop("Error with argument refine: ", e)
  }
  if(length(refine) != 1 && length(refine) != dim) {
    stop("Error with argument refine: must be length 1 or dim")
  }
  if(length(refine) == 1) {
    refine <- rep(refine, dim)
  }

  if(!isTRUE(e <- checkmate::check_numeric(mult, finite = TRUE, any.missing = FALSE))) {
    stop("Error with argument mult: ", e)
  }
  if(any(mult <= 0)) {
    stop("Error with argument mult: must be strictly positive")
  }
  if(length(mult) != 1 && length(mult) != dim) {
    stop("Error with argument mult: must be length 1 or dim")
  }
  if(length(mult) == 1) {
    mult <- rep(mult, dim)
  }

  if(!isTRUE(e <- checkmate::check_character(prefer, any.missing = FALSE))) {
    stop("Error with argument prefer: ", e)
  }
  prefer <- tolower(prefer)
  if(any(!(prefer %in% c("bessel", "intersection")))) {
    stop("Error with argument prefer: must be one of 'bessel' or 'intersection'")
  }
  if(length(prefer) != 1 && length(prefer) != dim) {
    stop("Error with argument prefer: must be length 1 or dim")
  }
  if(length(prefer) == 1) {
    prefer <- rep(prefer, dim)
  }


  if(!identical(cov, matrix(1, 1, 1))) {
    bm <- new.env(parent = emptyenv())
    bm$cov <- cov
    bm$chol <- chol(bm$cov)
    bm$chol.inv <- solve(chol(bm$cov))
    bm$dim <- dim
    bm$Z.bm <- list()
    Z_t <- tcrossprod(W_t, bm$chol.inv)
    for(d in 1:dim) {
      bm$Z.bm[[d]] <- create.bm_(t, Z_t[,d], refine[d], mult[d], prefer[d])
    }
    # Dummy for accessing untransformed path
    bm$Z <- list(bm)
    class(bm$Z) <- "BrownianMotionNdZ"
    # OLD IDEA: Add storage for outer path and layers containing the post simulation
    # transformed values (decided against compute on demand for this)
    # NEW IDEA: compute on demand first and maybe add caching so don't recompute
    # when requested again
    class(bm) <- "BrownianMotionNd"
  } else {
    bm <- create.bm_(t, W_t[,1], refine[1], mult[1], prefer[1])
  }

  bm
}

create.bm_ <- function(t, W_t, refine, mult, prefer, nested = FALSE) {
  bm <- new.env(parent = emptyenv())
  bm$W_t <- W_t[order(t)]
  bm$W_tm <- W_t[order(t)]
  bm$t <- sort(t)
  bm$labels <- list()
  bm$labels[["start"]] <- min(t)
  bm$labels[["end"]] <- max(t)
  bm$labels[["seg.start"]] <- min(t)
  bm$labels[["seg.end"]] <- max(t)
  bm$labels[["forced"]] <- sort(t)
  bm$refine <- refine
  bm$mult <- mult
  bm$prefer <- prefer
  bm$layers <- tibble(type = factor(levels = c("localised", "localised-bb", "bessel", "bessel-bb", "intersection", "intersection-bb")),
                      t.l = numeric(),
                      t.u = numeric(),
                      Ld = numeric(),
                      Uu = numeric(),
                      Lu = numeric(),
                      Ud = numeric(),
                      Lu.hard = logical(),
                      Ud.hard = logical())
  bm$user.layers <- tibble(t.l = numeric(),
                           t.u = numeric(),
                           L = numeric(),
                           U = numeric()) # What was user specified vs what was actually computed to achieve the simulation
  if(!nested) {
    bm$bb.local <- create.bm_(t, W_t, refine, mult, prefer, nested = TRUE) # Auxilliary info for Brownian bridge localised layers
  }
  class(bm) <- "BrownianMotion"
  bm
}
