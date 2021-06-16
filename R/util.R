#' @importFrom tibble tibble add_row
#' @importFrom glue glue
#' @importFrom Rmpfr mpfr

is.realscalar <- function(x) {
  is.atomic(x) && (length(x) == 1L) && !is.character(x) && !is.raw(x) && !is.logical(x) && (Im(x) == 0)
}
is.intscalar <- function(x) {
  is.atomic(x) && (length(x) == 1L) && !is.character(x) && !is.raw(x) && !is.logical(x) && (Im(x) == 0) && all(x - round(x) == 0)
}
is.realvector <- function(x) {
  is.atomic(x) && !is.character(x) && !is.raw(x) && !is.logical(x) && all(Im(x) == 0)
}
is.bm <- function(x) {
  is.environment(x) && ("BrownianMotion" %in% class(bm))
}
is.timevector <- function(x) {
  checkmate::test_numeric(x,
                          lower = 0,
                          finite = TRUE,
                          any.missing = FALSE,
                          min.len = 1,
                          unique = TRUE)
}


test.timevector <- function(x) {
  checkmate::test_numeric(x,
                          lower = 0,
                          finite = TRUE,
                          any.missing = FALSE,
                          min.len = 1,
                          unique = TRUE)
}
assert.timevector <- function(x) {
  if(!test.timevector(x)) {
    stop("invalid vector of times")
  }
}

assert.bmlabel <- function(x, t) {
  if(is.null(x)) {
    return()
  }
  if(length(x)!=1 && length(x)!=length(t)) {
    stop("label is incorrect length")
  }
  if(!checkmate::test_character(x)) {
    stop("invalid labels provided")
  }
}

add.labels_ <- function(bm, labels, t) {
  if(is.null(labels)) {
    return()
  } else if(length(labels) == 1) {
    bm$labels[[labels]] <- unique(c(bm$labels[[labels]], unname(t)))
  } else if(length(labels) == length(t)) {
    for(l in unique(labels)) {
      bm$labels[[l]] <- unique(c(bm$labels[[l]], unname(t[labels==l])))
    }
  } else {
    stop("invalid specification to add.labels_")
  }
}