#' @importFrom tibble tibble add_row
#' @importFrom glue glue

is.realscalar <- function(x) {
  is.atomic(x) && (length(x) == 1L) && !is.character(x) && !is.raw(x) && !is.logical(x) && (Im(x) == 0)
}
is.intscalar <- function(x) {
  is.atomic(x) && (length(x) == 1L) && !is.character(x) && !is.raw(x) && !is.logical(x) && (Im(x) == 0) && all(x - round(x) == 0)
}
is.realvector <- function(x) {
  is.atomic(x) && !is.character(x) && !is.raw(x) && !is.logical(x) && all(Im(x) == 0)
}
