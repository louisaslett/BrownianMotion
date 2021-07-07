#' @export
localise.layers <- function(bm, s, t, refine = bm$refine, mult = bm$mult, prefer = bm$prefer, label = c(names(s), names(t))) {
  layers(bm, s, t, "localised", refine, mult, prefer, label)
}
