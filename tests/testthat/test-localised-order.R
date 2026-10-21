# Regression tests: adding localised layers *before* existing localised bridge
# intervals used to corrupt the bounds of the existing intervals.

for (prefer in c("bessel", "intersection")) {
  for (split in c(FALSE, TRUE)) {
    test_that(paste0("localised layers added before existing bridge layers ",
                     "preserve them (prefer = ", prefer, ", split = ", split, ")"), {
      withr::local_seed(1)
      bm <- create_bm(c(0, 1, 2), c(0, 0, 0), prefer = prefer)
      layers(bm, 1, 2, type = "localised")
      if (split) sim(bm, 1.2)
      before <- bridge_bounds(bm)
      old.t <- bm$t
      old.values <- bm$W_t
      layers(bm, 0, 1, type = "localised")
      check_preserved(bm, before, old.t, old.values)

      sim(bm, c(0.25, 0.75, 1.75))
      refine(bm, 0, 2)
      check_offsets(bm)
    })
  }
}

test_that("filling a gap between bridge intervals preserves both neighbours", {
  # Ordinary layers to the right must not be adjusted, and the bridge interval
  # to the left must be preserved.
  withr::local_seed(1)
  bm <- create_bm(0:4, rep(0, 5))
  layers(bm, 0, 1, type = "localised")
  layers(bm, 2, 3, type = "localised")
  layers(bm, 3, 4, type = "intersection")
  ordinary <- bm$layers[bm$layers$t.l >= 3, ]
  before <- bridge_bounds(bm)
  old.t <- bm$t
  old.values <- bm$W_t
  layers(bm, 1, 2, type = "localised")
  check_preserved(bm, before, old.t, old.values)
  expect_identical(bm$layers[bm$layers$t.l >= 3, ], ordinary)
  sim(bm, c(1.25, 2.75))
  refine(bm, 0, 3)
  check_offsets(bm)
})
