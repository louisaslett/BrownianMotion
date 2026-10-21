layer_width <- function(lyr) pmax(lyr$Uu - lyr$Ud, lyr$Lu - lyr$Ld)

test_that("refine() returns bm invisibly and keeps the skeleton valid", {
  withr::local_seed(21)
  bm <- create_bm(refine = FALSE, mult = 3)
  layers(bm, 0, 10)
  res <- withVisible(refine(bm, 0, 10, mult = 1))
  expect_false(res$visible)
  expect_identical(res$value, bm)
  expect_valid_bm(bm)
})

test_that("refine() narrows bessel layers to the mult * sqrt(dt) criterion", {
  withr::local_seed(22)
  bm <- create_bm(refine = FALSE, mult = 3)
  layers(bm, 0, 10)
  sim(bm, seq(0, 10, by = 0.5))
  before <- bm$layers
  refine(bm, 5, 10, mult = 1)
  in.range <- bm$layers$t.l >= 5
  w <- layer_width(bm$layers[in.range, ])
  expect_true(all(w <= sqrt(bm$layers$t.u[in.range] - bm$layers$t.l[in.range]) + 1e-8))
  # Layers outside [s, t] are untouched
  expect_equal(bm$layers[!in.range, ], before[before$t.l < 5, ])
  # Refinement does not add skeleton points
  expect_equal(bm$t, seq(0, 10, by = 0.5), ignore_attr = TRUE)
  expect_valid_bm(bm)
})

test_that("refine() works on intersection layers", {
  withr::local_seed(23)
  bm <- create_bm(refine = FALSE, mult = 3, prefer = "intersection")
  layers(bm, 0, 4, type = "intersection")
  refine(bm, 0, 4, mult = 1)
  expect_equal(bm$layers$type, "intersection")
  w <- layer_width(bm$layers)
  expect_true(all(w <= 2 + 1e-8))
  expect_valid_bm(bm)
})

test_that("refine() works on localised and bridge layers", {
  withr::local_seed(24)
  bm <- create_bm(c(0, 1, 2), c(0, 0.5, 0), refine = FALSE)
  layers(bm, 0, 2, type = "localised")
  sim(bm, c(0.3, 1.5))
  refine(bm, 0, 2, mult = 0.5)
  expect_valid_bm(bm)
})

test_that("refine() with no layers in range is a no-op", {
  withr::local_seed(25)
  bm <- local_bm_1d()
  W <- bm$W_t
  refine(bm, 0, 3)
  expect_equal(bm$W_t, W)
  expect_equal(nrow(bm$layers), 0)
})
