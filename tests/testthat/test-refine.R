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

test_that("refine() validates its arguments", {
  # Previously a non-positive mult made refine() loop forever, and invalid s or
  # t silently did nothing.
  withr::local_seed(81)
  bm <- create_bm()
  sim(bm, 1:3)
  layers(bm, 0, 3, type = "bessel")
  before <- bm$layers

  expect_error(refine(bm, "a", 2), "s, t and mult must be real scalars")
  expect_error(refine(bm, 0, c(1, 2)), "s, t and mult must be real scalars")
  expect_error(refine(bm, NA, 2), "s, t and mult must be real scalars")
  expect_error(refine(bm, 0, 3, mult = "a"), "s, t and mult must be real scalars")
  expect_error(refine(bm, 2, 1), "s must be less than t")
  expect_error(refine(bm, 0, 3, mult = 0), "mult must be strictly positive and finite")
  expect_error(refine(bm, 0, 3, mult = -1), "mult must be strictly positive and finite")
  expect_error(refine(bm, 0, 3, mult = Inf), "mult must be strictly positive and finite")
  expect_identical(bm$layers, before)

  # Infinite limits are allowed and refine every layer
  expect_no_error(refine(bm, -Inf, Inf, mult = 0.3))
  expect_valid_bm(bm)

  nd <- create_bm(dim = 2)
  sim(nd, 1:2)
  layers(nd, 0, 2)
  expect_error(refine(nd, 0, 2, mult = 0), "mult must be strictly positive and finite")
})
