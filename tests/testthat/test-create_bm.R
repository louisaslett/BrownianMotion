test_that("create_bm() defaults give a 1D path at W_0 = 0", {
  bm <- create_bm()
  expect_s3_class(bm, "BrownianMotion")
  expect_true(is.environment(bm))
  expect_identical(bm$t, 0)
  expect_equal(bm$W_t, 0)
  expect_equal(bm$W_tm, 0)
  expect_equal(nrow(bm$layers), 0)
  expect_equal(nrow(bm$user.layers), 0)
  expect_true(bm$refine)
  expect_equal(bm$mult, 1)
  expect_equal(bm$prefer, "bessel")
  expect_valid_bm(bm)
})

test_that("create_bm() sorts supplied times and keeps values paired", {
  bm <- create_bm(c(1, -1, 0), c(99, 100, 101))
  expect_equal(bm$t, c(-1, 0, 1))
  expect_equal(bm$W_t, c(100, 101, 99))
  expect_equal(bm$W_tm, bm$W_t)
  expect_equal(bm$labels$start, -1)
  expect_equal(bm$labels$end, 1)
  expect_equal(bm$labels$forced, c(-1, 0, 1))
  expect_valid_bm(bm)
})

test_that("create_bm() stores refine, mult and prefer settings", {
  bm <- create_bm(refine = FALSE, mult = 3, prefer = "Intersection")
  expect_false(bm$refine)
  expect_equal(bm$mult, 3)
  expect_equal(bm$prefer, "intersection")
})

test_that("create_bm() has an auxiliary bridge path with the same start", {
  bm <- create_bm(-1, 100)
  expect_s3_class(bm$bb.local, "BrownianMotion")
  expect_equal(nrow(bm$bb.local$layers), 0)
  expect_identical(parent.env(bm), emptyenv())
})

test_that("create_bm() validates t and W_t", {
  expect_error(create_bm(t = c(0, 0)), "Error with argument t")
  expect_error(create_bm(t = c(0, NA)), "Error with argument t")
  expect_error(create_bm(t = c(0, 1), W_t = 0), "Error with argument W_t")
  expect_error(create_bm(dim = 2, W_t = c(1, 2, 3)),
               "W_t: must be of length 1 or of length dim")
  expect_error(create_bm(dim = 2, t = c(0, 1), W_t = c(1, 2)),
               "Error with argument W_t")
})

test_that("create_bm() validates dim and cov", {
  expect_error(create_bm(dim = 0), "Error with argument dim")
  expect_error(create_bm(dim = 1.5), "Error with argument dim")
  expect_error(create_bm(cov = -1), "scalar cov must be strictly positive")
  expect_error(create_bm(dim = 2, cov = matrix(c(1, 2, 3, 4), 2)), "not symmetric")
  expect_error(create_bm(dim = 2, cov = matrix(c(1, 2, 2, 1), 2)),
               "not positive definite")
  expect_error(create_bm(dim = 2, cov = diag(3)), "Error with argument cov")
})

test_that("create_bm() validates refine, mult and prefer", {
  expect_error(create_bm(refine = "yes"), "Error with argument refine")
  expect_error(create_bm(refine = c(TRUE, FALSE)), "refine: must be length 1 or dim")
  expect_error(create_bm(mult = 0), "mult: must be strictly positive")
  expect_error(create_bm(mult = Inf), "Error with argument mult")
  expect_error(create_bm(prefer = "foo"), "prefer: must be one of")
  expect_error(create_bm(prefer = NA_character_), "Error with argument prefer")
})

test_that("create_bm() with dim > 1 or non-unit cov returns a BrownianMotionNd", {
  bm <- create_bm(dim = 3)
  expect_s3_class(bm, "BrownianMotionNd")
  expect_equal(bm$dim, 3L)
  expect_equal(bm$cov, diag(3))
  expect_s3_class(bm$Z, "BrownianMotionNdZ")
  expect_valid_bm_nd(bm)

  # A scalar cov other than 1 also uses the N-d representation
  expect_s3_class(create_bm(cov = 2), "BrownianMotionNd")
})

test_that("create_bm() N-d starting values: scalar, vector and matrix forms", {
  expect_equal(create_bm(dim = 3)[t = 0, "W_t"], matrix(0, 1, 3))
  expect_equal(create_bm(W_t = 3, dim = 3)[t = 0, "W_t"], matrix(3, 1, 3))
  expect_equal(create_bm(W_t = c(0, 1, 2), dim = 3)[t = 0, "W_t"],
               matrix(c(0, 1, 2), 1, 3))
  W <- matrix(c(0, 1, 2, 2, 0, 1), nrow = 2, ncol = 3, byrow = TRUE)
  bm <- create_bm(t = c(1, 2), W_t = W, dim = 3)
  expect_equal(bm[t = c(1, 2), "W_t"], W)
})

test_that("create_bm() N-d with correlated cov stores Cholesky factors", {
  cov <- diag(x = 0.5, nrow = 3) + 0.5
  bm <- create_bm(dim = 3, cov = cov, W_t = c(1, -1, 2))
  expect_equal(bm$cov, cov)
  expect_equal(crossprod(bm$chol), cov)
  expect_equal(bm$chol %*% bm$chol.inv, diag(3))
  # The stored (untransformed) start point round-trips through Z
  expect_equal(bm[t = 0, "W_t"], matrix(c(1, -1, 2), 1, 3))
})
