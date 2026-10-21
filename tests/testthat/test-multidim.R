test_that("sim() on N-d objects returns bm invisibly and fills every dimension", {
  withr::local_seed(121)
  bm <- local_bm_nd(3)
  res <- withVisible(sim(bm, c(0.5, 5)))
  expect_false(res$visible)
  expect_identical(res$value, bm)
  expect_equal(bm[, "t"], c(0, 0.5, 5))
  expect_valid_bm_nd(bm)
})

test_that("layers() on N-d objects gives all dimensions the same skeleton", {
  withr::local_seed(122)
  bm <- local_bm_nd(3)
  sim(bm, c(0.5, 5))
  layers(bm, 0, 10)
  expect_valid_bm_nd(bm)
  for (d in 1:3) {
    lyr <- bm$Z.bm[[d]]$layers
    expect_equal(min(lyr$t.l), 0)
    expect_equal(max(lyr$t.u), 10)
  }
})

test_that("N-d outer layer bounds contain the path", {
  withr::local_seed(123)
  bm <- local_bm_nd(2)
  sim(bm, 1:3)
  layers(bm, 0, 3)
  sim(bm, c(0.5, 2.5))
  W <- bm[, "W_t"]
  t <- bm[, "t"]
  lyr <- bm[, "layers"]
  for (i in seq_len(nrow(lyr))) {
    inside <- t >= lyr$t.l[i] & t <= lyr$t.u[i]
    for (d in 1:2) {
      expect_true(all(W[inside, d] >= lyr$outer.L[i, d] - 1e-8 &
                        W[inside, d] <= lyr$outer.U[i, d] + 1e-8))
    }
  }
})

test_that("N-d workflow from the vignette keeps a valid skeleton", {
  withr::local_seed(124)
  bm <- create_bm(dim = 3, cov = diag(x = 0.5, nrow = 3) + 0.5)
  sim(bm, c(0.5, 5))
  layers(bm, 0, 10)
  refine(bm, s = 0, t = 10, mult = 0.1)
  delete_skeleton(bm, l = 5)
  add_segment(bm, l = 5, W_t = c(12, 10, 10))
  layers(bm, 5, 10)
  expect_equal(bm[t = 5, "W_t"], matrix(c(12, 10, 10), 1))
  expect_equal(max(bm[, "t"]), 10)
  expect_valid_bm_nd(bm)
})

test_that("independent N-d dimensions are each Brownian", {
  skip_on_cran()
  withr::local_seed(125)
  n <- 300
  bm <- create_bm(dim = 2, cov = diag(c(1, 4)))
  sim(bm, seq_len(n))
  inc <- diff(bm[, "W_t"])
  expect_gt(stats::ks.test(inc[, 1], "pnorm")$p.value, 0.001)
  expect_gt(stats::ks.test(inc[, 2], "pnorm", sd = 2)$p.value, 0.001)
})

test_that("N-d multi-dim objects reject per-dim settings of the wrong length", {
  expect_error(create_bm(dim = 3, refine = c(TRUE, FALSE)), "must be length 1 or dim")
  expect_error(create_bm(dim = 3, mult = c(1, 2)), "must be length 1 or dim")
  expect_error(create_bm(dim = 3, prefer = c("bessel", "bessel")),
               "must be length 1 or dim")
})

test_that("per-dimension refine, mult and prefer are passed to each component", {
  bm <- create_bm(dim = 2, refine = c(TRUE, FALSE), mult = c(1, 2),
                  prefer = c("bessel", "intersection"))
  expect_true(bm$Z.bm[[1]]$refine)
  expect_false(bm$Z.bm[[2]]$refine)
  expect_equal(bm$Z.bm[[2]]$mult, 2)
  expect_equal(bm$Z.bm[[2]]$prefer, "intersection")
})
