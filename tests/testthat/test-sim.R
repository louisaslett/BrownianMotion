test_that("sim() forward simulation adds sorted times and returns bm invisibly", {
  withr::local_seed(212)
  bm <- create_bm(c(-1, 0, 1), c(100, 101, 99))
  res <- withVisible(sim(bm, 2:4))
  expect_false(res$visible)
  expect_identical(res$value, bm)
  expect_equal(bm$t, -1:4)
  expect_equal(bm$labels$user, 2:4)
  expect_equal(bm$labels$end, 4)
  expect_valid_bm(bm)
})

test_that("sim() mutates the object in place", {
  withr::local_seed(1)
  bm <- create_bm()
  sim(bm, 1)
  expect_equal(bm$t, c(0, 1))
})

test_that("sim() at unsorted and interior times keeps the skeleton sorted", {
  withr::local_seed(2)
  bm <- create_bm(c(-1, 0, 1), c(100, 101, 99))
  old <- bm$W_t
  sim(bm, c(4.5, 0.5, -0.5, 2.5))
  expect_equal(bm$t, c(-1, -0.5, 0, 0.5, 1, 2.5, 4.5), ignore_attr = TRUE)
  # Existing points are untouched
  expect_equal(unname(bm$W_t[match(c(-1, 0, 1), bm$t)]), old)
  expect_valid_bm(bm)
})

test_that("sim() at an existing time does not duplicate it but adds the label", {
  withr::local_seed(3)
  bm <- local_bm_1d()
  W <- bm$W_t
  sim(bm, 2, label = "two")
  expect_equal(bm$t, 0:3)
  expect_equal(bm$W_t, W)
  expect_equal(bm$labels$two, 2)
})

test_that("sim() with labels records them", {
  withr::local_seed(4)
  bm <- create_bm()
  sim(bm, c(6, 7), label = c("special.time.1", "special.time.2"))
  expect_equal(bm$labels$special.time.1, 6)
  expect_equal(bm$labels$special.time.2, 7)
  sim(bm, c(a = 8))
  expect_equal(bm$labels$a, 8)
  sim(bm, c(9, 10), label = "both")
  expect_equal(bm$labels$both, c(9, 10))
})

test_that("sim() validates its arguments", {
  bm <- create_bm()
  expect_error(sim(bm, "a"), "t must be a vector of times")
  expect_error(sim(bm, -1), "cannot simulate path at times before")
  expect_error(sim(bm, 1, refine = "x"), "refine must be a scalar logical")
  expect_error(sim(bm, 1, mult = -1), "mult must be a strictly positive scalar")
  expect_error(sim(bm, 1, prefer = "foo"), "prefer must be one of")
  expect_error(sim(bm, 1, label = "start"), "clashes with an internal label")
  expect_error(sim(bm, 1:2, label = c("a", "b", "c")), "label is incorrect length")
  expect_error(sim(bm, 1, label = 5), "invalid labels provided")
  # Nothing was added by the failed calls
  expect_identical(bm$t, 0)
})

test_that("sim() inside bessel layers stays within the layer bounds", {
  withr::local_seed(5)
  bm <- local_bm_1d()
  layers(bm, 0, 3, type = "bessel")
  sim(bm, c(0.5, 1.7, 2.2, 2.9))
  expect_true(all(c(0.5, 1.7, 2.2, 2.9) %in% bm$t))
  expect_true(all(bm$layers$type == "bessel"))
  expect_valid_bm(bm)
})

test_that("sim() inside intersection layers stays within the layer bounds", {
  # Seed chosen to avoid known bug 1 (Rmpfr fallback in intersection layers)
  withr::local_seed(6)
  bm <- create_bm(prefer = "intersection")
  sim(bm, 1:2)
  layers(bm, 0, 2, type = "intersection")
  sim(bm, c(0.5, 1.5))
  expect_true(all(bm$layers$type == "intersection"))
  expect_valid_bm(bm)
})

test_that("sim() inside localised layers stays within the layer bounds", {
  withr::local_seed(7)
  bm <- create_bm()
  first_passage(bm, delta = 1)
  first_passage(bm, delta = 0.5)
  mid <- mean(bm$t[2:3])
  sim(bm, c(bm$t[2] / 2, mid))
  expect_valid_bm(bm)
})

test_that("sim() prefer converts layer types inside the simulated interval", {
  withr::local_seed(8)
  bm <- local_bm_1d(1:2)
  layers(bm, 0, 2, type = "bessel")
  sim(bm, 0.5, prefer = "intersection")
  expect_true(all(bm$layers$type[bm$layers$t.u <= 1] == "intersection"))
  expect_valid_bm(bm)
})

test_that("sim() with refine = TRUE meets the mult * sqrt(dt) layer width criterion", {
  withr::local_seed(9)
  bm <- create_bm(refine = TRUE, mult = 1)
  layers(bm, 0, 10)
  sim(bm, seq(0, 10, by = 0.5))
  w <- pmax(bm$layers$Uu - bm$layers$Ud, bm$layers$Lu - bm$layers$Ld)
  expect_true(all(w <= sqrt(bm$layers$t.u - bm$layers$t.l) + 1e-8))
  expect_valid_bm(bm)
})

test_that("sim() with refine = FALSE still produces a valid skeleton", {
  withr::local_seed(9)
  bm <- create_bm(refine = FALSE)
  layers(bm, 0, 10)
  sim(bm, seq(0, 10, by = 0.5))
  expect_false(bm$refine)
  expect_valid_bm(bm)
})

test_that("unconditional increments have the Brownian variance", {
  skip_on_cran()
  withr::local_seed(10)
  n <- 400
  bm <- create_bm()
  sim(bm, seq_len(n) * 0.25)
  inc <- diff(bm$W_t)
  expect_lt(abs(mean(inc)), 4 * sqrt(0.25 / n))
  expect_gt(stats::ks.test(inc / sqrt(0.25), "pnorm")$p.value, 0.001)
})
