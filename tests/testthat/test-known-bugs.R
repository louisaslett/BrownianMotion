# Reproductions of known bugs. Each is skipped so the suite passes; remove the
# skip() once the bug is fixed and the test should then pass unchanged.

test_that("intersection layer simulation never falls back to the broken Rmpfr path", {
  skip("Known bug: earpbrute uses pnorm() on Rmpfr numbers in intersection layers")
  # Tutorial steps; with this seed the second sim() hits the fallback and
  # errors with "Non-numeric argument to mathematical function".
  withr::local_seed(21)
  bm <- create_bm(prefer = "intersection")
  layers(bm, 1, 3, type = "bessel")
  sim(bm, 2, prefer = "bessel")
  expect_no_error(sim(bm, c(1.25, 2.5)))
  expect_valid_bm(bm)
})

# Suspected bugs found while writing this test suite. Each test asserts the
# behaviour I believe is intended; they are skipped pending a decision.

test_that("intersection layers on an interval with equal endpoint values terminate", {
  skip("Suspected bug: eabetaC_() loop in intersection.layers.sim_() never ends when W_s == W_t")
  # About 1 in 5 seeds loops forever in the `while(deind == 0)` loop because
  # eabetaC_() returns a negative lower bound whatever m is. Seed 20 hangs.
  # The same pattern in refine_intersection_() hangs the vignette workflow
  # for some seeds (e.g. set.seed(18) then the bm2 steps of the tutorial).
  withr::local_seed(20)
  bm <- create_bm(c(0, 1), c(0, 0))
  layers(bm, 0, 1, type = "intersection")
  expect_valid_bm(bm)
})

test_that("delete_skeleton(type = 'label') leaves layers untouched", {
  skip("Suspected bug: delete_skeleton() type = 'label' falls through and deletes layers")
  withr::local_seed(73)
  bm <- create_bm(mult = 3)
  layers(bm, 0, 10)
  sim(bm, seq(0, 10, by = 0.5))
  lyr <- bm$layers
  delete_skeleton(bm, l = 5, r = 7.5, type = "label")
  expect_equal(bm$layers, lyr)
})

test_that("[ by three or more unsorted indices returns values in request order", {
  skip("Suspected bug: `[` uses order(x) instead of order(order(x)) for index reordering")
  bm <- create_bm(c(-1, 0, 1), c(100, 101, 99))
  expect_equal(bm[c(2, 3, 1), "t"], c(0, 1, -1))
  expect_equal(bm[c(2, 3, 1)]$t, c(0, 1, -1))
  expect_equal(bm[c(2, 3, 1), "W_t"], c(101, 99, 100))
})

test_that("plot() of an N-d path without layers builds", {
  skip("Suspected bug: plot.BrownianMotionNd/NdZ error when there are no layers")
  withr::local_seed(1)
  bm <- local_bm_nd(2)
  sim(bm, 1:2)
  expect_plot_ok(plot(bm))
  expect_plot_ok(plot(bm$Z))
})

test_that("concat_bm() of N-d objects returns a fully formed BrownianMotionNd", {
  skip("Suspected bug: concat_bm.BrownianMotionNd() omits cov, chol.inv and Z")
  withr::local_seed(1)
  bm <- local_bm_nd(2)
  sim(bm, 1:2)
  cc <- concat_bm(bm, bm)
  expect_equal(cc$cov, bm$cov)
  expect_s3_class(cc$Z, "BrownianMotionNdZ")
  expect_no_error(add_points(cc, 10, c(0, 0)))
})
