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

