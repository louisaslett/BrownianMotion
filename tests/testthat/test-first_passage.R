test_that("first_passage() with delta exits a symmetric band", {
  withr::local_seed(31)
  bm <- create_bm()
  res <- withVisible(first_passage(bm, delta = 1))
  expect_false(res$visible)
  expect_identical(res$value, bm)
  expect_length(bm$t, 2)
  expect_true(abs(bm$W_t[2]) == 1)
  expect_equal(bm$layers$type, "localised")
  expect_equal(bm$user.layers$L, -1)
  expect_equal(bm$user.layers$U, 1)
  expect_equal(bm$labels$fpt, bm$t[2])
  expect_equal(bm$labels$end, bm$t[2])
  expect_valid_bm(bm)
})

test_that("first_passage() defaults to delta = 1", {
  withr::local_seed(32)
  bm <- create_bm(0, 5)
  first_passage(bm)
  expect_true(unname(bm$W_t[2]) %in% c(4, 6))
})

test_that("first_passage() with asymmetric l and u ends on a boundary", {
  withr::local_seed(33)
  bm <- create_bm()
  first_passage(bm, l = -1, u = 3)
  end <- unname(bm$W_t[length(bm$W_t)])
  expect_true(isTRUE(all.equal(end, -1)) || isTRUE(all.equal(end, 3)))
  expect_true(all(bm$W_t >= -1 - 1e-12 & bm$W_t <= 3 + 1e-12))
  expect_true(all(bm$layers$type == "localised"))
  expect_equal(bm$user.layers$L, -1)
  expect_equal(bm$user.layers$U, 3)
  expect_valid_bm(bm)
})

test_that("first_passage() with delta.l and delta.u is relative to the endpoint", {
  withr::local_seed(34)
  bm <- create_bm(0, 10)
  first_passage(bm, delta.l = -1, delta.u = 3)
  expect_equal(bm$user.layers$L, 9)
  expect_equal(bm$user.layers$U, 13)
  expect_valid_bm(bm)
})

test_that("first_passage() appends after existing path and applies label", {
  withr::local_seed(35)
  bm <- local_bm_1d()
  first_passage(bm, delta = 0.5, label = "hit")
  expect_gt(max(bm$t), 3)
  expect_equal(bm$labels$hit, max(bm$t))
  expect_equal(bm$layers$t.l[1], 3)
  expect_valid_bm(bm)
})

test_that("sim() inside a first passage layer respects its bounds", {
  withr::local_seed(36)
  bm <- create_bm()
  first_passage(bm, l = -1, u = 2)
  sim(bm, seq(0, max(bm$t), length.out = 7)[2:6])
  expect_valid_bm(bm)
})

test_that("first_passage() validates its arguments", {
  bm <- local_bm_1d(1)
  expect_error(first_passage(bm, delta = -1), "must be a positive real scalar")
  expect_error(first_passage(bm, delta = 1, l = -5), "other limit arguments")
  expect_error(first_passage(bm, delta.l = -1, u = 5), "other limit arguments")
  expect_error(first_passage(bm, delta.l = 1), "must be a negative real scalar")
  expect_error(first_passage(bm, delta.u = -1), "must be a positive real scalar")
  x <- bm$W_t[2]
  expect_error(first_passage(bm, l = x + 1, u = x + 2), "l must be less than")
  expect_error(first_passage(bm, l = x - 2, u = x - 1), "u must be greater than")
  expect_error(first_passage(bm, l = x + 1, u = x - 1), "l cannot be greater")
  expect_error(first_passage(bm, delta = 1, label = "end"), "clashes")
  expect_error(first_passage(list(), delta = 1), "must be a BrownianMotion object")
  expect_length(bm$t, 2)
})

test_that("first_passage() errors on multi-dimensional objects", {
  bm <- create_bm(dim = 2)
  expect_error(first_passage(bm, delta = 1), "only 1-dimensional")
})

test_that("first passage exit side is fair for a symmetric band", {
  skip_on_cran()
  withr::local_seed(37)
  n <- 400
  up <- vapply(seq_len(n), function(i) {
    bm <- create_bm()
    first_passage(bm, delta = 1)
    bm$W_t[2] > 0
  }, logical(1))
  expect_gt(stats::binom.test(sum(up), n)$p.value, 0.001)
})

test_that("first_passage() names the offending bound in its error", {
  bm <- create_bm()
  expect_error(first_passage(bm, l = c(-1, -2), u = 1), "^l must specify")
  expect_error(first_passage(bm, l = -1, u = c(1, 2)), "^u must specify")
})
