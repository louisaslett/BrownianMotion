test_that("add_labels() adds labels at existing times", {
  withr::local_seed(61)
  bm <- local_bm_1d()
  res <- withVisible(add_labels(bm, t = c(1, 3), label = "mine"))
  expect_false(res$visible)
  expect_identical(res$value, bm)
  expect_equal(bm$labels$mine, c(1, 3))
  add_labels(bm, t = 2, label = "mine")
  expect_setequal(bm$labels$mine, 1:3)
  add_labels(bm, t = c(a = 0))
  expect_equal(bm$labels$a, 0)
})

test_that("add_labels() pairs labels with times when both are vectors", {
  withr::local_seed(62)
  bm <- local_bm_1d()
  add_labels(bm, t = c(1, 2), label = c("one", "two"))
  expect_equal(bm$labels$one, 1)
  expect_equal(bm$labels$two, 2)
})

test_that("add_labels() validates its arguments", {
  withr::local_seed(63)
  bm <- local_bm_1d()
  expect_error(add_labels(bm, 0.5, "x"), "existing observed time points")
  expect_error(add_labels(bm, 1, "user"), "clashes")
  expect_error(add_labels(bm, 1:2, c("a", "b", "c")), "incorrect length")
})

test_that("delete_labels() by label removes it everywhere", {
  withr::local_seed(64)
  bm <- local_bm_1d()
  add_labels(bm, t = 1:3, label = "x")
  res <- withVisible(delete_labels(bm, label = "x"))
  expect_false(res$visible)
  expect_null(bm$labels$x)
  expect_equal(bm$labels$user, 1:3)
})

test_that("delete_labels() by time and label removes only those entries", {
  withr::local_seed(65)
  bm <- local_bm_1d(1:7)
  add_labels(bm, t = c(2, 3, 6, 7), label = "special")
  delete_labels(bm, t = c(7, 2, 3), label = "special")
  expect_equal(bm$labels$special, 6)
  delete_labels(bm, t = 6, label = "special")
  expect_null(bm$labels$special)
})

test_that("delete_labels() by time only removes all user labels there", {
  withr::local_seed(66)
  bm <- local_bm_1d()
  add_labels(bm, t = 1, label = "a")
  add_labels(bm, t = c(1, 2), label = "b")
  delete_labels(bm, t = 1)
  expect_null(bm$labels$a)
  expect_equal(bm$labels$b, 2)
  # Internal labels are kept
  expect_equal(bm$labels$user, 1:3)
})

test_that("delete_labels() validates its arguments", {
  withr::local_seed(67)
  bm <- local_bm_1d()
  expect_error(delete_labels(bm), "must specify at least one of t and label")
  expect_error(delete_labels(bm, label = "user"), "clashes")
  expect_error(delete_labels(bm, t = 0.5, label = "x"), "existing observed time points")
})

test_that("labels work for N-d objects across all dimensions", {
  withr::local_seed(68)
  bm <- local_bm_nd(2)
  sim(bm, 1:2)
  add_labels(bm, t = 2, label = "two")
  for (d in 1:2) expect_equal(bm$Z.bm[[d]]$labels$two, 2)
  expect_equal(bm["two", "t"], 2)
  delete_labels(bm, label = "two")
  for (d in 1:2) expect_null(bm$Z.bm[[d]]$labels$two)
})
