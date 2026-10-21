test_that("add_points() forces points into the skeleton (vignette example)", {
  bm <- create_bm(-1, 100)
  res <- withVisible(add_points(bm, c(0, 1), c(101, 99)))
  expect_false(res$visible)
  expect_identical(res$value, bm)
  expect_equal(bm$t, c(-1, 0, 1))
  expect_equal(bm$W_t, c(100, 101, 99))
  expect_equal(bm$W_tm, bm$W_t)
  expect_equal(bm$labels$forced, c(-1, 0, 1))
  expect_equal(bm$labels$user, c(0, 1))
  expect_valid_bm(bm)
})

test_that("add_points() before the start updates start and end labels", {
  withr::local_seed(41)
  bm <- local_bm_1d()
  add_points(bm, -2, 5)
  expect_equal(bm$t[1], -2)
  expect_equal(bm$labels$start, -2)
  expect_equal(bm$labels$seg.start[1], -2)
  add_points(bm, 10, 0)
  expect_equal(bm$labels$end, 10)
  expect_valid_bm(bm)
})

test_that("add_points() inserts single and multiple interior points in order", {
  bm <- create_bm(c(0, 10), c(0, 0))
  add_points(bm, 5, 1)
  add_points(bm, c(7, 2), c(3, 2))
  expect_equal(bm$t, c(0, 2, 5, 7, 10))
  expect_equal(bm$W_t, c(0, 2, 1, 3, 0))
  expect_valid_bm(bm)
})

test_that("add_points() applies labels", {
  bm <- create_bm()
  add_points(bm, c(1, 2), c(0, 0), label = "pt")
  expect_equal(bm$labels$pt, c(1, 2))
  add_points(bm, c(three = 3), 1)
  expect_equal(bm$labels$three, 3)
})

test_that("add_points() validates its arguments", {
  withr::local_seed(42)
  bm <- local_bm_1d(1)
  expect_error(add_points(bm, 1, 0), "already a time in the skeleton")
  expect_error(add_points(bm, c(5, 6), 0), "must be the same length")
  expect_error(add_points(bm, "a", 0), "t must be a vector of times")
  expect_error(add_points(bm, 5, "a"), "W_t must be a vector")
  expect_error(add_points(bm, 5, 0, label = "end"), "clashes")
  layers(bm, 0, 1)
  expect_error(add_points(bm, 0.5, 0), "within a layer")
  expect_equal(bm$t, c(0, 1))
})

test_that("add_points() on N-d objects updates every dimension", {
  bm <- create_bm(t = c(1, 2), W_t = matrix(c(0, 1, 2, 2, 0, 1), 2, 3, byrow = TRUE),
                  dim = 3)
  add_points(bm, t = 0, W_t = c(-1, -1, -1))
  expect_equal(bm[t = c(0, 1, 2), "W_t"],
               matrix(c(-1, -1, -1, 0, 1, 2, 2, 0, 1), 3, 3, byrow = TRUE))
  expect_valid_bm_nd(bm)
  expect_error(add_points(bm, t = 5, W_t = c(1, 2)), "length 1 or of length dim")
  expect_error(add_points(bm, t = 5:6, W_t = c(1, 2, 3)), "Error with argument W_t")
})

test_that("add_points() on N-d objects with correlated cov round-trips", {
  cov <- matrix(c(1, 0.5, 0.5, 1), 2)
  bm <- create_bm(dim = 2, cov = cov)
  add_points(bm, t = c(1, 2), W_t = matrix(c(1, 2, 3, 4), 2, 2))
  expect_equal(bm[t = c(1, 2), "W_t"], matrix(c(1, 2, 3, 4), 2, 2))
  expect_valid_bm_nd(bm)
})
