local_bm_grid <- function() {
  bm <- create_bm()
  sim(bm, seq(0, 10, length = 11))
  bm
}

test_that("add_segment() shifts a segment and creates jumps at its ends", {
  withr::local_seed(51)
  bm <- local_bm_grid()
  W <- bm$W_t
  res <- withVisible(add_segment(bm, l = 1, r = 3, delta = 1, label = "path.move"))
  expect_false(res$visible)
  expect_identical(res$value, bm)
  # W_t is shifted on [1, 3); W_tm is shifted on (1, 3]
  expect_equal(bm$W_t, W + c(0, 1, 1, rep(0, 8)))
  expect_equal(bm$W_tm, W + c(0, 0, 1, 1, rep(0, 7)))
  expect_equal(which(bm$W_t != bm$W_tm), c(2L, 4L))
  expect_equal(bm$labels$seg.start, c(0, 1, 3))
  expect_equal(bm$labels$seg.end, c(1, 3, 10))
  expect_equal(bm$labels$path.move, c(1, 3))
  expect_valid_bm(bm)
})

test_that("add_segment() with W_t sets the value at l", {
  withr::local_seed(52)
  bm <- local_bm_grid()
  add_segment(bm, l = 10, W_t = 0)
  expect_equal(unname(bm$W_t[bm$t == 10]), 0)
  sim(bm, 11)
  expect_valid_bm(bm)
})

test_that("add_segment() with only r shifts everything before r", {
  withr::local_seed(53)
  bm <- local_bm_grid()
  W <- bm$W_t
  add_segment(bm, r = 5, W_t = 4)
  expect_equal(unname(bm$W_t[1]), 4)
  expect_equal(bm$W_t[bm$t > 5], W[bm$t > 5])
  expect_equal(bm$W_tm[6], W[6] + 4 - W[1])
})

test_that("add_segment() over the whole path shifts everything, including layers", {
  withr::local_seed(54)
  bm <- local_bm_grid()
  layers(bm, 0, 10)
  W <- bm$W_t
  lyr <- bm$layers
  add_segment(bm, delta = -1)
  expect_equal(bm$W_t, W - 1)
  expect_equal(bm$W_tm, W - 1)
  expect_equal(bm$layers$Ld, lyr$Ld - 1)
  expect_equal(bm$layers$Uu, lyr$Uu - 1)
  expect_valid_bm(bm)
})

test_that("add_segment() then delete_skeleton() at the jumps removes them", {
  withr::local_seed(55)
  bm <- local_bm_grid()
  add_segment(bm, l = 1, r = 3, delta = 1)
  delete_skeleton(bm, l = 1, r = 3)
  expect_equal(bm$W_t, bm$W_tm)
  expect_equal(bm$labels$seg.start, 0)
  expect_valid_bm(bm)
})

test_that("sim() on a path with jumps respects the segments", {
  withr::local_seed(56)
  bm <- local_bm_grid()
  add_segment(bm, l = 1, r = 3, delta = 1)
  sim(bm, seq(0, 4, length = 21))
  expect_equal(sum(bm$W_t != bm$W_tm), 2)
  layers(bm, 0, 4)
  expect_valid_bm(bm)
})

test_that("add_segment() validates its arguments", {
  withr::local_seed(57)
  bm <- local_bm_1d(1:2)
  expect_error(add_segment(bm), "One of delta or W_t must be specified")
  expect_error(add_segment(bm, delta = 1, W_t = 1), "Can't specify both")
  expect_error(add_segment(bm, l = 0.5, delta = 1), "l must be a sampled path point")
  expect_error(add_segment(bm, r = 0.5, delta = 1), "r must be a sampled path point")
  expect_error(add_segment(bm, l = 2, r = 1, delta = 1), "l must be less than r")
  expect_error(add_segment(bm, l = "a", delta = 1), "l must be a scalar")
  expect_error(add_segment(bm, l = 1, delta = 1, label = "start"), "clashes")
})

test_that("add_segment() on N-d objects shifts every coordinate", {
  withr::local_seed(58)
  bm <- local_bm_nd(3)
  sim(bm, 1:4)
  W <- bm[, "W_t"]
  add_segment(bm, l = 2, W_t = c(12, 10, 10))
  expect_equal(bm[t = 2, "W_t"], matrix(c(12, 10, 10), 1))
  expect_equal(bm[t = 4, "W_t"], W[5, , drop = FALSE] - W[3, , drop = FALSE] +
                 matrix(c(12, 10, 10), 1))
  expect_equal(bm[t = 1, "W_t"], W[2, , drop = FALSE])
  expect_valid_bm_nd(bm)
  expect_error(add_segment(bm, l = 2, delta = c(1, 2)), "length 1 or of length dim")
})
