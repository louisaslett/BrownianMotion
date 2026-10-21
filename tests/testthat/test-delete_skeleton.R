local_bm_layered <- function() {
  bm <- create_bm(mult = 3)
  layers(bm, 0, 10)
  sim(bm, seq(0, 10, by = 0.5))
  bm
}

test_that("delete_skeleton() type = 'all' removes points, layers and labels inside", {
  withr::local_seed(71)
  bm <- local_bm_layered()
  add_labels(bm, t = 3, label = "three")
  res <- withVisible(delete_skeleton(bm, l = 2.5, r = 5))
  expect_false(res$visible)
  expect_identical(res$value, bm)
  expect_false(any(bm$t > 2.5 & bm$t < 5))
  expect_true(all(c(2.5, 5) %in% bm$t))
  expect_false(any(bm$layers$t.l >= 2.5 & bm$layers$t.u <= 5))
  expect_null(bm$labels$three)
  expect_valid_bm(bm)
})

test_that("delete_skeleton() type = 'layer' keeps the points", {
  withr::local_seed(72)
  bm <- local_bm_layered()
  t <- bm$t
  delete_skeleton(bm, l = 5, r = 7.5, type = "layer")
  expect_equal(bm$t, t)
  expect_false(any(bm$layers$t.l >= 5 & bm$layers$t.u <= 7.5))
  expect_true(any(bm$layers$t.l >= 7.5))
  expect_valid_bm(bm)
})

test_that("delete_skeleton() type = 'label' removes user labels but keeps points", {
  withr::local_seed(73)
  bm <- local_bm_layered()
  add_labels(bm, t = c(5.5, 8), label = "lab")
  t <- bm$t
  delete_skeleton(bm, l = 5, r = 7.5, type = "label")
  expect_equal(bm$labels$lab, 8)
  expect_equal(bm$t, t)
  # Layers should also be kept: see the suspected bug in test-known-bugs.R
})

test_that("delete_skeleton() with only l or r trims the ends", {
  withr::local_seed(74)
  bm <- local_bm_layered()
  delete_skeleton(bm, l = 8)
  expect_equal(max(bm$t), 8)
  expect_equal(bm$labels$end, 8)
  delete_skeleton(bm, r = 1)
  expect_equal(min(bm$t), 1)
  expect_equal(bm$labels$start, 1)
  expect_valid_bm(bm)
})

test_that("delete_skeleton() splits a user layer that straddles the range", {
  withr::local_seed(75)
  bm <- create_bm()
  first_passage(bm, delta = 2)
  tt <- bm$t[2]
  sim(bm, tt * c(0.2, 0.4, 0.6, 0.8))
  delete_skeleton(bm, l = bm$t[3], r = bm$t[4])
  ul <- bm$user.layers[order(bm$user.layers$t.l), ]
  expect_equal(nrow(ul), 2)
  expect_equal(ul$t.u[1], bm$t[3])
  expect_equal(ul$t.l[2], bm$t[4])
})

test_that("delete_skeleton() works on localised bridge layers", {
  withr::local_seed(76)
  bm <- create_bm(c(0, 1, 2), c(0, 0.5, 0))
  layers(bm, 0, 2, type = "localised")
  delete_skeleton(bm, l = 1, r = 2)
  expect_equal(bm$t[bm$t >= 1], c(1, 2))
  expect_false(any(bm$bb.local$t > 1 & bm$bb.local$t < 2))
  expect_valid_bm(bm)
})

test_that("delete_skeleton() validates its arguments", {
  withr::local_seed(77)
  bm <- local_bm_1d()
  expect_error(delete_skeleton(bm), "at least one of l or r must be specified")
  expect_error(delete_skeleton(bm, l = 0.5), "l must be a sampled path point")
  expect_error(delete_skeleton(bm, r = 0.5), "r must be a sampled path point")
  expect_error(delete_skeleton(bm, l = 2, r = 1), "l must be less than r")
  expect_error(delete_skeleton(bm, l = 1, type = "foo"), "type must be either")
  expect_error(delete_skeleton(bm, l = "a"), "l must be a scalar")
  expect_equal(bm$t, 0:3)
})

test_that("delete_skeleton() on N-d objects applies to all dimensions", {
  withr::local_seed(78)
  bm <- local_bm_nd(3)
  sim(bm, c(0.5, 5))
  layers(bm, 0, 10)
  delete_skeleton(bm, l = 5)
  expect_equal(max(bm[, "t"]), 5)
  expect_valid_bm_nd(bm)
})
