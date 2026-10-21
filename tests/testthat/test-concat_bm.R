test_that("assigning a BrownianMotion aliases the same object", {
  withr::local_seed(81)
  bm <- local_bm_1d()
  bm.x <- bm
  sim(bm.x, 100)
  expect_equal(max(bm$t), 100)
  expect_identical(bm, bm.x)
})

test_that("concat_bm() of a single object makes an independent deep copy", {
  withr::local_seed(82)
  bm <- create_bm(c(0, 1, 2), c(0, 0.5, 0))
  layers(bm, 0, 2, type = "localised")
  t <- bm$t
  aux.t <- bm$bb.local$t
  bm.x <- concat_bm(bm)
  expect_s3_class(bm.x, "BrownianMotion")
  expect_false(identical(bm, bm.x))
  expect_false(identical(bm$bb.local, bm.x$bb.local))
  expect_equal(bm.x$t, bm$t)
  expect_equal(bm.x$layers, bm$layers)
  sim(bm.x, c(0.5, 10))
  expect_equal(bm$t, t)
  expect_equal(bm$bb.local$t, aux.t)
  expect_valid_bm(bm)
  expect_valid_bm(bm.x)
})

test_that("concat_bm() with t0 shifts all times, layers and labels", {
  withr::local_seed(83)
  bm <- create_bm(-1, 100)
  sim(bm, 0:2)
  layers(bm, 0, 2)
  add_labels(bm, t = 1, label = "one")
  shifted <- concat_bm(bm, t0 = 0)
  expect_equal(shifted$t, bm$t + 1)
  expect_equal(shifted$W_t, bm$W_t)
  expect_equal(shifted$layers$t.l, bm$layers$t.l + 1)
  expect_equal(shifted$labels$one, 2)
  expect_equal(shifted$labels$start, 0)
  expect_valid_bm(shifted)
})

test_that("concat_bm() joins paths end to end with a jump at the join", {
  withr::local_seed(84)
  a <- local_bm_1d(1:2)
  b <- create_bm(0, 10)
  sim(b, 1)
  ab <- concat_bm(a, b)
  expect_equal(ab$t, c(0, 1, 2, 3))
  expect_equal(unname(ab$W_t), unname(c(a$W_t[1:2], b$W_t)))
  expect_equal(unname(ab$W_tm), unname(c(a$W_tm, b$W_tm[2])))
  # The join at t = 2 is a jump
  expect_equal(ab$labels$seg.start, c(0, 2))
  expect_equal(ab$labels$seg.end, c(2, 3))
  expect_valid_bm(ab)
})

test_that("concat_bm() of layered paths keeps layers valid and simulable", {
  withr::local_seed(85)
  a <- local_bm_1d()
  layers(a, 0, 3)
  ab <- concat_bm(a, a, t0 = 1)
  expect_equal(nrow(ab$layers), 2 * nrow(a$layers))
  expect_equal(range(ab$t), c(1, 7))
  sim(ab, c(1.5, 4.5))
  expect_valid_bm(ab)
})

test_that("concat_bm() accepts a list of objects", {
  withr::local_seed(86)
  a <- local_bm_1d(1)
  ab <- concat_bm(list(a, a, a))
  expect_equal(ab$t, 0:3)
})

test_that("concat_bm() warns when settings differ and uses the first", {
  a <- create_bm()
  expect_warning(ab <- concat_bm(a, create_bm(mult = 2)), "same mult parameter")
  expect_equal(ab$mult, 1)
  expect_warning(concat_bm(a, create_bm(prefer = "intersection")),
                 "same prefer parameter")
  expect_warning(concat_bm(a, create_bm(refine = FALSE)), "same refine parameter")
})

test_that("concat_bm() validates its arguments", {
  expect_error(concat_bm(create_bm(), create_bm(dim = 2)),
               "All arguments must be BrownianMotion objects")
  expect_error(concat_bm(create_bm(dim = 2), create_bm(dim = 3)),
               "same dimension")
  expect_error(concat_bm(create_bm(dim = 2), local_bm_nd(2)),
               "same covariance")
})

test_that("concat_bm() on N-d objects concatenates every dimension", {
  withr::local_seed(87)
  bm <- local_bm_nd(2)
  sim(bm, 1:2)
  layers(bm, 0, 2)
  cc <- concat_bm(bm, bm, t0 = 5)
  expect_s3_class(cc, "BrownianMotionNd")
  expect_equal(cc[, "t"], c(bm[, "t"] + 5, bm[, "t"][-1] + 7))
  expect_valid_bm_nd(cc)
  # Copy, not alias
  sim(cc, 20)
  expect_equal(max(bm[, "t"]), 2)
})
