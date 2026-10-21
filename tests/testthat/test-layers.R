test_that("layers() bessel creates one layer per skeleton interval", {
  withr::local_seed(11)
  bm <- local_bm_1d()
  res <- withVisible(layers(bm, 0, 3, type = "bessel"))
  expect_false(res$visible)
  expect_identical(res$value, bm)
  expect_equal(bm$layers$type, rep("bessel", 3))
  expect_equal(bm$layers$t.l, 0:2)
  expect_equal(bm$layers$t.u, 1:3)
  expect_named(bm$layers, c("type", "t.l", "t.u", "Ld", "Uu", "Lu", "Ud",
                            "Lu.hard", "Ud.hard"))
  expect_valid_bm(bm)
})

test_that("layers() simulates endpoints that are not yet in the skeleton", {
  withr::local_seed(12)
  bm <- local_bm_1d()
  layers(bm, 10, 12)
  expect_true(all(c(10, 12) %in% bm$t))
  expect_equal(bm$layers$t.l, 10)
  expect_equal(bm$layers$t.u, 12)
  expect_valid_bm(bm)
})

test_that("layers() intersection creates intersection layers", {
  withr::local_seed(13)
  bm <- local_bm_1d()
  layers(bm, 0, 3, type = "intersection")
  expect_equal(bm$layers$type, rep("intersection", 3))
  expect_valid_bm(bm)
})

test_that("layers() can convert between bessel and intersection", {
  withr::local_seed(14)
  bm <- local_bm_1d(1:2)
  layers(bm, 0, 2, type = "bessel")
  layers(bm, 0, 2, type = "intersection")
  expect_equal(bm$layers$type, rep("intersection", 2))
  expect_valid_bm(bm)
  # Converting back is only possible for intersection layers whose inner
  # bounds coincide with the endpoint values (i.e. a hard Bessel layer)
  layers(bm, 0, 2, type = "bessel")
  x <- bm$W_t[match(bm$layers$t.l, bm$t)]
  y <- bm$W_tm[match(bm$layers$t.u, bm$t)]
  hard <- bm$layers$Lu == pmin(x, y) & bm$layers$Ud == pmax(x, y)
  expect_true(any(hard))
  expect_equal(bm$layers$type[hard], rep("bessel", sum(hard)))
  expect_equal(bm$layers$type[!hard], rep("intersection", sum(!hard)))
  expect_valid_bm(bm)
})

test_that("layers() localised between existing points creates bridge layers", {
  withr::local_seed(15)
  bm <- create_bm(c(0, 1, 2), c(0, 0.5, 0))
  layers(bm, 0, 2, type = "localised")
  expect_true(all(bm$layers$type %in% c("localised-bb", "bessel-bb")))
  expect_true(all(c(0, 1, 2) %in% bm$t))
  expect_equal(unname(bm$W_t[match(c(0, 1, 2), bm$t)]), c(0, 0.5, 0))
  # Bridge layers tile [0, 2]
  lyr <- bm$layers[order(bm$layers$t.l), ]
  expect_equal(min(lyr$t.l), 0)
  expect_equal(max(lyr$t.u), 2)
  expect_equal(utils::head(lyr$t.u, -1), utils::tail(lyr$t.l, -1))
  expect_valid_bm(bm)
})

test_that("layers() localised beyond the end of the path extends it", {
  withr::local_seed(16)
  bm <- create_bm()
  layers(bm, 0, 2, type = "localised")
  expect_equal(max(bm$t), 2)
  expect_true(any(bm$layers$type == "localised"))
  expect_valid_bm(bm)
})

test_that("layers() localised then sim inside bridge layers keeps bounds", {
  withr::local_seed(17)
  bm <- create_bm(c(0, 1, 2), c(0, 0.5, 0))
  layers(bm, 0, 2, type = "localised")
  sim(bm, c(0.3, 1.5))
  expect_true(all(c(0.3, 1.5) %in% bm$t))
  expect_true(all(c(0.3, 1.5) %in% bm$bb.local$t))
  expect_valid_bm(bm)
})

test_that("layers() localised with a larger mult uses fewer, wider layers", {
  withr::local_seed(18)
  bm1 <- create_bm(c(0, 1), c(0, 0))
  layers(bm1, 0, 1, type = "localised", mult = 0.25)
  withr::local_seed(18)
  bm5 <- create_bm(c(0, 1), c(0, 0))
  layers(bm5, 0, 1, type = "localised", mult = 5)
  expect_lt(nrow(bm5$layers), nrow(bm1$layers))
  expect_valid_bm(bm1)
  expect_valid_bm(bm5)
})

test_that("layers() validates its arguments", {
  bm <- local_bm_1d(1)
  expect_error(layers(bm, 2, 1), "s must be less than t")
  expect_error(layers(bm, 1, 1), "s must be less than t")
  expect_error(layers(bm, "a", 1), "must be real scalars")
  expect_error(layers(bm, 0, 1, type = "foo"), "layer type must be one of")
  expect_error(layers(bm, 0, 1, mult = 0), "mult must be strictly positive")
  expect_error(layers(bm, 0, 1, prefer = "localised"),
               "must be one of 'bessel' or 'intersection'")
  expect_error(layers(bm, 0, 1, label = "end"), "clashes with an internal label")
  expect_equal(nrow(bm$layers), 0)
})

test_that("layers() labels both endpoints", {
  withr::local_seed(19)
  bm <- local_bm_1d()
  layers(bm, 0.5, 2.5, label = c("left", "right"))
  expect_equal(bm$labels$left, 0.5)
  expect_equal(bm$labels$right, 2.5)
})

test_that("Brownian bridge marginal at 0.5 is N(0.25, 1/4) for every layer type", {
  # Bridge from W_0 = 0 to W_1 = 0.5: W_0.5 ~ N(0.25, 0.25) whatever layer is
  # used. The endpoints differ because intersection layers on intervals with
  # equal endpoint values can loop forever (see test-known-bugs.R).
  skip_on_cran()
  n <- 150
  for (type in c("bessel", "localised", "intersection")) {
    # Fixed seed chosen to avoid known bug 1 for intersection layers
    withr::local_seed(20)
    x <- vapply(seq_len(n), function(i) {
      bm <- create_bm(c(0, 1), c(0, 0.5))
      layers(bm, 0, 1, type = type)
      sim(bm, 0.5)
      unname(bm$W_t[bm$t == 0.5])
    }, numeric(1))
    expect_lt(abs(mean(x) - 0.25), 4 * 0.5 / sqrt(n), label = paste(type, "mean"))
    expect_gt(stats::ks.test(x, "pnorm", mean = 0.25, sd = 0.5)$p.value, 0.001,
              label = paste(type, "KS p-value"))
  }
})
