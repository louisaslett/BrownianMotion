test_that("plot() of a 1D path builds without warnings", {
  withr::local_seed(101)
  bm <- create_bm(c(-1, 0, 1), c(100, 101, 99))
  expect_plot_ok(plot(bm))
  sim(bm, 2:4)
  expect_plot_ok(plot(bm))
})

test_that("plot() of a 1D path with every layer type builds", {
  withr::local_seed(102)
  bm <- create_bm(prefer = "intersection")
  layers(bm, 1, 3, type = "bessel")
  layers(bm, 3, 4, type = "intersection")
  layers(bm, 4, 5, type = "localised")
  first_passage(bm, delta = 1)
  expect_plot_ok(plot(bm))
  expect_plot_ok(plot(bm, hide.user = TRUE))
  expect_plot_ok(plot(bm, t.lim = c(1, 4.5)))
})

test_that("plot() of a 1D path with bridge layers of each kind builds", {
  withr::local_seed(103)
  bm <- create_bm(c(0, 1, 2, 3), c(0, 0, 0, 0), prefer = "intersection")
  layers(bm, 0, 1, type = "localised")
  layers(bm, 2, 3, type = "localised", prefer = "bessel")
  expect_true(any(grepl("-bb$", bm$layers$type)))
  expect_plot_ok(plot(bm))
  expect_plot_ok(plot(bm, t.lim = c(0, 1), hide.user = TRUE))
})

test_that("plot() of a 1D path with jumps builds", {
  withr::local_seed(104)
  bm <- create_bm()
  sim(bm, 0:5)
  add_segment(bm, l = 1, r = 3, delta = 1)
  expect_plot_ok(plot(bm))
})

test_that("points() adds a 1D path to a base graphics plot", {
  withr::local_seed(105)
  bm <- local_bm_1d()
  withr::local_pdf(NULL)
  plot(bm$t, bm$W_t, type = "n")
  expect_no_error(points(bm))
})

test_that("plot() of an N-d path with layers builds in each layers.2d mode", {
  withr::local_seed(106)
  bm <- local_bm_nd(3)
  sim(bm, c(0.5, 5))
  layers(bm, 0, 10)
  expect_plot_ok(plot(bm))
  expect_plot_ok(plot(bm, layers.2d = "convex"))
  expect_plot_ok(plot(bm, layers.2d = "detailed"))
  expect_plot_ok(plot(bm, t.lim = c(0, 0.5)))
  expect_plot_ok(plot(bm, dims = 2))
  expect_plot_ok(plot(bm, dims = c(3, 2), pairs = FALSE))
})

test_that("plot() of an N-d path builds with layers.2d = 'union'", {
  skip_if_not_installed("sf")
  withr::local_seed(107)
  bm <- local_bm_nd(2)
  sim(bm, 1)
  layers(bm, 0, 2)
  expect_plot_ok(plot(bm, layers.2d = "union"))
  expect_plot_ok(plot(bm, layers.2d = "union", dims = c(2, 1), pairs = FALSE))
})

test_that("plot() of bm$Z builds", {
  withr::local_seed(108)
  bm <- local_bm_nd(3)
  sim(bm, c(0.5, 5))
  layers(bm, 0, 10)
  expect_plot_ok(plot(bm$Z))
  expect_plot_ok(plot(bm$Z, t.lim = c(0, 0.5)))
  expect_plot_ok(plot(bm$Z, layers.2d = "detailed"))
  expect_plot_ok(plot(bm$Z, dims = 1))
})

test_that("plot() of an N-d path warns when layers do not cover the time range", {
  withr::local_seed(109)
  bm <- local_bm_nd(2)
  sim(bm, 1)
  layers(bm, 0, 1)
  sim(bm, 2)
  expect_warning(plot(bm), "no layer information")
})

test_that("plot() of an N-d path rejects unknown layers.2d settings", {
  withr::local_seed(110)
  bm <- local_bm_nd(2)
  layers(bm, 0, 1)
  expect_error(build_plot(plot(bm, layers.2d = "foo")), "Unrecognised setting")
})

test_that("plot() of a concatenated N-d path builds", {
  withr::local_seed(111)
  bm <- local_bm_nd(3)
  sim(bm, 1:2)
  layers(bm, 0, 2)
  expect_plot_ok(plot(concat_bm(bm, bm, t0 = 5)))
})
