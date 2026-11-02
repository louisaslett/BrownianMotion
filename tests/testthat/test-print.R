test_that("print() shows the skeleton, layers, jumps and labels", {
  withr::local_seed(81)
  bm <- create_bm()
  expect_snapshot(print(bm))

  sim(bm, c(1, peak = 2, 3, 4))
  layers(bm, 0, 2, type = "bessel")
  layers(bm, 2, 3, type = "intersection")
  add_segment(bm, l = 3, delta = 1)
  expect_snapshot(print(bm))
  expect_snapshot(print(bm, digits = 2))
})

test_that("print() shows the first and last points of long paths", {
  withr::local_seed(82)
  bm <- local_bm_1d(1:30)
  expect_snapshot(print(bm))
  expect_snapshot(print(bm, n = 3))
  expect_length(format(bm, n = Inf), 30 + 1 + 4)
})

test_that("print() shows localised and bridge layers with their bounds", {
  withr::local_seed(83)
  bm <- local_bm_1d(1:2)
  layers(bm, 0, 1, type = "localised")
  layers(bm, 2, 3, type = "localised")
  expect_snapshot(print(bm, n = Inf))
})

test_that("print() returns its argument invisibly", {
  withr::local_seed(84)
  bm <- local_bm_1d()
  bm2 <- local_bm_nd()
  for(x in list(bm, bm2, bm2$Z)) {
    expect_output(res <- withVisible(print(x)), "Brownian motion")
    expect_false(res$visible)
    expect_identical(res$value, x)
    expect_type(format(x), "character")
  }
  expect_output(print(bm), "Brownian motion")
})

test_that("print() shows multi-dimensional paths and their transformed view", {
  withr::local_seed(85)
  bm <- local_bm_nd(dim = 2, rho = 0.5)
  sim(bm, c(1, last = 2))
  layers(bm, 0, 2)
  expect_snapshot(print(bm))
  expect_snapshot(print(bm$Z))

  bm3 <- create_bm(dim = 3, cov = c(1, 2, 3))
  sim(bm3, 1)
  add_segment(bm3, l = 1, delta = 1)
  expect_snapshot(print(bm3))
})

test_that("print() fits the output to the width", {
  withr::local_seed(86)
  bm <- create_bm(dim = 4)
  sim(bm, 1:2)
  out <- format(bm, width = 40)
  expect_true(all(cli::ansi_nchar(out) <= 40))
  expect_match(out, "components? not shown", all = FALSE)

  bm1 <- local_bm_1d()
  add_labels(bm1, t = 1:3, label = strrep("x", 40))
  out <- format(bm1, width = 50)
  expect_true(all(cli::ansi_nchar(out) <= 50))
  out <- format(bm1, width = 30)
  expect_match(out, "Labels not shown", all = FALSE, fixed = TRUE)
})

test_that("print() checks its arguments", {
  bm <- create_bm()
  expect_error(print(bm, n = 0), "n must be")
  expect_error(print(bm, n = "a"), "n must be")
  expect_error(print(bm, digits = 0), "digits must be")
  expect_error(print(summary(bm), digits = 1.5), "digits must be")
})

test_that("summary() of a 1-dimensional path", {
  withr::local_seed(87)
  bm <- create_bm()
  sim(bm, c(1, 2, peak = 3, 4))
  layers(bm, 0, 2, type = "bessel")
  layers(bm, 3, 4, type = "intersection")
  first_passage(bm, delta = 1)
  add_segment(bm, l = 2, r = 3, delta = 1)

  s <- summary(bm)
  expect_s3_class(s, "summary.BrownianMotion")
  expect_equal(s$n.points, length(bm$t))
  expect_equal(s$t.range, range(bm$t))
  expect_equal(sum(s$origin), length(bm$t))
  expect_equal(s$origin[["fixed"]], 1)
  expect_equal(s$origin[["passage"]], 1)
  expect_equal(s$path["min", "W_t"], min(bm$W_t, bm$W_tm))
  expect_equal(s$path["max", "W_t"], max(bm$W_t, bm$W_tm))
  expect_equal(s$jumps$t, c(2, 3))
  expect_equal(s$jumps$size, c(1, -1))
  expect_equal(s$n.layers, nrow(bm$layers))
  expect_equal(sum(s$layers$intervals), nrow(bm$layers))
  expect_equal(s$gaps, data.frame(from = 2, to = 3))
  expect_equal(s$coverage, 1 - 1 / diff(range(bm$t)))
  expect_lte(s$bounds[["L"]], min(bm$W_t[bm$t <= 2 | bm$t >= 3]))
  expect_gte(s$bounds[["U"]], max(bm$W_t[bm$t <= 2 | bm$t >= 3]))
  expect_equal(s$labels, c(peak = 1L))
  expect_equal(s$passages, 1)
  expect_equal(s$settings, list(refine = TRUE, mult = 1, prefer = "bessel"))

  expect_snapshot(print(s))
  expect_output(res <- withVisible(print(s)), "summary")
  expect_false(res$visible)
  expect_identical(res$value, s)
})

test_that("summary() of a path with no layers or a single point", {
  bm <- create_bm()
  s <- summary(bm)
  expect_null(s$spacing)
  expect_null(s$bounds)
  expect_true(is.na(s$coverage))
  expect_equal(s$n.layers, 0)
  expect_snapshot(print(s))
})

test_that("summary() bounds contain localised bridge layers", {
  withr::local_seed(88)
  bm <- local_bm_1d(1:3)
  layers(bm, 0, 3, type = "localised")
  expect_true(all(grepl("-bb$", bm$layers$type)))
  s <- summary(bm)
  expect_equal(s$coverage, 1)
  expect_lte(s$bounds[["L"]], min(bm$W_t[bm$t <= 3]))
  expect_gte(s$bounds[["U"]], max(bm$W_t[bm$t <= 3]))

  # Every bridge layer bounds the path at its ends
  lyrs <- layer_bounds_(bm)
  for(i in seq_len(nrow(lyrs))) {
    x <- interval_values(bm, lyrs$t.l[i], lyrs$t.u[i])
    expect_true(all(x >= lyrs$L[i] - 1e-8 & x <= lyrs$U[i] + 1e-8))
  }
})

test_that("summary() of a multi-dimensional path and its transformed view", {
  withr::local_seed(89)
  bm <- create_bm(dim = 3, cov = matrix(c(1, 0.5, 0.2, 0.5, 2, 0.3, 0.2, 0.3, 1), 3))
  sim(bm, c(1, 2, mid = 3, 4))
  layers(bm, 0, 3, type = "intersection")

  s <- summary(bm)
  expect_s3_class(s, "summary.BrownianMotionNd")
  expect_equal(s$view, "W")
  expect_equal(s$cov, bm$cov)
  expect_equal(s$components$start, unname(bm[1, "W_t"][1, ]))
  expect_equal(s$components$max, unname(apply(bm[, "W_t"], 2, max)))
  expect_equal(s$layers["intersection", ], c("Z[1]" = 3L, "Z[2]" = 3L, "Z[3]" = 3L))
  expect_equal(s$coverage, 3 / 4)
  expect_equal(s$labels, c(mid = 1L))

  # The layer bounds are the outer box of the layers given by indexing
  lyrs <- bm[, "layers"]
  expect_equal(s$components$L, apply(lyrs$outer.L, 2, min))
  expect_equal(s$components$U, apply(lyrs$outer.U, 2, max))
  expect_snapshot(print(s))

  z <- summary(bm$Z)
  expect_equal(z$view, "Z")
  expect_equal(z$cov, diag(3))
  expect_equal(z$components$end, unname(bm$Z[5, "Z_t"][1, ]))
  expect_snapshot(print(z))
})

test_that("summary() of a multi-dimensional path bounds localised bridge layers", {
  withr::local_seed(90)
  bm <- local_bm_nd(dim = 2, rho = 0.5)
  sim(bm, 1:2)
  layers(bm, 0, 2, type = "localised")
  s <- summary(bm)
  W <- bm[, "W_t"]
  expect_true(all(s$components$L <= apply(W, 2, min) + 1e-8))
  expect_true(all(s$components$U >= apply(W, 2, max) - 1e-8))
})
