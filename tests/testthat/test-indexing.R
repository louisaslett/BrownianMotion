local_bm_indexed <- function() {
  bm <- create_bm(c(-1, 0, 1), c(100, 101, 99))
  add_points(bm, c(2, 3), c(98, 97))
  layers(bm, 1, 3)
  add_labels(bm, t = 2, label = "two")
  bm
}

test_that("[ with no index returns every field", {
  withr::local_seed(91)
  bm <- local_bm_indexed()
  res <- bm[]
  expect_named(res, c("t", "W_t", "W_tm", "layers", "labels"))
  expect_equal(res$t, bm$t)
  expect_equal(res$W_t, unname(bm$W_t))
})

test_that("[ extracts t, W_t, W_tm, layers and labels by name", {
  withr::local_seed(92)
  bm <- local_bm_indexed()
  expect_equal(bm[, "t"], c(-1, 0, 1, 2, 3))
  expect_equal(bm[, "W_t"], c(100, 101, 99, 98, 97))
  expect_equal(bm[, "W_tm"], c(100, 101, 99, 98, 97))
  expect_equal(bm[, "layers"], bm$layers)
  expect_equal(bm[, "labels"]$two, 2)
  expect_equal(bm[, "labels"]$user, c(2, 3))
})

test_that("[ by time, by index and by unordered time", {
  withr::local_seed(93)
  bm <- local_bm_indexed()
  expect_equal(bm[t = 1, "W_t"], 99)
  expect_equal(bm[1, "W_t"], 100)
  expect_equal(bm[c(2, 1), "W_t"], c(101, 100))
  expect_equal(bm[t = c(2, 0), "W_t"], c(98, 101))
  expect_equal(bm[t = 0, "t"], 0)
})

test_that("[ by time not in the skeleton returns NA and covering layers", {
  withr::local_seed(94)
  bm <- local_bm_indexed()
  expect_equal(bm[t = 1.5, "W_t"], NA)
  lyr <- bm[t = 1.5, "layers"]
  expect_equal(nrow(lyr), 1)
  expect_true(lyr$t.l <= 1.5 && lyr$t.u > 1.5)
  expect_equal(bm[t = 0.5, "t"], 0.5)
  expect_equal(nrow(bm[t = 0.5, "layers"]), 0)
  expect_equal(bm[t = 0.5, "labels"], list())
  # Mixture of skeleton and non-skeleton times keeps the requested order
  expect_equal(bm[t = c(2, 1.5), "W_t"], c(98, NA))
})

test_that("[ by time range t = a:b includes all skeleton points in between", {
  withr::local_seed(95)
  bm <- local_bm_indexed()
  sim(bm, c(1.25, 1.75))
  res <- bm[t = 1:2]
  expect_equal(res$t, c(1, 1.25, 1.75, 2))
  expect_equal(res$W_t, unname(bm$W_t[match(res$t, bm$t)]))
  # Layers covering any of the returned times (including [2, 3) at t = 2)
  expect_true(all(res$layers$t.l <= 2 & res$layers$t.u > 1))
  expect_true(any(res$layers$t.l == 1))
})

test_that("[ by label strings", {
  withr::local_seed(96)
  bm <- local_bm_indexed()
  expect_equal(bm["two", "t"], 2)
  expect_equal(bm["two", "W_t"], 98)
  expect_equal(bm["user", "t"], c(2, 3))
  expect_equal(bm["end", "t"], 3)
  expect_equal(bm["start", "W_t"], 100)
})

test_that("[ rejects non-integer and out of bounds indices", {
  withr::local_seed(97)
  bm <- local_bm_indexed()
  expect_error(bm[1.5, "W_t"], "must be an integer index")
  expect_error(bm[10, "W_t"], "index out of bounds")
})

test_that("[ on N-d objects returns matrices for W_t", {
  withr::local_seed(98)
  bm <- local_bm_nd(3)
  sim(bm, c(0.5, 5))
  layers(bm, 0, 10)
  W <- bm[, "W_t"]
  expect_true(is.matrix(W))
  expect_equal(ncol(W), 3)
  expect_equal(nrow(W), length(bm[, "t"]))
  expect_equal(bm[t = 0.5, "W_t"], W[bm[, "t"] == 0.5, , drop = FALSE])
  expect_equal(bm[4, "W_t"], W[4, , drop = FALSE])
  lyr <- bm[, "layers"]
  expect_true(all(c("t.l", "t.u", "outer.L", "outer.U", "inner.cube") %in% names(lyr)))
  expect_equal(ncol(lyr$outer.L), 3)
  expect_true(all(lyr$outer.L <= lyr$outer.U))
  expect_equal(bm[t = 0.25, "W_t"], matrix(NA, 1, 3))
})

test_that("bm$Z indexing returns the transformed path", {
  withr::local_seed(99)
  bm <- local_bm_nd(3)
  sim(bm, c(0.5, 5))
  layers(bm, 0, 10)
  Z <- bm$Z[, "Z_t"]
  W <- bm[, "W_t"]
  expect_equal(dim(Z), dim(W))
  expect_equal(bm$Z[, "t"], bm[, "t"])
  expect_equal(bm$Z[t = 0.5, "Z_t"], Z[bm[, "t"] == 0.5, , drop = FALSE])
  expect_equal(bm$Z[4, "Z_t"], Z[4, , drop = FALSE])
  expect_equal(unname(bm$Z[, "Z_t"][, 2]), unname(bm$Z.bm[[2]]$W_t))
  expect_s3_class(bm$Z[, "layers"], "tbl_df")
  expect_named(bm$Z[], c("t", "Z_t", "Z_tm", "layers", "labels"))
})

test_that("[ by unsorted indices returns every variable in request order", {
  # Index results were permuted with order(x) rather than its inverse
  # order(order(x)), which only coincide for two or fewer indices.
  withr::local_seed(41)
  idx <- c(4, 2, 5, 1, 3, 2)

  bm <- create_bm(c(0, 1), c(0, 0))
  sim(bm, c(2, 3, 4))
  add_segment(bm, l = 2, delta = 1)          # a jump so W_t and W_tm differ
  expect_equal(bm[idx, "t"], bm$t[idx])
  expect_equal(bm[idx, "W_t"], bm$W_t[idx])
  expect_equal(bm[idx, "W_tm"], bm$W_tm[idx])
  full <- bm[idx]
  expect_equal(full$t, bm$t[idx])
  expect_equal(full$W_t, bm$W_t[idx])
  expect_equal(full$W_tm, bm$W_tm[idx])
  # Same answer as asking by (unsorted) time
  expect_equal(bm[unique(idx), "W_t"], bm[t = bm$t[unique(idx)], "W_t"])

  nd <- create_bm(dim = 2, cov = matrix(c(1, 0.5, 0.5, 2), 2))
  sim(nd, 1:4)
  all.W <- nd[, "W_t"]
  expect_equal(nd[idx, "t"], (0:4)[idx])
  expect_equal(nd[idx, "W_t"], all.W[idx, , drop = FALSE])
  expect_equal(nd[idx]$W_t, all.W[idx, , drop = FALSE])
  expect_equal(nd[idx, "W_tm"], nd[, "W_tm"][idx, , drop = FALSE])
  all.Z <- nd$Z[, "Z_t"]
  expect_equal(nd$Z[idx, "Z_t"], all.Z[idx, , drop = FALSE])
  expect_equal(nd$Z[idx]$Z_t, all.Z[idx, , drop = FALSE])
})

test_that("[ evaluates var, so it can be held in a variable", {
  # var used to be taken from the call unevaluated, so this returned NULL
  withr::local_seed(42)
  bm <- create_bm()
  sim(bm, 1:3)
  what <- "W_t"
  expect_equal(bm[2, what], bm[2, "W_t"])
  expect_equal(bm[t = 2, what], bm[t = 2, "W_t"])
  get_var <- function(b, v) b[, v]
  expect_equal(get_var(bm, "W_tm"), bm[, "W_tm"])
  expect_equal(get_var(bm, "labels"), bm[, "labels"])

  nd <- create_bm(dim = 2)
  sim(nd, 1:2)
  expect_equal(nd[t = 1, what], nd[t = 1, "W_t"])
  zwhat <- "Z_t"
  expect_equal(nd$Z[2, zwhat], nd$Z[2, "Z_t"])
})

test_that("[ rejects indices that do not match bm[i, var] or bm[t = times, var]", {
  # These used to be silently ignored or return NULL
  withr::local_seed(43)
  bm <- create_bm()
  sim(bm, 1:3)
  expect_error(bm[1, "W_t", t = 2], "too many indices")
  expect_error(bm[1, "W_t", 2], "too many indices")
  expect_error(bm[1, t = 2], "times must be given by name as the first index")
  expect_error(bm["W_t", t = 2], "times must be given by name as the first index")
  expect_error(bm[i = 1, "W_t"], "other indices must not be named")
  expect_error(bm[t = 2, 1], "var must be one of")
  expect_error(bm[2, "foo"], 'var must be one of "t", "W_t", "W_tm", "layers", "labels"')
  expect_error(bm[2, c("t", "W_t")], "var must be one of")

  nd <- create_bm(dim = 2)
  sim(nd, 1:2)
  expect_error(nd[1, "W_t", t = 1], "too many indices")
  expect_error(nd[2, "Z_t"], "var must be one of")
  expect_error(nd$Z[1, "W_t"], 'var must be one of "t", "Z_t", "Z_tm", "layers", "labels"')
  expect_error(nd$Z[1, t = 1], "times must be given by name as the first index")

  # Valid forms with no var still return everything
  expect_named(bm[2], c("t", "W_t", "W_tm", "layers", "labels"))
  expect_named(bm[t = 2], c("t", "W_t", "W_tm", "layers", "labels"))
  expect_named(bm[2, ], c("t", "W_t", "W_tm", "layers", "labels"))
})
