# Helpers for checking localised bridge (-bb) layers, migrated from the
# original tests/localised-order.R regression script.

# Recover the four sloping auxiliary bounds at both ends of every bridge layer.
# These describe the same physical bounds regardless of the auxiliary origin.
bridge_bounds <- function(bm) {
  master <- bm$layers[grepl("-bb$", bm$layers$type), ]
  aux <- bm$bb.local$layers
  idx <- match(master$t.l, aux$t.l)
  expect_false(anyNA(idx), label = "every bridge layer has an auxiliary layer")
  expect_equal(master$t.u, aux$t.u[idx])
  bounds <- as.matrix(aux[idx, c("Ld", "Lu", "Ud", "Uu")])
  list(t = master$t.l,
       left = sweep(bounds, 1, master$Ld, "+"),
       right = sweep(bounds, 1, master$Uu, "+"))
}

# The stored offsets must equal B - W_aux at both ends of each bridge layer,
# and the auxiliary path must have consistent vectors.
check_offsets <- function(bm) {
  master <- bm$layers[grepl("-bb$", bm$layers$type), ]
  aux <- bm$bb.local
  expect_length(aux$W_t, length(aux$t))
  expect_length(aux$W_tm, length(aux$t))
  expect_false(anyNA(aux$W_tm))
  expect_equal(master$Ld,
               bm$W_t[match(master$t.l, bm$t)] -
                 aux$W_t[match(master$t.l, aux$t)],
               tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(master$Uu,
               bm$W_tm[match(master$t.u, bm$t)] -
                 aux$W_t[match(master$t.u, aux$t)],
               tolerance = 1e-12, ignore_attr = TRUE)
  invisible(bm)
}

# Physical bridge bounds recorded in `before` are unchanged, as are the path
# values at the previously simulated times.
check_preserved <- function(bm, before, old.t, old.values) {
  after <- bridge_bounds(bm)
  idx <- match(before$t, after$t)
  expect_false(anyNA(idx), label = "earlier bridge layers still present")
  expect_equal(after$left[idx, , drop = FALSE], before$left, tolerance = 1e-12)
  expect_equal(after$right[idx, , drop = FALSE], before$right, tolerance = 1e-12)
  expect_identical(bm$W_t[match(old.t, bm$t)], old.values)
  check_offsets(bm)
}
