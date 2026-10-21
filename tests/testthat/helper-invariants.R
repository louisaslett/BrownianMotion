# Shared structural checks for BrownianMotion objects.
#
# The checks are written as expectations so a failure points at the broken
# invariant. They are deliberately conservative: they only assert properties
# that must hold for every valid skeleton.

bb_types <- c("localised-bb", "bessel-bb", "intersection-bb")

# Path values used to bound the interval [t.l, t.u]: W_t at every skeleton time
# in [t.l, t.u) and the left limit W_tm at t.u.
interval_values <- function(bm, t.l, t.u) {
  inside <- bm$t >= t.l & bm$t < t.u
  c(unname(bm$W_t[inside]), unname(bm$W_tm[match(t.u, bm$t)]))
}

# Basic shape of the skeleton: sorted unique times, equal lengths, no NAs.
expect_valid_skeleton <- function(bm, tol = 1e-8) {
  expect_true(is.numeric(bm$t))
  expect_false(anyNA(bm$t))
  expect_false(is.unsorted(bm$t, strictly = TRUE))
  expect_length(bm$W_t, length(bm$t))
  expect_length(bm$W_tm, length(bm$t))
  expect_false(anyNA(bm$W_t))
  expect_false(anyNA(bm$W_tm))
}

# Non-bridge layers must tile skeleton intervals and bound the path.
expect_layers_bound_path <- function(bm, tol = 1e-8) {
  lyrs <- bm$layers[!(bm$layers$type %in% bb_types), ]
  if (nrow(lyrs) == 0) {
    return(invisible(bm))
  }
  expect_true(all(lyrs$t.l < lyrs$t.u))
  expect_true(all(lyrs$t.l %in% bm$t))
  expect_true(all(lyrs$t.u %in% bm$t))
  expect_true(all(lyrs$Ld <= lyrs$Lu + tol))
  expect_true(all(lyrs$Ud <= lyrs$Uu + tol))
  for (i in seq_len(nrow(lyrs))) {
    x <- interval_values(bm, lyrs$t.l[i], lyrs$t.u[i])
    expect_true(all(x >= lyrs$Ld[i] - tol & x <= lyrs$Uu[i] + tol),
                label = paste0("path inside ", lyrs$type[i], " layer [",
                               lyrs$t.l[i], ", ", lyrs$t.u[i], "]"))
  }
  invisible(bm)
}

# Bridge (-bb) layers store offsets B - W_aux at each end of the interval,
# relative to the auxiliary path in bm$bb.local.
expect_bridge_offsets <- function(bm, tol = 1e-10) {
  master <- bm$layers[bm$layers$type %in% bb_types, ]
  if (nrow(master) == 0) {
    return(invisible(bm))
  }
  aux <- bm$bb.local
  idx <- match(paste(master$t.l, master$t.u),
               paste(aux$layers$t.l, aux$layers$t.u))
  expect_false(anyNA(idx), label = "every bridge layer has an auxiliary layer")
  expect_equal(paste0(aux$layers$type[idx], "-bb"), master$type)
  expect_true(all(is.na(master$Lu)) && all(is.na(master$Ud)))
  expect_equal(master$Ld,
               bm$W_t[match(master$t.l, bm$t)] - aux$W_t[match(master$t.l, aux$t)],
               tolerance = tol, ignore_attr = TRUE)
  # Right ends use left limits on both paths, which differ from W_t at jumps
  expect_equal(master$Uu,
               bm$W_tm[match(master$t.u, bm$t)] - aux$W_tm[match(master$t.u, aux$t)],
               tolerance = tol, ignore_attr = TRUE)
  invisible(bm)
}

# Layers (bridge or not) should not overlap one another.
expect_layers_disjoint <- function(bm) {
  lyrs <- bm$layers[order(bm$layers$t.l), ]
  if (nrow(lyrs) > 1) {
    expect_true(all(utils::head(lyrs$t.u, -1) <= utils::tail(lyrs$t.l, -1)))
  }
  invisible(bm)
}

# Every labelled time is a skeleton time; start/end track the range.
expect_labels_consistent <- function(bm) {
  expect_equal(bm$labels$start, min(bm$t))
  expect_equal(bm$labels$end, max(bm$t))
  expect_true(all(unlist(bm$labels) %in% bm$t))
  invisible(bm)
}

# All of the above, for a 1D object and its auxiliary bridge path.
expect_valid_bm <- function(bm, check_offsets = TRUE, tol = 1e-8) {
  expect_s3_class(bm, "BrownianMotion")
  expect_valid_skeleton(bm)
  expect_layers_bound_path(bm, tol = tol)
  expect_layers_disjoint(bm)
  expect_labels_consistent(bm)
  if (!is.null(bm$bb.local)) {
    expect_valid_skeleton(bm$bb.local)
    expect_layers_bound_path(bm$bb.local, tol = tol)
    if (check_offsets) expect_bridge_offsets(bm)
  }
  invisible(bm)
}

# Multi-dimensional objects: every component path is valid and all share the
# same skeleton times.
expect_valid_bm_nd <- function(bm, ...) {
  expect_s3_class(bm, "BrownianMotionNd")
  expect_length(bm$Z.bm, bm$dim)
  for (d in seq_len(bm$dim)) {
    expect_valid_bm(bm$Z.bm[[d]], ...)
    expect_equal(bm$Z.bm[[d]]$t, bm$Z.bm[[1]]$t)
  }
  invisible(bm)
}
