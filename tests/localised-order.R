library(BrownianMotion)

# Recover the four sloping auxiliary bounds at both ends of every bridge layer.
# These describe the same physical bounds regardless of the auxiliary origin.
bridge_bounds <- function(bm) {
  master <- bm$layers[grepl("-bb$", bm$layers$type), ]
  aux <- bm$bb.local$layers
  idx <- match(master$t.l, aux$t.l)
  stopifnot(!anyNA(idx), all(master$t.u == aux$t.u[idx]))
  bounds <- as.matrix(aux[idx, c("Ld", "Lu", "Ud", "Uu")])
  list(t = master$t.l,
       left = sweep(bounds, 1, master$Ld, "+"),
       right = sweep(bounds, 1, master$Uu, "+"))
}

check_offsets <- function(bm) {
  master <- bm$layers[grepl("-bb$", bm$layers$type), ]
  aux <- bm$bb.local
  stopifnot(
    length(aux$t) == length(aux$W_t),
    length(aux$t) == length(aux$W_tm),
    !anyNA(aux$W_tm),
    isTRUE(all.equal(master$Ld,
                    bm$W_t[match(master$t.l, bm$t)] -
                      aux$W_t[match(master$t.l, aux$t)],
                    tolerance = 1e-12, check.attributes = FALSE)),
    isTRUE(all.equal(master$Uu,
                    bm$W_tm[match(master$t.u, bm$t)] -
                      aux$W_t[match(master$t.u, aux$t)],
                    tolerance = 1e-12, check.attributes = FALSE))
  )
}

check_preserved <- function(bm, before, old.t, old.values) {
  after <- bridge_bounds(bm)
  idx <- match(before$t, after$t)
  stopifnot(
    !anyNA(idx),
    isTRUE(all.equal(before$left, after$left[idx, , drop = FALSE],
                    tolerance = 1e-12)),
    isTRUE(all.equal(before$right, after$right[idx, , drop = FALSE],
                    tolerance = 1e-12)),
    identical(old.values, bm$W_t[match(old.t, bm$t)])
  )
  check_offsets(bm)
}

# Reverse construction, including bridge layers already split by simulation.
for (prefer in c("bessel", "intersection")) {
  for (split in c(FALSE, TRUE)) {
    set.seed(1)
    bm <- create_bm(c(0, 1, 2), c(0, 0, 0), prefer = prefer)
    layers(bm, 1, 2, type = "localised")
    if (split) sim(bm, 1.2)
    before <- bridge_bounds(bm)
    old.t <- bm$t
    old.values <- bm$W_t
    layers(bm, 0, 1, type = "localised")
    check_preserved(bm, before, old.t, old.values)

    sim(bm, c(0.25, 0.75, 1.75))
    refine(bm, 0, 2)
    check_offsets(bm)
  }
}

# Fill a gap between existing bridge intervals. Ordinary layers to the right
# must not be adjusted, and the bridge interval to the left must be preserved.
set.seed(1)
bm <- create_bm(0:4, rep(0, 5))
layers(bm, 0, 1, type = "localised")
layers(bm, 2, 3, type = "localised")
layers(bm, 3, 4, type = "intersection")
ordinary <- bm$layers[bm$layers$t.l >= 3, ]
before <- bridge_bounds(bm)
old.t <- bm$t
old.values <- bm$W_t
layers(bm, 1, 2, type = "localised")
check_preserved(bm, before, old.t, old.values)
stopifnot(identical(ordinary, bm$layers[bm$layers$t.l >= 3, ]))
sim(bm, c(1.25, 2.75))
refine(bm, 0, 3)
check_offsets(bm)
