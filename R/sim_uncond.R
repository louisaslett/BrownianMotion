sim.uncond_ <- function(bm, t, label) {
  t_delta <- c(t[1] - tail(bm$t, 1), diff(t))
  bm$t <- c(bm$t, unname(t))

  W_delta <- rnorm(length(t_delta), sd = sqrt(t_delta))
  W_delta[1] <- W_delta[1] + tail(bm$W_t, 1)
  bm$W_t <- c(bm$W_t, cumsum(W_delta))
  bm$W_tm <- c(bm$W_tm, cumsum(W_delta))

  add.labels_(bm, "user", t)
  add.labels_(bm, label, t)
  bm$labels[["end"]] <- max(t)
  bm
}
