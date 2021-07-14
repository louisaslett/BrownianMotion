sim.bb_ <- function(bm, s_idx, q, t_idx, label) {
  n <- length(q)

  W_t <- rep(0, n)
  W_t[1] <- bm$W_t[s_idx]
  W_t[n] <- bm$W_t[t_idx]

  for(i in 2:(n-1)) {
    W_t[i] <- rnorm(1,
                    mean = W_t[i-1] + (q[i]-q[i-1])/(q[n]-q[i-1])*(W_t[n]-W_t[i-1]),
                    sd = sqrt((q[n]-q[i])*(q[i]-q[i-1])/(q[n]-q[i-1])))
  }
  bm$t <- c(bm$t[1:s_idx],
            head(tail(q, -1), -1),
            bm$t[t_idx:length(bm$t)])
  bm$W_t <- c(bm$W_t[1:s_idx],
              head(tail(W_t, -1), -1),
              bm$W_t[t_idx:length(bm$W_t)])

  add.labels_(bm, "user", q[2:(n-1)])
  add.labels_(bm, label, q[2:(n-1)])

  bm
}
