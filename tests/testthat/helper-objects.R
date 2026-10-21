# Small reusable fixtures. Each call creates a fresh object, so tests never
# share mutable state. Callers set the seed (withr::local_seed()) first.

local_bm_1d <- function(times = 1:3) {
  bm <- create_bm()
  sim(bm, times)
  bm
}

local_bm_nd <- function(dim = 2, rho = 0.5) {
  cov <- matrix(rho, dim, dim)
  diag(cov) <- 1
  create_bm(dim = dim, cov = cov)
}
