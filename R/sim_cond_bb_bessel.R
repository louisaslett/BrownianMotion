#' Simulate Brownian motion conditional on Brownian bridge Bessel layer
#'
#' Simulates Brownian motion
#'
#'
#' @param bm a Brownian motion object from which simulation should continue.
#'   Note the object is updated in place.
#'
#' @return the Brownian motion object which was passed in argument \code{bm} is
#'   updated in place and returned, enabling chaining of commands with
#'   dplyr (and other) style pipes.
#'
#' @export
sim.condbbbessel <- function(bm, s, t, q = NULL, q.grid = NULL, label = names(q)) {
  # Arg types & combos
  if(!("BrownianMotion" %in% class(bm))) {
    stop("bm argument must be a BrownianMotion object.")
  }
  if(as.integer(is.null(q)) + as.integer(is.null(q.grid)) != 1) {
    stop("Exactly one of q or q.grid must be specified")
  }
  if(!is.null(q.grid) && !is.intscalar(q.grid)) {
    stop("q.grid must be either NULL or an integer number of grid points.")
  }
  if(!is.null(q) && !is.realvector(q)) {
    stop("q must be either NULL or a vector of times.")
  }
  if(!is.realscalar(t)) {
    stop("t must be a scalar.")
  }
  if(!is.realscalar(s)) {
    stop("s must be a scalar.")
  }

  # Checks
  if(!(s %in% bm$t)) {
    stop("s not found in bm path.")
  }
  if(!(t %in% bm$t)) {
    stop("t not found in bm path.")
  }
  bessel.bb <- bm$layers[bm$layers$type == "bessel-bb",]
  if(!(t %in% bessel.bb$t.u)) {
    stop("t is not the time of a realised maximum/minimum in a bessel-bb layer.") # "Dr Bessel is turning is his grave as we speak!  Don't associate my name with any of this!" - M. Pollock (2021)
  }
  # Do we then need to check that s is the start of this same bessel layer?
  # Actually this is taken care of by the following which allows no intermediate obs
  s_idx <- match(s, bm$t)
  t_idx <- match(t, bm$t)
  if(t_idx != s_idx+1) {
    stop("Cannot have other brownian motion observations between times s and t")
  }

  # Make grid if it wasn't supplied
  if(is.null(q)) {
    q <- head(tail(seq(s, t, length.out = q.grid+2), -1), -1)
  }
  if(any(q < s) || any(q > t)) {
    stop("All q must lie between s and t")
  }
  # Label checks now we have q
  assert.bmlabel(label, q)
  if(!is.null(label) && length(label) == 1) {
    label <- rep(label, length(q))
  }
  if(is.unsorted(q)) {
    o <- order(q)
    label <- label[o]
    q <- q[o]
  }
  # Eliminate times we know
  q2 <- setdiff(q, bm$t)
  label <- label[q2 %in% q]
  q <- q2

  for(qq in q) {
    bm.res <- sim.condbbbessel_(bm, s_idx, qq, t_idx, label[qq==q])
    s_idx <- s_idx+1
    t_idx <- t_idx+1
  }

  bm.res$layers <- bm.res$layers[order(bm$layers$t.l),]
  invisible(bm.res)
}

sim.condbbbessel_ <- function(bm, s_idx, q, t_idx, label) {
  bb.s_idx <- match(bm$t[s_idx], bm$bb.local$t)
  bb.t_idx <- match(bm$t[t_idx], bm$bb.local$t)

  Bs <- bm$W_t[s_idx]
  Bt <- bm$W_t[t_idx]
  Ws <- bm$bb.local$W_t[bb.s_idx]
  Wt <- bm$bb.local$W_t[bb.t_idx]
  s <- bm$t[s_idx]
  t <- bm$t[t_idx]

  sim.condbessel_(bm$bb.local, bb.s_idx, q, bb.t_idx, label)
  Wq <- bm$bb.local$W_t[bb.s_idx+1]
  Bq <- Bs + (Wq-Ws) + (q-s)/(t-s)*((Bt-Bs)-(Wt-Ws))

  bm$t <- c(bm$t[1:s_idx],
            q,
            bm$t[t_idx:length(bm$t)])
  bm$W_t <- c(bm$W_t[1:s_idx],
              Bq,
              bm$W_t[t_idx:length(bm$W_t)])

  cur.layer <- which(bm$layers$t.l == s & bm$layers$t.u == t)
  i <- nrow(bm$bb.local$layers)-1
  bm$layers <- add_row(bm$layers,
                       type = "bessel-bb",
                       t.l = bm$bb.local$layers$t.l[i],
                       t.u = bm$bb.local$layers$t.u[i],
                       # Ld = Bs + bm$bb.local$layers$Ld[i] - Ws,
                       # Uu = Bs + bm$bb.local$layers$Uu[i] - Ws + ((Bq-Bs)-(Wq-Ws)),
                       # Lu = Bs + bm$bb.local$layers$Ld[i] - Ws + ((Bq-Bs)-(Wq-Ws)),
                       # Ud = Bs + bm$bb.local$layers$Uu[i] - Ws,
                       Ld = Bs - Ws,
                       Uu = Bs - Ws + ((Bq-Bs)-(Wq-Ws)),
                       Lu = NA,
                       Ud = NA,
                       Lu.hard = bm$bb.local$layers$Lu.hard[i], #ifelse(head(c(0, x.new), -1) > x.new, TRUE, FALSE)
                       Ud.hard = bm$bb.local$layers$Ud.hard[i])
  bm$layers <- add_row(bm$layers,
                       type = "bessel-bb",
                       t.l = bm$bb.local$layers$t.l[i+1],
                       t.u = bm$bb.local$layers$t.u[i+1],
                       # Ld = Bq + bm$bb.local$layers$Ld[i+1] - Wq,
                       # Uu = Bq + bm$bb.local$layers$Uu[i+1] - Wq + ((Bt-Bq)-(Wt-Wq)),
                       # Lu = Bq + bm$bb.local$layers$Ld[i+1] - Wq + ((Bt-Bq)-(Wt-Wq)),
                       # Ud = Bq + bm$bb.local$layers$Uu[i+1] - Wq,
                       Ld = Bq - Wq,
                       Uu = Bq - Wq + ((Bt-Bq)-(Wt-Wq)),
                       Lu = NA,
                       Ud = NA,
                       Lu.hard = bm$bb.local$layers$Lu.hard[i+1], #ifelse(head(c(0, x.new), -1) > x.new, TRUE, FALSE)
                       Ud.hard = bm$bb.local$layers$Ud.hard[i+1])
  bm$layers <- bm$layers[-cur.layer,]

  # Sort the aux path since sorting usually happens in the wrapper which was
  # never called from the aux path's perspective ... our main path will be sorted
  # externally
  bm$bb.local$layers <- bm$bb.local$layers[order(bm$bb.local$layers$t.l),]

  add.labels_(bm, "user", q)
  add.labels_(bm, label, q)

  bm
}
