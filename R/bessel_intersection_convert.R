# These functions should accept any layer indexes and convert only if possible
# otherwise leaving everything in tact

bessel.to.intersection_ <- function(bm, lyr_idxs, mult) {
  res <- c()

  bm$layers <- bm$layers[order(bm$layers$t.l),]

  lyr_idxs_bessel <- lyr_idxs[which(bm$layers$type[lyr_idxs] == "bessel")]
  for(lyr_idx in lyr_idxs_bessel) {
    intersection.layers.IL_(bm,
                            bm$layers$t.l[lyr_idx],
                            bm$layers$t.u[lyr_idx],
                            mult)
    bm$layers <- bm$layers[order(bm$layers$t.l),]
    res <- c(res, lyr_idx)
  }

  lyr_idxs_bessel_bb <- lyr_idxs[which(bm$layers$type[lyr_idxs] == "bessel-bb")]
  for(lyr_idx in lyr_idxs_bessel_bb) {
    lyr_idx_aux <- which(bm$bb.local$layers$t.l == bm$layers$t.l[lyr_idx])
    intersection.layers.IL_(bm$bb.local,
                            bm$bb.local$layers$t.l[lyr_idx_aux],
                            bm$bb.local$layers$t.u[lyr_idx_aux],
                            mult)
    bm$bb.local$layers <- bm$bb.local$layers[order(bm$bb.local$layers$t.l),]
    bm$layers$type[lyr_idx] <- "intersection-bb"
    bm$layers$Lu.hard[lyr_idx] <- bm$bb.local$layers$Lu.hard[lyr_idx_aux]
    bm$layers$Ud.hard[lyr_idx] <- bm$bb.local$layers$Ud.hard[lyr_idx_aux]

    res <- c(res, lyr_idx)
  }

  res
}

intersection.to.bessel_ <- function(bm, lyr_idxs) {
  res <- c()

  lyr_idxs_intersection <- lyr_idxs[which(bm$layers$type[lyr_idxs] == "intersection")]
  for(lyr_idx in lyr_idxs_intersection) {
    s <- bm$layers$t.l[lyr_idx]
    t <- bm$layers$t.u[lyr_idx]
    x <- bm$W_t[which(bm$t==s)]
    y <- bm$W_t[which(bm$t==t)]
    if(bm$layers$Ud[lyr_idx] == max(x,y) && bm$layers$Lu[lyr_idx] == min(x,y)) {
      bm$layers$type[lyr_idx] <- "bessel"
      res <- c(res, lyr_idx)
    }
  }

  lyr_idxs_intersection_bb <- lyr_idxs[which(bm$layers$type[lyr_idxs] == "intersection-bb")]
  for(lyr_idx in lyr_idxs_intersection_bb) {
    s <- bm$layers$t.l[lyr_idx]
    t <- bm$layers$t.u[lyr_idx]
    lyr_idx_aux <- which(bm$bb.local$layers$t.l == s)
    x <- bm$bb.local$W_t[which(bm$bb.local$t==s)]
    y <- bm$bb.local$W_t[which(bm$bb.local$t==t)]
    if(bm$bb.local$layers$Ud[lyr_idx_aux] == max(x,y) && bm$bb.local$layers$Lu[lyr_idx_aux] == min(x,y)) {
      bm$bb.local$layers$type[lyr_idx_aux] <- "bessel"
      bm$layers$type[lyr_idx] <- "bessel-bb"

      res <- c(res, lyr_idx)
    }
  }

  res
}
