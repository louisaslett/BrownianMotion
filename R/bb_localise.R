bb.localise_ <- function(bm, s, t, refine, mult, prefer) {
  theta <- mult*sqrt(t-s)

  x <- bm$W_t[match(s, bm$t)]
  y <- bm$W_t[match(t, bm$t)]

  # Force points in auxiliary path
  new.bm <- create.bm(bm$t[1],
                      bm$W_t[1],
                      refine,
                      mult,
                      prefer)
  new.bm$t <- c(s)
  new.bm$W_t <- c(0)

  # Sim aux first passages until exceed t-s, then infill @ t-s
  first.passage(new.bm, delta = theta/2)
  while(max(new.bm$t) < t) {
    first.passage(new.bm, delta = theta/2)
  }
  if(!(t %in% new.bm$t)) {
    # t is not the end point, so conditionally simulate at t and then chop off the end
    sim(new.bm, t, refine, mult, prefer)
    delete.skeleton(new.bm, l = t)
  }

  # Copy the aux path into the master path bb localisation store
  # First check is there already a realisation at s?  If so, we need to shift
  # our new path up/down to match
  if(s %in% bm$bb.local$t) {
    new.bm$layers[,4:7] <- new.bm$layers[,4:7] + bm$bb.local$W_t[match(s, bm$bb.local$t)]
    new.bm$W_t <- new.bm$W_t + bm$bb.local$W_t[match(s, bm$bb.local$t)]
    # remove s from new bm so we don't now duplicate this obs
    new.bm$t <- new.bm$t[-1]
    new.bm$W_t <- new.bm$W_t[-1]
  }
  if(t %in% bm$bb.local$t) {
    bm$bb.local$layers[bm$bb.local$layers$t.l >= t,4:7] <- bm$bb.local$layers[bm$bb.local$layers$t.l >= t,4:7] + (new.bm$W_t[match(t, new.bm$t)] - bm$bb.local$W_t[match(t, bm$bb.local$t)])
    bm$bb.local$W_t[bm$bb.local$t >= t] <- bm$bb.local$W_t[bm$bb.local$t >= t] + (new.bm$W_t[match(t, new.bm$t)] - bm$bb.local$W_t[match(t, bm$bb.local$t)])
    # remove t from new bm so we don't now duplicate this obs
    new.bm$t <- head(new.bm$t, -1)
    new.bm$W_t <- head(new.bm$W_t, -1)
  }
  # Now do the copy, maintaining temporal ordering of all objects
  bm$bb.local$W_t <- c(bm$bb.local$W_t[bm$bb.local$t<=s],
                       new.bm$W_t,
                       bm$bb.local$W_t[bm$bb.local$t>=t])
  bm$bb.local$t <- c(bm$bb.local$t[bm$bb.local$t<=s],
                     new.bm$t,
                     bm$bb.local$t[bm$bb.local$t>=t])
  bm$bb.local$layers <- rbind(bm$bb.local$layers, new.bm$layers)
  bm$bb.local$layers <- bm$bb.local$layers[order(bm$bb.local$layers$t.l),]

  # The end points are removed before proceeding to transform if they were not removed for the above
  if(s %in% new.bm$t) {
    new.bm$t <- new.bm$t[-1]
    new.bm$W_t <- new.bm$W_t[-1]
  }
  if(t %in% new.bm$t) {
    new.bm$t <- head(new.bm$t, -1)
    new.bm$W_t <- head(new.bm$W_t, -1)
  }

  # Transform path back to master path
  trans.B_t <- rep(0, length(new.bm$t))
  B.last <- x
  t.last <- s
  W.last <- bm$bb.local$W_t[match(s, bm$bb.local$t)]
  W.end <- bm$bb.local$W_t[match(t, bm$bb.local$t)]
  if(length(new.bm$t)>0) {
    for(i in 1:length(new.bm$t)) {
      trans.B_t[i] <- B.last +
        (new.bm$W_t[i] - W.last) +
        (new.bm$t[i]-t.last)/(t-t.last)*((y-B.last)-(W.end-W.last))

      bm$layers <- add_row(bm$layers,
                           type = "localised-bb",
                           t.l = new.bm$layers$t.l[i],
                           t.u = new.bm$layers$t.u[i],
                           # Ld = B.last + new.bm$layers$Ld[i] - W.last,
                           # Uu = B.last + new.bm$layers$Uu[i] - W.last + ((trans.B_t[i]-B.last)-(new.bm$W_t[i]-W.last)),
                           # Lu = B.last + new.bm$layers$Ld[i] - W.last + ((trans.B_t[i]-B.last)-(new.bm$W_t[i]-W.last)),
                           # Ud = B.last + new.bm$layers$Uu[i] - W.last,
                           Ld = B.last - W.last,
                           Uu = B.last - W.last + ((trans.B_t[i]-B.last)-(new.bm$W_t[i]-W.last)),
                           Lu = NA,
                           Ud = NA,
                           Lu.hard = new.bm$layers$Lu.hard[i], #ifelse(head(c(0, x.new), -1) > x.new, TRUE, FALSE)
                           Ud.hard = new.bm$layers$Ud.hard[i])

      B.last <- trans.B_t[i]
      t.last <- new.bm$t[i]
      W.last <- new.bm$W_t[i]
    }
  } else {
    i <- 0
  }
  # Final layer segment will be missing so add
  bm$layers <- add_row(bm$layers,
                       type = paste0(new.bm$layers$type[i+1], "-bb"),
                       t.l = new.bm$layers$t.l[i+1],
                       t.u = new.bm$layers$t.u[i+1],
                       # Ld = B.last + new.bm$layers$Ld[i+1] - W.last,
                       # Uu = B.last + new.bm$layers$Uu[i+1] - W.last + ((y-B.last)-(W.end-W.last)),
                       # Lu = B.last + new.bm$layers$Ld[i+1] - W.last + ((y-B.last)-(W.end-W.last)),
                       # Ud = B.last + new.bm$layers$Uu[i+1] - W.last,
                       Ld = B.last - W.last,
                       Uu = B.last - W.last + ((y-B.last)-(W.end-W.last)),
                       Lu = NA,
                       Ud = NA,
                       Lu.hard = new.bm$layers$Lu.hard[i+1], #ifelse(head(c(0, x.new), -1) > x.new, TRUE, FALSE)
                       Ud.hard = new.bm$layers$Ud.hard[i+1])
  bm$W_t <- c(bm$W_t[bm$t<=s],
              trans.B_t,
              bm$W_t[bm$t>=t])
  bm$t <- c(bm$t[bm$t<=s],
            new.bm$t,
            bm$t[bm$t>=t])
  bm$layers <- bm$layers[order(bm$layers$t.l),]
  bm$labels[['fpt']] <- c(bm$labels[['fpt']], new.bm$labels[['fpt']][new.bm$labels[['fpt']]>=s & new.bm$labels[['fpt']]<=t])
  bm$labels[['internal']] <- c(bm$labels[['internal']], new.bm$labels[['internal']][new.bm$labels[['internal']]>=s & new.bm$labels[['internal']]<=t])

  rm(new.bm)
  invisible(bm)
}
