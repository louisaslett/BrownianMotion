#' Create localised layer information around a Brownian bridge
#'
#' @param s start time
#' @param t end time
#'
#' The function will create localised layers across all currently layer-free Brownian bridges between times s and t
#'
#' @export
bb.localise <- function(bm, s, t, mult = 1) {
  # Arg types & combos
  if(!("BrownianMotion" %in% class(bm))) {
    stop("bm argument must be a BrownianMotion object.")
  }
  if(!is.realscalar(t)) {
    stop("t must be a scalar.")
  }
  if(!is.realscalar(s)) {
    stop("s must be a scalar.")
  }
  if(s >= t) {
    stop("s must be less than t.")
  }
  if(mult <= 0) {
    stop("mult must be strictly positive.")
  }


  # It t is currently beyond the end of the simulated path then first extend the path to create a bridgeable endpoint
  if(t > max(bm$t)) {
    # Maybe do first passage instead until past t and then simulate at t?
    sim(bm, t)
  }
  # Get all observed times between limits
  tt <- bm$t[bm$t>=s & bm$t<=t]
  # Determine all pairs not currently in a layer
  if(length(tt) < 2) {
    stop("there are no pairs of observed times between the specified s and t.")
  }
  all.s <- c()
  all.t <- c()
  for(i in 1:(length(tt)-1)) {
    # Only add if the interval is not in a layer
    mid <- mean(tt[i:(i+1)])
    if(!any(mid >= bm$layers$t.l & mid <= bm$layers$t.u)) {
      all.s <- c(all.s, tt[i])
      all.t <- c(all.t, tt[i+1])
    }
  }

  # Now perform localisation to all un-layered Brownian bridges we've discovered
  for(i in 1:length(all.s)) {
    bb.localise_(bm, all.s[i], all.t[i], mult)
  }

  invisible(bm)
}

bb.localise_ <- function(bm, s, t, mult) {
  theta <- mult*sqrt(t-s)

  x <- bm$W_t[match(s, bm$t)]
  y <- bm$W_t[match(t, bm$t)]

  # Force points in auxiliary path
  new.bm <- create.bm()
  new.bm$t <- c(s)
  new.bm$W_t <- c(0)

  # Sim aux first passages until exceed t-s, then infill @ t-s
  first.passage(new.bm, delta = theta/2)
  while(max(new.bm$t) < t) {
    first.passage(new.bm, delta = theta/2)
  }
  if(!(t %in% new.bm$t)) {
    # t is not the end point, so conditionally simulate at t and then chop off the end
    sim(new.bm, t)
    new.bm$t <- head(new.bm$t, -1)
    new.bm$W_t <- head(new.bm$W_t, -1)
    new.bm$layers <- new.bm$layers[-nrow(new.bm$layers),]
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

  # Transform path back to master path
  # for(i in 2:length(new.bm$t[new.bm$t<(t-s)])) {
  #   bm$t <- c(bm$t,
  #             s+new.bm$t[i])
  #   bm$W_t <- c(bm$W_t,
  #               x +
  #                 new.bm$W_t[i] - new.bm$W_t[i-1] +
  #                 (new.bm$t[i]-new.bm$t[i-1])/(t-s-new.bm$t[i-1])*((y-x)-W_T) + x)
  # }
  #
  #
  # new.bb.local <- cbind(new.bm$layers,
  #                       W_t.l = head(new.bm$W_t,-1),
  #                       W_t.u = tail(new.bm$W_t,-1))
  #
  #
  # bm$layers <- add_row(bm$layers,
  #                      type = "localised-bb",
  #                      t.l = s,
  #                      t.u = t,
  #                      Ld = NA,
  #                      Uu = NA,
  #                      Lu = NA,
  #                      Ud = NA,
  #                      Lu.hard = TRUE,
  #                      Ud.hard = TRUE)
  # bm$bb.local <- rbind(bm$bb.local, new.bb.local[new.bb.local$t.u <= t-s,])

  rm(new.bm)
  invisible(bm)
}
