#' Simulate an intersection layer
#'
#'
#'
#' @export
intersection.layers <- function(bm, s, t, refine = bm$refine, mult = bm$mult, prefer = bm$prefer) {
  ## NOTE TO LOUIS: This shares a lot of setup with Bessel Layers -- pull out into utility code (except the find Bessel layers part)
  if(!("BrownianMotion" %in% class(bm))) {
    stop("bm argument must be a BrownianMotion object.")
  }
  if(!is.realscalar(s) || !is.realscalar(t) || !is.realscalar(mult)) {
    stop("s, t and mult must be real scalars.")
  }
  if(s >= t) {
    stop("s must be less than t.")
  }
  if(mult <= 0) {
    stop("mult must be strictly positive.")
  }

  # Are these endpoints part of the skeleton?
  # If not, and the point is beyond the end, simulate it.  Otherwise, error
  if(!(s %in% bm$t)) {
    sim(bm, s, refine, mult, prefer)
  }
  if(!(t %in% bm$t)) {
    sim(bm, t, refine, mult, prefer)
  }

  # Find all intervals between s and t not currently inside any other layer or inside a Bessel layer
  s.i <- which(s == bm$t)
  t.i <- which(t == bm$t)
  all.pairs.noL <- matrix(c(bm$t[s.i:(t.i-1)], bm$t[(s.i+1):t.i]), byrow = FALSE, ncol = 2)
  all.pairs.BL  <- matrix(c(bm$t[s.i:(t.i-1)], bm$t[(s.i+1):t.i]), byrow = FALSE, ncol = 2)
  # No layer
  mid <- (all.pairs.noL[,2]+all.pairs.noL[,1])/2
  incl <- rep(TRUE, nrow(all.pairs.noL))
  for(i in 1:length(mid)) {
    if(any(mid[i] >= bm$layers$t.l & mid[i] <= bm$layers$t.u) ||
       any(mid[i] >= bm$user.layers$t.l & mid[i] <= bm$user.layers$t.u)) {
      incl[i] <- FALSE
    }
  }
  all.pairs.noL <- all.pairs.noL[incl,,drop = FALSE]
  # cat("No layer:"); print(all.pairs.noL)
  # Bessel layer
  mid <- (all.pairs.BL[,2]+all.pairs.BL[,1])/2
  incl <- rep(FALSE, nrow(all.pairs.BL))
  for(i in 1:length(mid)) {
    if(any(mid[i] >= bm$layers$t.l & mid[i] <= bm$layers$t.u & bm$layers$type == "bessel")) {
      incl[i] <- TRUE
    }
  }
  all.pairs.BL <- all.pairs.BL[incl,,drop = FALSE]
  # cat("Bessel layer:"); print(all.pairs.BL)

  if(nrow(all.pairs.noL) > 0) {
    for(i in 1:nrow(all.pairs.noL)) {
      intersection.layers.BLIL_(bm, all.pairs.noL[i,1], all.pairs.noL[i,2], mult)
    }
  }
  if(nrow(all.pairs.BL) > 0) {
    for(i in 1:nrow(all.pairs.BL)) {
      intersection.layers.IL_(bm, all.pairs.BL[i,1], all.pairs.BL[i,2], mult)
    }
  }

  bm$layers <- bm$layers[order(bm$layers$t.l),]
  invisible(bm)
}



intersection.layers.IL_ <- function(bm, s, t, mult) {
  BL <- which(s == bm$layers$t.l & t == bm$layers$t.u & bm$layers$type == "bessel")
  if(length(BL) != 1) {
    stop("Unable to identify bessel layer.")
  }

  res <- intersection.layers.sim_(bm, s, t, bm$W_t[match(s, bm$t)], bm$W_t[match(t, bm$t)],
                                  Ll = bm$layers[BL,"Ld",drop=TRUE],
                                  Lu = bm$layers[BL,"Lu",drop=TRUE],
                                  Ul = bm$layers[BL,"Ud",drop=TRUE],
                                  Uu = bm$layers[BL,"Uu",drop=TRUE],
                                  act = 2) # Just need anything >1 to indicate it is not already intersection

  bm$layers <- bm$layers[-BL,]
  bm$layers <- add_row(bm$layers,
                       type = "intersection",
                       t.l = s,
                       t.u = t,
                       Ld = res$Ll,
                       Uu = res$Uu,
                       Lu = res$Lu,
                       Ud = res$Ul,
                       Lu.hard = TRUE,
                       Ud.hard = TRUE)

  bm
}



intersection.layers.BLIL_ <- function(bm, s, t, mult) {
  # Simulate bessel layer before additional layer simulation
  BL <- bessel.layers.sim_(bm, s, t, bm$W_t[match(s, bm$t)], bm$W_t[match(t, bm$t)], mult)

  res <- intersection.layers.sim_(bm, s, t, bm$W_t[match(s, bm$t)], bm$W_t[match(t, bm$t)],
                                  Ll = BL$xb - BL$au,
                                  Lu = BL$xb - BL$al,
                                  Ul = BL$yb + BL$al,
                                  Uu = BL$yb + BL$au,
                                  act = BL$act)

  bm$layers <- add_row(bm$layers,
                       type = "intersection",
                       t.l = s,
                       t.u = t,
                       Ld = res$Ll,
                       Uu = res$Uu,
                       Lu = res$Lu,
                       Ud = res$Ul,
                       Lu.hard = TRUE,
                       Ud.hard = TRUE)

  bm
}



intersection.layers.sim_ <- function(bm, s, t, x, y, Ll, Lu, Ul, Uu, act) {
  xb <- min(x,y)
  yb <- max(x,y)

  # Additional layer simulation
  Di <- 0
  if(act==1) {
    Di <- 1
  } else {
    bbind <- deind <- 0
    m2 <- m3 <- 5
    u2 <- runif(1, 0, 1)
    while(deind == 0) {
      debd <- eabetaC_(m2, s, t, x, y, Ll, xb, yb, Uu)[2:1]
      if(debd[2] <= 0) {
        m2 <- m2+2
      } else {
        deind <- 1
      }
    }
    while(bbind == 0) {
      nubd <- eabetaC_(m3, s, t, x, y, Ll, Lu, Ul, Uu)[1:2]
      if(u2 <= (nubd[1]/debd[1])) {
        bbind <- 1
        Di <- 1
      } else {
        if(u2 >= (nubd[2]/debd[2])) {
          bbind <- 1
          u3 <- runif(1, 0, 1)
          if(u3<1/2) {
            Uu <- Ul
            Ul <- yb
            Di <- 2
          } else {
            Ll <- Lu
            Lu <- xb
            Di <- 3
          }
        } else {
          m2 <- m2+2
          m3 <- m3+2
          debd <- eabetaC_(m2, s, t, x, y, Ll, xb, yb, Uu)[2:1]
        }
      }
    }
  }

  list(s = s, t = t,x = x, y = y, Ll = Ll, Lu = Lu, Ul = Ul, Uu = Uu)
}

