#' Simulate an intersection layer
#'
#'
#'
#' @export
intersection.layers <- function(bm, s, t, refine = bm$refine, mult = bm$mult, prefer = bm$prefer, label = c(names(s), names(t))) {
  layers(bm, s, t, "intersection", refine, mult, prefer, label)
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

