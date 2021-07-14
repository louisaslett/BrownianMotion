#' @export
plot.BrownianMotion <- function(x, y, ...) {
  # # TODO: update computing y.l and y.u to be over a restricted range is the ...
  # #       variadic argument contains xlim for a restricted time range plot
  # y.l <- min(x$W_t, x$bounds$l)
  # y.u <- max(x$W_t, x$bounds$u)
  # plot(x$t, x$W_t, type = "l", ylim = c(y.l, y.u), xlab = "Time", ylab = "State", ...)
  # if(nrow(x$bounds) > 0) {
  #   for(i in 1:nrow(x$bounds)) {
  #     points(unlist(x$bounds[i,1:2]), rep(unlist(x$bounds[i,3]), 2), type = "l", col = "red")
  #     points(unlist(x$bounds[i,1:2]), rep(unlist(x$bounds[i,4]), 2), type = "l", col = "red")
  #   }
  # }

  opts <- list(...)

  bm <- x

  p <- ggplot()

  jd <- c(1, which(bm$W_t != bm$W_tm), length(bm$t))
  for(j in 2:length(jd)) {
    if(j != length(jd))
      p <- p + geom_vline(xintercept = bm$t[jd[j]], linetype = "dotted")
    if(jd[j]-jd[j-1] > 0)
      p <- p + geom_line(aes(x = t, y = x), tibble(t = bm$t[jd[j-1]:jd[j]], x = c(bm$W_t[jd[j-1]:(jd[j]-1)], bm$W_tm[jd[j]])), colour = "grey")
  }

  p <- p + geom_point(aes(x = t, y = x), tibble(t = bm$t, x = bm$W_t), colour = "black", size = 0.1)
  jd <- jd[c(-1, -length(jd))]
  if(length(jd) > 0) {
    p <- p + geom_point(aes(x = t, y = x), tibble(t = bm$t[jd], x = bm$W_tm[jd]), colour = "black", size = 1, shape = 1)
  }

  localised <- bm$layers[bm$layers$type == "localised",]
  intersection <- bm$layers[bm$layers$type == "intersection",]
  bessel <- bm$layers[bm$layers$type == "bessel",]
  localised.bb <- bm$layers[bm$layers$type == "localised-bb",]
  intersection.bb <- bm$layers[bm$layers$type == "intersection-bb",]
  bessel.bb <- bm$layers[bm$layers$type == "bessel-bb",]

  if(is.null(opts[["hide.user"]]) || !opts[["hide.user"]]) {
    if(nrow(bm$user.layers) > 0) {
      p <- p +
        geom_segment(aes(x = t.l, xend = t.u, y = L, yend = L), bm$user.layers, colour = "green", size = 1.3) +
        geom_segment(aes(x = t.l, xend = t.u, y = U, yend = U), bm$user.layers, colour = "green", size = 1.3)
    }
  }

  if(nrow(localised) > 0) {
    p <- p +
      geom_segment(aes(x = t.l, xend = t.u, y = Ld, yend = Ld), localised, colour = "red") + #, linetype = "dashed") +
      geom_segment(aes(x = t.l, xend = t.u, y = Uu, yend = Uu), localised, colour = "red") + #, linetype = "dashed") +
      geom_segment(aes(x = t.l, xend = t.u, y = Ud, yend = Ud), localised, colour = "red", linetype = ifelse(localised$Ud.hard, "longdash", "dotted")) +
      geom_segment(aes(x = t.l, xend = t.u, y = Lu, yend = Lu), localised, colour = "red", linetype = ifelse(localised$Lu.hard, "longdash", "dotted"))
  }

  if(nrow(intersection) > 0) {
    p <- p +
      geom_segment(aes(x = t.l, xend = t.u, y = Ld, yend = Ld), intersection, colour = "blue") +
      geom_segment(aes(x = t.l, xend = t.u, y = Uu, yend = Uu), intersection, colour = "blue") +
      geom_segment(aes(x = t.l, xend = t.u, y = Ud, yend = Ud), intersection, colour = "blue", linetype = ifelse(intersection$Ud.hard, "longdash", "dotted")) +
      geom_segment(aes(x = t.l, xend = t.u, y = Lu, yend = Lu), intersection, colour = "blue", linetype = ifelse(intersection$Lu.hard, "longdash", "dotted"))
  }

  if(nrow(bessel) > 0) {
    p <- p +
      geom_segment(aes(x = t.l, xend = t.u, y = Ld, yend = Ld), bessel, colour = "purple") +
      geom_segment(aes(x = t.l, xend = t.u, y = Uu, yend = Uu), bessel, colour = "purple") +
      geom_segment(aes(x = t.l, xend = t.u, y = Ud, yend = Ud), bessel, colour = "purple", linetype = ifelse(bessel$Ud.hard, "longdash", "dotted")) +
      geom_segment(aes(x = t.l, xend = t.u, y = Lu, yend = Lu), bessel, colour = "purple", linetype = ifelse(bessel$Lu.hard, "longdash", "dotted"))
  }

  if(nrow(localised.bb) > 0) {
    localised.bb <- dplyr::left_join(localised.bb,
                                     bm$bb.local$layers,
                                     by = c("t.l", "t.u"),
                                     suffix = c("",".bb"))

    p <- p +
      geom_segment(aes(x = t.l, xend = t.u, y = Ld.bb+Ld, yend = Ld.bb+Uu), localised.bb, colour = "red") +
      geom_segment(aes(x = t.l, xend = t.u, y = Uu.bb+Ld, yend = Uu.bb+Uu), localised.bb, colour = "red") +
      geom_segment(aes(x = t.l, xend = t.u, y = Ud.bb+Ld, yend = Ud.bb+Uu), localised.bb, colour = "red", linetype = ifelse(localised.bb$Ud.hard, "longdash", "dotted")) +
      geom_segment(aes(x = t.l, xend = t.u, y = Lu.bb+Ld, yend = Lu.bb+Uu), localised.bb, colour = "red", linetype = ifelse(localised.bb$Lu.hard, "longdash", "dotted"))
    if(is.null(opts[["hide.user"]]) || !opts[["hide.user"]]) {
      p <- p +
        geom_segment(aes(x = t.l, xend = t.u, y = pmin(Ld.bb+Ld,Ld.bb+Uu), yend = pmin(Ld.bb+Ld,Ld.bb+Uu)), localised.bb, colour = "red", alpha = 0.75) +
        geom_segment(aes(x = t.l, xend = t.u, y = pmax(Uu.bb+Ld,Uu.bb+Uu), yend = pmax(Uu.bb+Ld,Uu.bb+Uu)), localised.bb, colour = "red", alpha = 0.75)
    }
  }

  if(nrow(intersection.bb) > 0) {
    intersection.bb <- dplyr::left_join(intersection.bb,
                                        bm$bb.local$layers,
                                        by = c("t.l", "t.u"),
                                        suffix = c("",".bb"))

    p <- p +
      geom_segment(aes(x = t.l, xend = t.u, y = Ld.bb+Ld, yend = Ld.bb+Uu), intersection.bb, colour = "blue") +
      geom_segment(aes(x = t.l, xend = t.u, y = Uu.bb+Ld, yend = Uu.bb+Uu), intersection.bb, colour = "blue") +
      geom_segment(aes(x = t.l, xend = t.u, y = Ud.bb+Ld, yend = Ud.bb+Uu), intersection.bb, colour = "blue", linetype = ifelse(intersection.bb$Ud.hard, "longdash", "dotted")) +
      geom_segment(aes(x = t.l, xend = t.u, y = Lu.bb+Ld, yend = Lu.bb+Uu), intersection.bb, colour = "blue", linetype = ifelse(intersection.bb$Lu.hard, "longdash", "dotted"))
    if(is.null(opts[["hide.user"]]) || !opts[["hide.user"]]) {
      p <- p +
        geom_segment(aes(x = t.l, xend = t.u, y = pmin(Ld.bb+Ld,Ld.bb+Uu), yend = pmin(Ld.bb+Ld,Ld.bb+Uu)), intersection.bb, colour = "blue", alpha = 0.75) +
        geom_segment(aes(x = t.l, xend = t.u, y = pmax(Uu.bb+Ld,Uu.bb+Uu), yend = pmax(Uu.bb+Ld,Uu.bb+Uu)), intersection.bb, colour = "blue", alpha = 0.75)
    }
  }

  if(nrow(bessel.bb) > 0) {
    bessel.bb <- dplyr::left_join(bessel.bb,
                                  bm$bb.local$layers,
                                  by = c("t.l", "t.u"),
                                  suffix = c("",".bb"))

    p <- p +
      geom_segment(aes(x = t.l, xend = t.u, y = Ld.bb+Ld, yend = Ld.bb+Uu), bessel.bb, colour = "purple") +
      geom_segment(aes(x = t.l, xend = t.u, y = Uu.bb+Ld, yend = Uu.bb+Uu), bessel.bb, colour = "purple") +
      geom_segment(aes(x = t.l, xend = t.u, y = Ud.bb+Ld, yend = Ud.bb+Uu), bessel.bb, colour = "purple", linetype = ifelse(bessel.bb$Ud.hard, "longdash", "dotted")) +
      geom_segment(aes(x = t.l, xend = t.u, y = Lu.bb+Ld, yend = Lu.bb+Uu), bessel.bb, colour = "purple", linetype = ifelse(bessel.bb$Lu.hard, "longdash", "dotted"))
    if(is.null(opts[["hide.user"]]) || !opts[["hide.user"]]) {
      p <- p +
        geom_segment(aes(x = t.l, xend = t.u, y = pmin(Ld.bb+Ld,Ld.bb+Uu), yend = pmin(Ld.bb+Ld,Ld.bb+Uu)), bessel.bb, colour = "purple", alpha = 0.75) +
        geom_segment(aes(x = t.l, xend = t.u, y = pmax(Uu.bb+Ld,Uu.bb+Uu), yend = pmax(Uu.bb+Ld,Uu.bb+Uu)), bessel.bb, colour = "purple", alpha = 0.75)
    }
  }

  if(!is.null(opts[["t.lim"]])) {
    t.lim <- opts[["t.lim"]]

    lyrs <- bm$layers[!(bm$layers$type %in% c("localised-bb","intersection-bb","bessel-bb")),]
    ymin <- min(lyrs[lyrs$t.u>t.lim[1] & lyrs$t.l<t.lim[2],]$Ld, bm$W_t[bm$t>t.lim[1] & bm$t<t.lim[2]])
    ymax <- max(lyrs[lyrs$t.u>t.lim[1] & lyrs$t.l<t.lim[2],]$Uu, bm$W_t[bm$t>t.lim[1] & bm$t<t.lim[2]])

    bb.lyrs <- rbind(localised.bb, intersection.bb, bessel.bb)
    if(nrow(bb.lyrs) > 0) {
      bb.lyrs$min <- pmin(bb.lyrs$Ld.bb+bb.lyrs$Ld,bb.lyrs$Ld.bb+bb.lyrs$Uu)
      bb.lyrs$max <- pmax(bb.lyrs$Uu.bb+bb.lyrs$Ld,bb.lyrs$Uu.bb+bb.lyrs$Uu)
      ymin <- min(ymin,
                  bb.lyrs[bb.lyrs$t.u>t.lim[1] & bb.lyrs$t.l<t.lim[2],]$min)
      ymax <- max(ymax,
                  bb.lyrs[bb.lyrs$t.u>t.lim[1] & bb.lyrs$t.l<t.lim[2],]$max)
    }

    print(p + coord_cartesian(xlim = opts[["t.lim"]], ylim = c(ymin,ymax)))
  } else {
    print(p)
  }
}

#' @export
points.BrownianMotion <- function(x, y, ...) {
  points(x$t, x$W_t, type = "l", ...)
}
