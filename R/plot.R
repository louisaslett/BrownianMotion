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

  p <- ggplot() +
    geom_line(aes(x = t, y = x), tibble(t = bm$t, x = bm$W_t), colour = "grey") +
    geom_point(aes(x = t, y = x), tibble(t = bm$t, x = bm$W_t), colour = "black", size = 0.1)

  localised <- bm$layers[bm$layers$type == "localised",]
  intersection <- bm$layers[bm$layers$type == "intersection",]
  bessel <- bm$layers[bm$layers$type == "bessel",]
  localised.bb <- bm$layers[bm$layers$type == "localised-bb",]
  intersection.bb <- bm$layers[bm$layers$type == "intersection-bb",]

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
    p <- p +
      geom_segment(aes(x = t.l, xend = t.u, y = Ld, yend = Lu), localised.bb, colour = "red") +
      geom_segment(aes(x = t.l, xend = t.u, y = Ud, yend = Uu), localised.bb, colour = "red") +
      geom_segment(aes(x = t.l, xend = t.u, y = pmin(Ld,Lu), yend = pmin(Ld,Lu)), localised.bb, colour = "red", alpha = 0.75) +
      geom_segment(aes(x = t.l, xend = t.u, y = pmax(Ud,Uu), yend = pmax(Ud,Uu)), localised.bb, colour = "red", alpha = 0.75)
  }

  if(nrow(intersection.bb) > 0) {
    p <- p +
      geom_segment(aes(x = t.l, xend = t.u, y = Ld, yend = Lu), intersection.bb, colour = "blue") +
      geom_segment(aes(x = t.l, xend = t.u, y = Ud, yend = Uu), intersection.bb, colour = "blue") +
      geom_segment(aes(x = t.l, xend = t.u, y = pmin(Ld,Lu), yend = pmin(Ld,Lu)), intersection.bb, colour = "blue", alpha = 0.75) +
      geom_segment(aes(x = t.l, xend = t.u, y = pmax(Ud,Uu), yend = pmax(Ud,Uu)), intersection.bb, colour = "blue", alpha = 0.75)
  }

  if(!is.null(opts[["t.lim"]])) {
    print(p + xlim(opts[["t.lim"]]))
  } else {
    print(p)
  }
}

#' @export
points.BrownianMotion <- function(x, y, ...) {
  points(x$t, x$W_t, type = "l", ...)
}
