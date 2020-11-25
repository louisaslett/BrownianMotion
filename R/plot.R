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
  p <- ggplot() +
    geom_line(aes(x = t, y = x), tibble(t = bm$t, x = bm$W_t), colour = "grey") +
    geom_point(aes(x = t, y = x), tibble(t = bm$t, x = bm$W_t), colour = "black", size = 0.1)

  if(nrow(bm$bounds) > 0) {
    p <- p +
      geom_segment(aes(x = t.l, xend = t.u, y = L, yend = L), bm$bounds, colour = "red") + #, linetype = "dashed") +
      geom_segment(aes(x = t.l, xend = t.u, y = U, yend = U), bm$bounds, colour = "red") + #, linetype = "dashed") +
      geom_segment(aes(x = t.l, xend = t.u, y = u, yend = u), bm$bounds, colour = "red") +
      geom_segment(aes(x = t.l, xend = t.u, y = l, yend = l), bm$bounds, colour = "red")
  }

  if(nrow(bm$bessel.layers) > 0) {
    p <- p +
      geom_segment(aes(x = t.l, xend = t.u, y = L, yend = L), bm$bessel.layers, colour = "blue", linetype = "dashed") +
      geom_segment(aes(x = t.l, xend = t.u, y = U, yend = U), bm$bessel.layers, colour = "blue", linetype = "dashed") +
      geom_segment(aes(x = t.l, xend = t.u, y = u, yend = u), bm$bessel.layers, colour = "blue") +
      geom_segment(aes(x = t.l, xend = t.u, y = l, yend = l), bm$bessel.layers, colour = "blue")
  }

  print(p)
}

#' @export
points.BrownianMotion <- function(x, y, ...) {
  points(x$t, x$W_t, type = "l", ...)
}
