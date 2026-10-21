# Build a plot fully (all patchwork panels) without drawing it to a device.
build_plot <- function(p) {
  # Use a null device so no Rplots.pdf is left behind
  withr::local_pdf(NULL)
  if (inherits(p, "patchwork")) {
    patchwork::patchworkGrob(p)
  } else {
    ggplot2::ggplot_build(p)
  }
}

# A plot call returns a ggplot that builds without errors or warnings.
expect_plot_ok <- function(expr) {
  expect_no_warning(p <- expr)
  expect_s3_class(p, "ggplot")
  expect_no_warning(build_plot(p))
  invisible(p)
}
