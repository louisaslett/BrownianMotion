#' BrownianMotion: exact simulation of (conditioned) Brownian motion
#'
#' Simulation of one- and multi-dimensional Brownian motion, including
#' Brownian bridges, first passage times, and layers (spatial bounds which
#' contain the path over an interval of time). Simulation is exact: each new
#' point is drawn conditionally on everything already known about the path,
#' including any layers.
#'
#' A new path can be created with \code{\link{create_bm}} and is then extended
#' or filled in with \code{\link{sim}}. The object is updated in place, so calls
#' can be chained with pipes. Its contents can be extracted by indexing (see
#' \code{\link[=BrownianMotion-indexing]{[.BrownianMotion}}) and plotted
#' directly using \code{\link[=plot.BrownianMotion]{plot}}.
#'
#' @section Main functions:
#' \describe{
#'   \item{Creating and simulating}{\code{\link{create_bm}}, \code{\link{sim}},
#'     \code{\link{add_points}}}
#'   \item{Layers}{\code{\link{layers}}, \code{\link{refine}},
#'     \code{\link{first_passage}}}
#'   \item{Editing a path}{\code{\link{add_segment}},
#'     \code{\link{delete_skeleton}}, \code{\link{concat_bm}}}
#'   \item{Labels}{\code{\link{add_labels}}, \code{\link{delete_labels}}}
#' }
#'
#' @section Layer types:
#' Three types of layer are available (see \code{\link{layers}}), which follow
#' (broadly) the approaches in the references below:
#' \describe{
#'   \item{Bessel}{Beskos, Papaspiliopoulos and Roberts (2008).}
#'   \item{Localised}{Devroye (2009), following the implementation in Pollock,
#'     Fearnhead, Johansen and Roberts (2020).}
#'   \item{Intersection}{Pollock, Johansen and Roberts (2016).}
#' }
#'
#' A full introduction is given in the tutorial vignette:
#' \code{vignette("tutorial", package = "BrownianMotion")}.
#'
#' @references
#' Beskos, A., Papaspiliopoulos, O. and Roberts, G.O. (2008). A factorisation
#' of diffusion measure and finite sample path constructions.
#' \emph{Methodology and Computing in Applied Probability}, 10, 85--104.
#'
#' Devroye, L. (2009). On exact simulation algorithms for some distributions
#' related to Jacobi theta functions. \emph{Statistics and Probability
#' Letters}, 79, 2251--2259.
#'
#' Pollock, M., Fearnhead, P., Johansen, A.M. and Roberts, G.O. (2020).
#' Quasi-stationary Monte Carlo methods and the ScaLE algorithm (with
#' discussion). \emph{Journal of the Royal Statistical Society, Series B
#' (Statistical Methodology)}, 82, 1--59.
#'
#' Pollock, M., Johansen, A.M. and Roberts, G.O. (2016). On the exact and
#' \eqn{\varepsilon}-strong simulation of (jump) diffusions. \emph{Bernoulli},
#' 22(2), 794--856.
#'
#' @examples
#' set.seed(1)
#' bm <- create_bm()
#' bm |> sim(1:4) |> layers(0, 2, type = "bessel") |> first_passage(delta = 1)
#' bm[, "W_t"]
#' plot(bm)
"_PACKAGE"
