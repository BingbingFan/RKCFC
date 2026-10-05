#' @title Initialise the K-centers iteration
#'
#' @description
#' Produces the starting partition for \code{\link{RKCFC}}, either by k-means on
#' the pooled functional principal component scores or by random assignment.
#'
#' @details
#' The K-centers iteration is a Lloyd-type scheme converging to a local optimum,
#' so the starting partition selects the basin and the basin largely selects the
#' answer. Two consequences follow.
#'
#' First, the estimator behind the start must match the estimator behind the
#' iteration. When \code{robust = TRUE} the scores come from
#' \code{\link{KFPCA}}, not from \code{fdapace::FPCA}. The breakdown point of a
#' composite estimator is bounded above by the minimum over its components, so a
#' non-robust start caps the robustness of the whole procedure however robust the
#' later stages are; and PACE requires finite fourth moments for its empirical
#' covariance surface to be consistent, whereas Kendall's tau does not, so on a
#' process without them a non-robust start has no reliable population target at
#' all. Versions of this package before 0.2.0 called \code{fdapace::FPCA}
#' unconditionally here and never received \code{robust}; results produced with
#' those versions are affected.
#'
#' Second, no single start should be trusted. Run several and keep the partition
#' attaining the smallest total squared reconstruction error.
#'
#' \code{nstart} is left at the \code{stats::kmeans} default of 1, faithfully to
#' the manuscript, which specifies only "using K-means".
#'
#' @param Ly A list of \emph{n} numeric vectors of measurements, one per subject.
#' @param Lt A list of \emph{n} numeric vectors of observation times.
#' @param K Number of classes.
#' @param init.method One of \code{"kmeans"}, \code{"random"} or
#'   \code{"two-staged"}. The last applies k-means and falls back to random
#'   assignment when k-means returns a class with fewer than ten curves or fewer
#'   than \code{K} non-empty classes, which is the condition that would otherwise
#'   make \code{\link{RKCFC}} refuse the fit.
#' @param robust Logical; \code{TRUE} scores by Kendall's tau, \code{FALSE} by
#'   the PACE covariance estimator.
#' @param optns Estimation options, as for \code{\link{RKCFC}}.
#'
#' @return An integer vector of length \emph{n} holding the starting class label
#'   of each subject.
#'
#' @seealso \code{\link{RKCFC}}
#'
#' @examples
#' \donttest{
#' data(pbc2)
#' pbc <- pbc2[!is.na(pbc2$serBilir), c("id", "year", "serBilir")]
#' dt <- fdapace::MakeFPCAInputs(IDs = pbc$id, tVec = pbc$year, yVec = pbc$serBilir)
#' start <- KCFC_init(dt$Ly, dt$Lt, K = 2, init.method = "random",
#'                    robust = FALSE, optns = list(dataType = "Sparse"))
#' table(start)
#' }
#'
#' @export
KCFC_init <- function(Ly, Lt, K, init.method, robust, optns) {
  init.method <- match.arg(init.method, c("kmeans", "random", "two-staged"))
  n <- length(Ly)

  random_labels <- function()
    sample.int(K, n, replace = TRUE, prob = rep(1 / K, K))

  if (init.method == "random") return(random_labels())

  ref <- univar_FPCA_k(Ly, Lt, 1L, rep(1L, n), robust, optns)[[1]]
  scores <- if (robust) ref$score else ref$xiEst
  clust0 <- stats::kmeans(scores, centers = K)$cluster

  # "or randomly assign the class memberships": a last resort, fired only when
  # k-means still returns a degenerate split, which is the condition that later
  # makes the fit be refused for having a group under ten curves.
  if (init.method == "two-staged" &&
      (min(table(clust0)) < 10L || length(unique(clust0)) < K)) {
    clust0 <- random_labels()
  }
  clust0
}
