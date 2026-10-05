# Robust location and scatter for a score matrix.
#
# The multivariate extension, MRKCFC, uses the same estimators, with
# robust_scatter() named robust_stage2() there after the Happ and Greven second
# stage that combines per-variable scores. The numerical bodies are identical
# and must stay identical; only the name differs.
#
# Used by whiten_scores() for the score-space K-selection criteria, and by the
# robust Mahalanobis separation measure the simulation study reports alongside
# the raw mean shift.

#' @title Spatial median of a score matrix
#'
#' @description
#' The multivariate \eqn{L_1} median, computed by Weiszfeld's iteration from a
#' coordinatewise-median start. Unlike the mean it has a population target
#' without any moment assumption, which is what makes it usable on heavy-tailed
#' principal component scores.
#'
#' @param X A numeric matrix, one row per subject.
#' @param tol Convergence tolerance on the change in the estimate.
#' @param max_it Maximum number of iterations.
#'
#' @return A numeric vector of length \code{ncol(X)}.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(stats::rt(200, df = 2), ncol = 2)
#' spatial_median(X)
#'
#' @export
spatial_median <- function(X, tol = 1e-8, max_it = 200L) {
  m <- apply(X, 2, stats::median)
  for (i in seq_len(max_it)) {
    d <- sqrt(rowSums((X - rep(m, each = nrow(X)))^2))
    d[d < 1e-12] <- 1e-12
    m_new <- colSums(X / d) / sum(1 / d)
    if (sqrt(sum((m_new - m)^2)) < tol) return(m_new)
    m <- m_new
  }
  m
}

#' @title Spatial sign covariance matrix
#'
#' @description
#' The covariance of the unit vectors obtained by projecting deviations from the
#' spatial median onto the unit sphere. Under an elliptical law its eigenvectors
#' are consistent for those of the scatter matrix, which is the same argument
#' Kendall's tau uses at the operator level, so no moment assumption is required.
#'
#' @details
#' Eigen*vectors* are recovered; eigen*values* are not. The matrix has trace 1 by
#' construction, so its eigenvalues carry no information about the spread along
#' each direction. \code{\link{robust_scatter}} recovers that separately.
#'
#' @param X A numeric matrix, one row per subject.
#'
#' @return A \code{ncol(X)} by \code{ncol(X)} symmetric matrix.
#'
#' @references
#' Visuri, S., Koivunen, V. and Oja, H. (2000). Sign and rank covariance
#' matrices. \emph{Journal of Statistical Planning and Inference}, \bold{91}(2),
#' 557-575.
#'
#' @seealso \code{\link{robust_scatter}}
#'
#' @export
sscm <- function(X) {
  U <- X - rep(spatial_median(X), each = nrow(X))
  nrm <- sqrt(rowSums(U^2))
  keep <- nrm > 1e-12
  U <- U[keep, , drop = FALSE] / nrm[keep]
  crossprod(U) / nrow(U)
}

#' @title Robust eigenstructure of a score matrix
#'
#' @description
#' Eigenvectors from the spatial sign covariance matrix, with the spread along
#' each estimated direction recovered by a robust scale rather than by a second
#' moment. With \code{robust = FALSE} it reduces to the ordinary sample
#' second-moment eigendecomposition, which is the ablation arm.
#'
#' @details
#' The spatial sign covariance matrix has trace 1, so inverting or ranking by its
#' eigenvalues would put the result on an arbitrary scale. The spread along each
#' estimated direction is therefore recovered as \eqn{(\mathrm{IQR}/1.349)^2} of
#' the projected scores, which is Fisher-consistent for the variance at the
#' normal and exists without a moment assumption.
#'
#' \strong{The components are ordered by that recovered scale, not by the sign
#' covariance eigenvalue.} In the population the two orderings agree, since the
#' sign eigenvalues are monotone in the scatter eigenvalues, but in a finite
#' sample under heavy tails they do not: on a simulated arm with infinite
#' variance the sign ordering returned spreads of 8.54, 8.78 and 7.87, making
#' component 2 more important than component 1. Any use that truncates to the
#' leading components would then silently work in the wrong subspace.
#'
#' @param Xi A numeric \emph{n} by \emph{M} matrix of principal component scores.
#' @param M_keep Number of components to retain; defaults to all.
#' @param robust Logical; \code{TRUE} uses the spatial sign covariance matrix,
#'   \code{FALSE} the sample second-moment matrix.
#'
#' @return A list with components
#'   \item{vectors}{an \emph{M} by \code{M_keep} matrix of eigenvectors in columns}
#'   \item{values}{the spread along each retained direction, on the score scale}
#'
#' @seealso \code{\link{sscm}}, \code{\link{whiten_scores}}
#'
#' @examples
#' set.seed(1)
#' Xi <- matrix(stats::rt(300, df = 2), ncol = 3) %*% diag(c(3, 2, 1))
#' robust_scatter(Xi, M_keep = 2)$values
#'
#' @export
robust_scatter <- function(Xi, M_keep = ncol(Xi), robust = TRUE) {
  stopifnot(is.matrix(Xi), M_keep >= 1L, M_keep <= ncol(Xi))
  S <- if (robust) sscm(Xi) else crossprod(Xi) / nrow(Xi)
  e <- eigen(S, symmetric = TRUE)

  if (robust) {
    nu_all <- apply(Xi %*% e$vectors, 2, function(z) (stats::IQR(z) / 1.349)^2)
    keep <- order(nu_all, decreasing = TRUE)[seq_len(M_keep)]
    C <- e$vectors[, keep, drop = FALSE]
    nu <- nu_all[keep]
  } else {
    # eigen() already returns decreasing eigenvalues, so the moment branch is
    # ordered by construction.
    C <- e$vectors[, seq_len(M_keep), drop = FALSE]
    nu <- e$values[seq_len(M_keep)]
  }
  list(vectors = C, values = nu)
}
