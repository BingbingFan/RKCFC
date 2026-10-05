#' @title Reconstruct a subject's trajectory under a given class fit
#'
#' @description
#' Evaluates Equation (11) of the accompanying manuscript: the class mean
#' function plus the subject's estimated principal component scores against the
#' class eigenfunctions, at the subject's own observation times.
#'
#' @details
#' This is reconstruction, not forecasting. The subject's scores are estimated
#' from the very observations being reconstructed, so the resulting error is an
#' in-sample quantity and will flatter a more flexible fit. For predictive
#' performance, withhold observations from the fit and score them here, or in
#' simulation score against the true latent trajectory.
#'
#' In the robust branch the scores come from \code{\link{predict.KFPCA}}. In the
#' non-robust branch they are the conditional expectation of Yao, Mueller and
#' Wang (2005), computed by \code{.rkcfc_nonrobust_scores}, which solves the
#' model covariance system directly and falls back to a ridge and then to a
#' generalised inverse rather than always forming a pseudo-inverse.
#'
#' @param FPCA_k A fitted class object: a \code{KFPCA} object when
#'   \code{robust = TRUE}, otherwise an \code{fdapace::FPCA} object.
#' @param ith_yy Numeric vector of the subject's measurements.
#' @param target_t Numeric vector of the subject's observation times.
#' @param robust Logical; the estimator that produced \code{FPCA_k}.
#' @param optns Estimation options, as for \code{\link{RKCFC}}.
#'
#' @return A list with components
#'   \item{fittedX}{the reconstructed values at \code{target_t}}
#'   \item{xi_Est_i}{the estimated principal component scores}
#'   \item{pi_i}{the eigenfunctions evaluated at \code{target_t}}
#'
#' @seealso \code{\link{univar_clust_update}}, \code{\link{rre}}
#'
#' @export
univar_predicted_yhat <- function(FPCA_k, ith_yy, target_t, robust, optns) {
  if (robust) {
    prd <- stats::predict(FPCA_k, newLt = list(target_t), newLy = list(ith_yy),
                          nK = optns$nK, more = TRUE)
    xi_Est_i <- prd$score_new
    pi_temp <- prd$FPC_dis_new
    fittedX <- as.numeric(prd$meanest_new + xi_Est_i %*% t(pi_temp))
  } else {
    n_M <- length(FPCA_k$lambda)
    lambda <- FPCA_k$lambda
    phi <- FPCA_k$phi
    mu <- FPCA_k$mu
    workGrid <- FPCA_k$workGrid
    sig2 <- FPCA_k$sigma2

    pi_temp <- vapply(seq_len(n_M), function(l) {
      stats::spline(workGrid, phi[, l], xout = target_t)$y
    }, numeric(length(ith_yy)))
    mu_ij <- stats::spline(workGrid, mu, xout = target_t)$y
    xi_Est_i <- .rkcfc_nonrobust_scores(
      pi_temp,
      lambda,
      sig2,
      ith_yy - mu_ij
    )
    fittedX <- as.numeric(mu_ij + pi_temp %*% xi_Est_i)
  }

  list(fittedX = fittedX, xi_Est_i = xi_Est_i, pi_i = pi_temp)
}
