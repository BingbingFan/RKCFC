#' @title Reassign one subject to its best-fitting class
#'
#' @description
#' Reconstructs subject \code{m}'s observations under each class-specific fit and
#' returns the class giving the smallest sum of squared residuals. This is
#' Equation (12) of the accompanying manuscript, and the reassignment half of the
#' K-centers iteration.
#'
#' @details
#' Because the class-specific fits are held fixed while every subject is
#' reassigned, this step minimises the objective
#' \eqn{J(c) = \sum_i \sum_j \{Y_{ij} - \hat{Y}_{ij,(c_i)}\}^2} exactly, and is
#' therefore monotone by construction. The re-estimation step
#' (\code{\link{univar_FPCA_k}}) is not; see \code{\link{RKCFC}}.
#'
#' With \code{LOO = TRUE} the fit for the class subject \code{m} currently
#' belongs to is re-estimated with \code{m} removed, so that the subject's own
#' data does not bias the decision to keep it where it already is. This removes a
#' real optimistic bias and costs one FPCA refit per subject per iteration, which
#' dominates the running time. Fits for the other classes are reused unchanged,
#' since \code{m} did not contribute to them.
#'
#' @param m Index of the subject to reassign.
#' @param ClustIds A list of length \code{K} giving the subject indices in each
#'   class, as returned internally from the current partition.
#' @param Ly A list of \emph{n} numeric vectors of measurements.
#' @param Lt A list of \emph{n} numeric vectors of observation times.
#' @param K Number of classes.
#' @param ref_hatyc The class-specific fits from \code{\link{univar_FPCA_k}}.
#' @param optns Estimation options, as for \code{\link{RKCFC}}.
#' @param robust Logical; the eigenspace estimator, as for \code{\link{RKCFC}}.
#' @param LOO Logical; leave subject \code{m} out when refitting its own class.
#'
#' @return A list with components
#'   \item{clust_i}{the index of the best-fitting class}
#'   \item{err}{the residual vector under the subject's current class}
#'   \item{hatyc}{the reconstruction under each class}
#'   \item{ss}{the residual sum of squares under each class}
#'
#' @seealso \code{\link{RKCFC}}, \code{\link{univar_predicted_yhat}}
#'
#' @export
univar_clust_update <- function(m, ClustIds, Ly, Lt, K, ref_hatyc, optns, robust, LOO) {
  ith_yy <- Ly[[m]]
  ith_tt <- Lt[[m]]
  hatyc <- list()
  err <- NULL
  for (u in 1:K) {                 ## index for cluster group
    uu <- ClustIds[[u]]            ## index of curves in cluster u
    if (LOO) {
      if (sum(m == uu) == 1) {     ## If m-th curve is in the group u
        yy <- lapply(uu, function(x) Ly[[x]])
        yy[[which(m == uu)]] <- NULL
        tt <- lapply(uu, function(x) Lt[[x]])
        tt[[which(m == uu)]] <- NULL
        # FPCA without subject m
        if (robust) {
          temp_fpca <- KFPCA(Lt = tt, Ly = yy, interval = optns$interval,
                             nK = optns$nK, bw = optns$bw, bwK = optns$bwK,
                             nRegGrid = optns$nRegGrid, fdParobj = optns$basis)
        } else {
          temp_fpca <- fdapace::FPCA(yy, tt, optns = optns)
        }
        yhatc_pred <- univar_predicted_yhat(temp_fpca, ith_yy, ith_tt, robust, optns)
      }

      if (sum(m == uu) == 0) {     ## If m-th curve is not in the group u
        yhatc_pred <- univar_predicted_yhat(ref_hatyc[[u]], ith_yy, ith_tt, robust, optns)
      }
    } else {
      yhatc_pred <- univar_predicted_yhat(ref_hatyc[[u]], ith_yy, ith_tt, robust, optns)
    }
    hatyc[[u]] <- yhatc_pred$fittedX
  }
  # calculate reconstruction error
  ss <- list()
  for (u in 1:K) {
    uu <- ClustIds[[u]]
    if (sum(m == uu) == 1) {
      err <- hatyc[[u]] - (ith_yy)
    }
    ss[[u]] <- sum((hatyc[[u]] - (ith_yy))^2)
  }
  clust_i <- which.min(ss)
  list(clust_i = clust_i, err = err, hatyc = hatyc, ss = ss)
}
