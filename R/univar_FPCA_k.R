#' @title Class-specific functional principal component fits
#'
#' @description
#' Fits one functional principal component model per class, using the subjects
#' currently assigned to it. This is the estimation half of the K-centers
#' iteration; \code{\link{univar_clust_update}} is the reassignment half.
#'
#' @details
#' The \code{robust} switch is owned here rather than duplicated at each call
#' site, so that every part of the iteration estimates the eigenspace the same
#' way. Note that the two engines name their scores differently:
#' \code{\link{KFPCA}} returns \code{$score} and \code{fdapace::FPCA} returns
#' \code{$xiEst}.
#'
#' @param Ly A list of \emph{n} numeric vectors of measurements.
#' @param Lt A list of \emph{n} numeric vectors of observation times.
#' @param K Number of classes.
#' @param clust_f An integer vector of length \emph{n} holding the current class
#'   label of each subject.
#' @param robust Logical; \code{TRUE} uses \code{\link{KFPCA}}, \code{FALSE} uses
#'   \code{fdapace::FPCA}.
#' @param optns Estimation options, as for \code{\link{RKCFC}}.
#'
#' @return A list of length \code{K}, each element the fitted object for one
#'   class: a \code{KFPCA} object when \code{robust = TRUE}, otherwise an
#'   \code{fdapace::FPCA} object.
#'
#' @seealso \code{\link{RKCFC}}, \code{\link{univar_clust_update}}
#'
#' @export
univar_FPCA_k <- function(Ly, Lt, K, clust_f, robust, optns) {
  ClustIds <- .rkcfc_cluster_ids(clust_f, K)
  ref_hatyc <- list()
  for (u in 1:K) {
    uu <- ClustIds[[u]]
    yy <- lapply(uu, function(u) Ly[[u]])
    tt <- lapply(uu, function(u) Lt[[u]])

    if (robust) {
      # `meanrobust` defaults to FALSE via isTRUE(NULL), so an optns list that
      # predates the robust mean option -- which is every optns the published results
      # were produced with -- selects the least-squares mean exactly as before.
      ref_hatyc[[u]] <- KFPCA(Lt = tt, Ly = yy, interval = optns$interval,
                              nK = optns$nK, bw = optns$bw, bwK = optns$bwK,
                              meanrobust = isTRUE(optns$meanrobust),
                              scatter = if (is.null(optns$scatter)) "kendall"
                                        else optns$scatter,
                              nRegGrid = optns$nRegGrid, fdParobj = optns$basis)
    } else {
      ref_hatyc[[u]] <- fdapace::FPCA(yy, tt, optns = optns)
    }
  }

  ref_hatyc
}
