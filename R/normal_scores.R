#' @title Normal-scores (rank) transform of sparse longitudinal measurements
#'
#' @description
#' Replaces the measurements by their normal scores: each observation is ranked
#' among all observations pooled across subjects and times, and the rank is
#' mapped through the standard normal quantile function. The resulting values are
#' marginally normal by construction, so a conventional PACE analysis applied to
#' them is not exposed to heavy marginal tails.
#'
#' @details
#' This is the cheapest of the robust options in this package and the one most
#' different in kind. It is a wrapper rather than an estimator: nothing inside
#' the fitting changes, and the transformed measurements are handed to whatever
#' method the caller chooses.
#'
#' \strong{It changes the estimand, not only the estimator, and that is not a
#' detail.} The eigenfunctions of the transformed process are the eigenfunctions
#' of a different process, so a mean function, a covariance surface and a
#' reconstruction error computed on normal scores are not on the scale of the
#' data and are not comparable with quantities computed on the raw scale. A
#' partition, by contrast, is invariant to any strictly monotone marginal
#' transform of this kind and \emph{is} comparable. Use this for clustering, and
#' report reconstruction on the original scale or not at all.
#'
#' The transform is applied to the POOLED observations rather than separately at
#' each time point. On a sparse irregular design each time point carries at most
#' a handful of measurements -- often exactly one -- so a per-time-point transform
#' would be estimated from almost no data and would mostly encode the sampling
#' pattern rather than the marginal law.
#'
#' @param Ly A \code{list} of \emph{n} numeric vectors, one per subject, holding
#'   the measurements. This is the same \code{Ly} the fitting functions take.
#' @param ties.method How ties are resolved, passed to \code{\link[base]{rank}}.
#'   The default \code{"average"} gives tied observations a common score.
#'
#' @return A \code{list} of the same length and shape as \code{Ly}, holding the
#'   normal scores. \code{NA} values are preserved in place and are excluded from
#'   the ranking rather than ranked.
#'
#' @examples
#' set.seed(1)
#' Ly <- lapply(1:20, function(i) stats::rt(sample(4:8, 1), df = 1.5))
#' Lz <- normal_scores(Ly)
#' # shape preserved, marginal tails removed
#' identical(lengths(Ly), lengths(Lz))
#' range(unlist(Lz))
#'
#' @seealso \code{\link{KFPCA}} for the estimator-side alternatives.
#' @export
normal_scores <- function(Ly, ties.method = "average") {
  if (!is.list(Ly)) stop("Ly must be a list of numeric vectors.", call. = FALSE)
  ni <- lengths(Ly)
  y  <- unlist(Ly, use.names = FALSE)
  if (!is.numeric(y)) stop("Ly must contain numeric vectors.", call. = FALSE)

  ok <- is.finite(y)
  N  <- sum(ok)
  if (N < 2L) stop("need at least two finite observations to rank.", call. = FALSE)

  z <- rep(NA_real_, length(y))
  # (N + 1) in the denominator, not N: it is the Blom-style plotting position that
  # keeps the largest observation off +Inf. With N it would map to qnorm(1).
  z[ok] <- stats::qnorm(rank(y[ok], ties.method = ties.method) / (N + 1))

  # Rebuild the list geometry exactly. split() on a factor of subject indices
  # would reorder by level; this cannot.
  out <- vector("list", length(Ly))
  end <- cumsum(ni); start <- end - ni + 1L
  for (i in seq_along(Ly)) {
    out[[i]] <- if (ni[i] > 0L) z[start[i]:end[i]] else numeric(0)
    names(out[[i]]) <- names(Ly[[i]])
  }
  names(out) <- names(Ly)
  out
}
