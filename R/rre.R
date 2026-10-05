#' @title Relative reconstruction error
#'
#' @description
#' The residual sum of squares divided by the total sum of squares about the
#' sample mean: 0 is perfect reconstruction, 1 is no better than the mean of the
#' observations, and values above 1 mean the fit reconstructs worse than a
#' constant.
#'
#' @details
#' This replaces the mean squared reconstruction error used in the original
#' submission, for three reasons. It is unit-free, so a simulation on one scale
#' and an application on another are directly comparable, which a raw mean
#' squared error is not. It carries its own reference point, so a single number
#' can be read without a competitor beside it. And it is the native language of
#' functional principal component analysis, where the proportion of variance
#' explained governs component selection, so nothing new has to be defended.
#'
#' The quantity is an \emph{in-sample} reconstruction error whenever
#' \code{yhat} comes from a fit that used \code{y}. For predictive performance,
#' withhold observations from the fit and pass only those here; in simulation,
#' pass the true latent trajectory as \code{y}.
#'
#' \code{trim} additionally reports the error over the middle
#' \code{1 - trim} of observations. If the trimmed and untrimmed values agree,
#' an advantage is broad-based and can be claimed as such; if they diverge, it is
#' concentrated in the tail, which is still a real robustness result but a
#' different claim. \strong{Trimming is on how extreme the observation is, not on
#' how large the error is.} Trimming by error would remove each method's own
#' worst cases and score two methods on different subsets, which is exactly the
#' comparison the trimmed value exists to protect against.
#'
#' @param y Numeric vector of observed values.
#' @param yhat Numeric vector of reconstructed values, the same length as
#'   \code{y}.
#' @param trim Proportion of the most extreme observations to exclude, by
#'   absolute deviation from the median. \code{0} uses all observations.
#'
#' @return A single numeric value.
#'
#' @seealso \code{\link{univar_predicted_yhat}}, \code{\link{robust_normalize}}
#'
#' @examples
#' set.seed(1)
#' y <- stats::rnorm(100)
#' rre(y, yhat = 0.8 * y)          # reconstructs most of the variation
#' rre(y, yhat = rep(mean(y), 100)) # exactly 1: no better than the mean
#'
#' @export
rre <- function(y, yhat, trim = 0) {
  stopifnot(length(y) == length(yhat), trim >= 0, trim < 1)
  i <- seq_along(y)
  if (trim > 0) {
    dev <- abs(y - stats::median(y))
    i <- i[dev <= stats::quantile(dev, 1 - trim, names = FALSE)]
  }
  sum((y[i] - yhat[i])^2) / sum((y[i] - mean(y[i]))^2)
}
