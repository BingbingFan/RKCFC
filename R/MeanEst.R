#' @title Local linear estimates of mean function
#'
#' @description Local linear estimates of mean function.
#'
#' @param Lt A \code{list} of \emph{n} vectors, where \emph{n} is the sample size. Each entry contains the observation time in ascending order for each subject.
#' @param Ly A \code{list} of \emph{n} vectors, where \emph{n} is the sample size. Each entry contains the measurements of each subject at the observation time correspond to \code{Lt}.
#' @param kern A \code{character} denoting the kernel type; 'epan'(Epanechnikov), 'unif'(Uniform), 'quar'(Quartic), 'gauss'(Gaussian).
#' @param bw A scalar denoting the bandwidth.
#' @param gridout A \code{vector} denoting the time points that the mean function need to be estimated.
#' @param robust A \code{logical}. If \code{FALSE} (the default) the mean is local-linear least squares, exactly as before. If \code{TRUE} it is a local-linear M-estimator fitted by iteratively reweighted least squares with Huber weights. The default MUST stay \code{FALSE}: published results were produced under it.
#' @param tuning A scalar Huber tuning constant, used only when \code{robust = TRUE}. The default 1.345 gives 95 percent asymptotic efficiency at the Gaussian.
#' @param maxit Maximum IRLS iterations, used only when \code{robust = TRUE}.
#' @param tol Relative convergence tolerance on the fitted curve, used only when \code{robust = TRUE}.
#'
#' @return A \code{list} containing the following components:
#' \item{Grid}{A \code{vector} denoting the time points that the mean function need to be estimated.}
#' \item{mean}{A \code{vector} containing the mean function estimates.}
#' @noRd
#' @examples
#' # Generate data
#' n <- 100
#' interval <- c(0, 10)
#' lambda_1 <- 9 #the first eigenvalue
#' lambda_2 <- 1.5 #the second eigenvalue
#' eigfun <- list()
#' eigfun[[1]] <- function(x){cos(pi * x/10)/sqrt(5)}
#' eigfun[[2]] <- function(x){sin(pi * x/10)/sqrt(5)}
#' score <- cbind(rnorm(n, 0, sqrt(lambda_1)), rnorm(n, 0, sqrt(lambda_2)))
#' DataNew <- GenDataKL(n, interval = interval, sparse = 6:8, regular = FALSE,
#'                      meanfun = function(x){x}, score = score,
#'                      eigfun = eigfun, sd = sqrt(0.1))
#' # Mean function estimate at all observation time points
#' bwOpt <- GetGCVbw1D(DataNew$Lt, DataNew$Ly, kern = "epan")
#' meanest <- MeanEst(DataNew$Lt, DataNew$Ly, kern = "epan", bw = bwOpt,
#'                    gridout = sort(unique(unlist(DataNew$Lt))))
#' plot(meanest$Grid, meanest$mean)
MeanEst <- function(Lt, Ly, kern, bw, gridout, robust = FALSE,
                    tuning = 1.345, maxit = 20L, tol = 1e-7){
  xinraw <- unlist(Lt)
  yinraw <- unlist(Ly)

  #time in ascending order
  xin <- sort(xinraw)
  yin <- yinraw[order(xinraw)]

  if (!isTRUE(robust)) {
    # THE DEFAULT PATH IS THE ORIGINAL LINE, UNCHANGED AND UNCONDITIONAL. It is
    # written out rather than routed through the weighted call below with unit
    # weights, because `Lwls1D`'s own default is `rep(1L, ...)` -- integer, not
    # double -- and reproducing 537.9 core-h of published results is not the
    # place to rely on that coercion being immaterial.
    meanest <- fdapace::Lwls1D(bw, kern, xin = xin, yin = yin, xout = gridout)
  } else {
    # ROBUST LOCAL-LINEAR M-ESTIMATOR, by iteratively reweighted
    # local linear least squares with Huber weights.
    #
    # WHY THIS EXISTS. Kendall's tau protects the estimator of the class-specific
    # EIGENSPACE and does nothing whatever for the mean, because the mean is
    # local-linear least squares. Measured in the estimation block 2026-09-07:
    # on the t(1.5) arm the eigenfunction angles hold at ~6 degrees for both the
    # robust and the non-robust estimator, while `mean_err` degrades IDENTICALLY
    # for both -- 0.609 against 0.652, indistinguishable, from ~0.24 on the
    # Gaussian arm. That gap is what this closes.
    #
    # Huber rather than a redescending psi: a redescending weight can reject a
    # whole local window at once on sparse data, where a window may hold only a
    # handful of points, and the estimate then has no support at all. Huber
    # downweights without ever zeroing.
    ux    <- unique(xin)
    idx   <- match(xin, ux)
    w     <- rep(1, length(xin))
    fit_u <- fdapace::Lwls1D(bw, kern, win = w, xin = xin, yin = yin, xout = ux)
    for (it in seq_len(maxit)) {
      r <- yin - fit_u[idx]
      # MAD, normal-consistent. If it collapses to zero the residuals are
      # degenerate and there is nothing to downweight; stop rather than divide.
      s <- stats::median(abs(r - stats::median(r))) / 0.6744897501960817
      if (!is.finite(s) || s <= 0) break
      w    <- pmin(1, tuning / pmax(abs(r / s), .Machine$double.eps))
      prev <- fit_u
      fit_u <- fdapace::Lwls1D(bw, kern, win = w, xin = xin, yin = yin, xout = ux)
      if (max(abs(fit_u - prev)) < tol * max(1, max(abs(prev)))) break
    }
    meanest <- fdapace::Lwls1D(bw, kern, win = w, xin = xin, yin = yin,
                               xout = gridout)
  }

  ret <- list()
  ret$Grid <- gridout
  ret$mean <- meanest

  return(ret)
}
