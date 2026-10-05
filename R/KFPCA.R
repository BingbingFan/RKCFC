#' @title Kendall Functional Principal Component Analysis (KFPCA) for sparse design
#'
#' @description KFPCA for non-Gaussian functional data with sparse design or longitudinal data.
#'
#' @param Lt A \code{list} of \emph{n} vectors, where \emph{n} is the sample size. Each entry contains the observation time in ascending order for each subject.
#' @param Ly A \code{list} of \emph{n} vectors, where \emph{n} is the sample size. Each entry contains the measurements of each subject at the observation time correspond to \code{Lt}.
#' @param interval A \code{vector} of length two denoting the supporting interval.
#' @param dataType A \code{character} denoting the data type; 'Sparse'-default, 'Dense'.
#' @param nK An integer denoting the number of FPCs.
#' @param kern A \code{character} denoting the kernel type for the Nadaraya-Watson estimators; 'epan'(Epanechnikov)-default, 'unif'(Uniform), 'quar'(Quartic), 'gauss'(Gaussian).
#' @param bw A scalar denoting the bandwidth for the Nadaraya-Watson estimators.
#' @param kernK A \code{character} denoting the kernel type for the estimation of the Kendall's tau function; 'epan'(Epanechnikov)-default, 'unif'(Uniform), 'quar'(Quartic), 'gauss'(Gaussian).
#' @param bwK The bandwidth for the estimation of the Kendall's tau function. If \code{is.numeric(bwK) == T}, \code{bwK} is exactly the bandwidth. If \code{bwK == "GCV"}, the bandwidth is chosen by GCV. (default: "GCV")
#' @param kernmean A \code{character} denoting the kernel type for the estimation of the mean function; 'epan'(Epanechnikov)-default, 'unif'(Uniform), 'quar'(Quartic), 'gauss'(Gaussian).
#' @param bwmean The bandwidth for the estimation of the mean function. If \code{is.numeric(bwmean) == T}, \code{bwmean} is exactly the bandwidth. If \code{bwmean == "GCV"}, the bandwidth is chosen by GCV. (default: "GCV")
#' @param meanrobust A \code{logical}. If \code{FALSE} (the default) the mean function is local-linear least squares. If \code{TRUE} it is a local-linear M-estimator with Huber weights, which is robust to outlying observations in the same way the Kendall's tau covariance is robust to outlying subjects. The choice is stored on the returned object and is inherited by \code{predict.KFPCA}, which re-estimates the mean from the fit's own training data rather than reading \code{$mean}. Default \code{FALSE}, and it must stay \code{FALSE}: published results were produced under it.
#' @param scatter A \code{character} choosing what the covariance estimator signs. \code{"kendall"} (the default) signs PAIRWISE DIFFERENCES between subjects, which is the multivariate Kendall's tau of Zhong et al.\ and requires no location estimate. \code{"sign"} signs each subject's DEVIATION FROM THE MEAN, giving the spherical or spatial-sign covariance; it uses the same mean this fit uses, so that changing this argument alone isolates the scatter estimator. Default \code{"kendall"}, and it must stay so: published results were produced under it.
#' @param nRegGrid An integer denoting the number of equally spaced time points in the supporting interval. The eigenfunctions and mean function are estimated at these equally spaced time points.
#' @param fdParobj A functional parameter object for the smoothing of the eigenfunctions. For more detail, see \code{\link[fda]{smooth.basis}}.
#' @param more Logical; If \code{FALSE}, estimates of FPC scores and predictions of trajectories are not returned.
#'
#' @return A \code{list} containing the following components:
#' \item{ObsGrid}{A \code{vector} containing all observation time points in ascending order.}
#' \item{RegGrid}{A \code{vector} of the equally spaced time points in the support interval.}
#' \item{bwmean}{A scalar denoting the bandwidth for the mean function estimate.}
#' \item{kernmean}{A \code{character} denoting the kernel type for the estimation of the mean function}
#' \item{bwK}{A scalar denoting the bandwidth for the Kendall's tau function estimate.}
#' \item{kernK}{A \code{character} denoting the kernel type for the estimation of the Kendall's tau function}
#' \item{mean}{A \code{vector} of length \code{nRegGrid} denoting the mean function estimate.}
#' \item{KendFun}{A \code{nRegGrid} by \code{nRegGrid} \code{matrix} denoting the Kendall's tau function estimate.}
#' \item{FPC_dis}{A \code{nRegGrid} by \code{nK} \code{matrix} containing the eigenfunction estimates at \code{RegGrid}.}
#' \item{values}{A \code{vector} of length \code{nRegGrid} holding the FULL eigenvalue
#'   spectrum of the estimated Kendall's tau surface, not only the leading \code{nK}.
#'   These are eigenvalues of the \emph{discretised} surface, so they carry a grid
#'   factor that cancels in any ratio but not otherwise. The smoothed surface is not
#'   guaranteed positive semi-definite, so some values are NEGATIVE and must be
#'   truncated at zero (\code{pmax(values, 0)}) before a cumulative fraction of
#'   variance explained is taken. Supplied so that \code{nK} can be selected from the
#'   data on the robust branch, as \code{fdapace::FPCA()$cumFVE} allows on the
#'   non-robust one.}
#' \item{FPC_smooth}{A functional data object for the eigenfunction estimates.}
#' \item{score}{A \emph{n} by \code{nK} \code{matrix} containing the estimates of the FPC scores, where \emph{n} is the sample size. The results are returned when \code{more = TRUE}.}
#' \item{X_fd}{A functional data object for the prediction of trajectories. The results are returned when \code{more = TRUE}.}
#' \item{Xest_ind}{A \code{list} containing the prediction of each trajectory at their own observation time points. The results are returned when \code{more = TRUE}.}
#' \item{Lt}{The input 'Lt'.}
#' \item{Ly}{The input 'Ly'.}
#' \item{CompTime}{A scalar denoting the computation time.}
#' @export
#'
#' @references
#' \cite{Rou Zhong, Shishi Liu, Haocheng Li, Jingxiao Zhang (2021). "Robust Functional Principal Component Analysis for Non-Gaussian Longitudinal Data." Journal of Multivariate Analysis, https://doi.org/10.1016/j.jmva.2021.104864.}
#'
#' @examples
#' # Sparse, irregular curves with two eigenfunctions and heavy-tailed scores.
#' \donttest{
#' set.seed(1)
#' n <- 100
#' interval <- c(0, 10)
#' eigfun <- list(function(x) cos(pi * x / 10) / sqrt(5),
#'                function(x) sin(pi * x / 10) / sqrt(5))
#' score <- cbind(stats::rt(n, df = 3) * 3, stats::rt(n, df = 3) * 1.2)
#' Lt <- lapply(seq_len(n), function(i) sort(stats::runif(sample(6:8, 1), 0, 10)))
#' Ly <- lapply(seq_len(n), function(i)
#'   score[i, 1] * eigfun[[1]](Lt[[i]]) + score[i, 2] * eigfun[[2]](Lt[[i]]) +
#'     stats::rnorm(length(Lt[[i]]), 0, sqrt(0.1)))
#'
#' basis <- fda::create.bspline.basis(interval, nbasis = 13, norder = 4,
#'                                    breaks = seq(0, 10, length.out = 11))
#' Klist <- KFPCA(Lt, Ly, interval, nK = 2, bw = 1, nRegGrid = 51,
#'                fdParobj = basis)
#' plot(Klist$FPC_smooth)
#' }
KFPCA <- function(Lt, Ly, interval, dataType = "Sparse", nK,
                  kern = "epan", bw, kernK = "epan", bwK = "GCV",
                  kernmean = "epan", bwmean = "GCV", meanrobust = FALSE,
                  scatter = c("kendall", "sign"),
                  nRegGrid, fdParobj, more = TRUE){
  scatter <- match.arg(scatter)
  startFPCA <- Sys.time() #Starting time for KFPCA

  # n: sample size
  n <- length(Lt)

  # Xlist: the list of n matrices, where n is the sample size. The i-th matrix contains
  # the observations of the i-th subject and the Nadaraya-Watson estimators of the other
  # subjects at the observation time for the i-th subject
  Xlist <- GetMatList(Lt, Ly, kern, bw)

  # For each matrix in Xlist, drop out the rows that contain NaNs
  naid <- which(unlist(lapply(Xlist, function(x){length(which(is.nan(x) == T))})) != 0)
  for(i in naid){
    narow <- which(apply(Xlist[[i]], 1, function(x){length(which(is.nan(x) == T))}) != 0)
    Xlist[[i]] <- Xlist[[i]][-narow,]
  }
  sub_id <- which(unlist(lapply(Xlist, is.matrix)) == T) 
  #the subjects that can be used to compute the raw Kendall's tau covariances

  # THE GRID AND THE MEAN ARE COMPUTED HERE, ahead of the scatter, because the
  # spatial-sign estimator needs a location to sign deviations from. Both are
  # deterministic functions of (Lt, Ly) and neither touches the RNG, so hoisting
  # them cannot change the default branch -- and the package's fingerprint test
  # asserts exactly that rather than trusting the argument.
  gridobs <- sort(unique(c(unlist(Lt))))
  gridreg <- seq(interval[1], interval[2], length.out = nRegGrid)

  if(is.numeric(bwmean)){
    bw3 <- bwmean
  }else if(bwmean == "GCV"){
    bw3 <- GetGCVbw1D(Lt, Ly, kernmean, dataType = dataType)
  }else{
    stop('bwmean must be a scalar or the character denoting GCV !')
  }
  meanest <- MeanEst(Lt, Ly, kern = kernmean, bw3, gridreg,
                     robust = meanrobust)$mean

  # Lk: the list of n matrices, where n is the sample size. The i-th matrix is the raw
  # covariance contribution obtained from the i-th subject.
  #
  # TWO ESTIMATORS SHARE THIS SLOT, and they differ in WHAT gets signed.
  #
  #   "kendall"  signs PAIRWISE DIFFERENCES. For subject i, every other subject's
  #              Nadaraya-Watson curve is evaluated at i's times and subtracted;
  #              each difference is projected to the unit sphere and averaged.
  #              This is the multivariate Kendall's tau of Zhong et al., it needs
  #              no location estimate at all, and it is the DEFAULT.
  #
  #   "sign"     signs DEVIATIONS FROM A LOCATION -- the spherical / spatial-sign
  #              covariance. One term per subject rather
  #              than n - 1, so it is also far cheaper.
  #
  # The location is the SAME mean this fit uses, `meanest`, so that "sign" against
  # "kendall" isolates the scatter estimator and nothing else. Using a separately
  # robustified location here would confound the two things the comparison exists
  # to separate, which is the mistake the E1/E1M contrast was built to avoid.
  Lk <- NULL
  if (scatter == "kendall") {
    for(i in sub_id){
      mi <- ncol(Xlist[[i]])
      ni <- nrow(Xlist[[i]])
      Lk[[i]] <- matrix(0, ncol = mi, nrow = mi)
      for(j in (2:ni)){
        XX <- Xlist[[i]][1,] - Xlist[[i]][j,]
        Lk[[i]] <- Lk[[i]] + (XX %*% t(XX))/(normq(XX, 2))^2 * length(XX)
      }
      Lk[[i]] <- Lk[[i]]/(ni - 1)
    }
  } else {
    for(i in sub_id){
      # The location at THIS subject's observation times, interpolated from the
      # mean estimated on the regular grid.
      mu_i <- stats::approx(gridreg, meanest, xout = Lt[[i]], rule = 2)$y
      XX <- Xlist[[i]][1,] - mu_i
      nx <- normq(XX, 2)
      # A subject sitting exactly at the location has no direction to contribute.
      # Signing it would divide by zero; contributing a zero matrix is the correct
      # limit and keeps the time-pair bookkeeping below aligned.
      Lk[[i]] <- if (nx > 0) (XX %*% t(XX)) / nx^2 * length(XX) else
                 matrix(0, ncol = length(XX), nrow = length(XX))
    }
  }

  # xintotal: a matrix of the assembled time pairs for the raw Kendall's tau covariance
  xin <- NULL
  for(i in sub_id){
    xin <- rbind(xin, cbind(i, t(utils::combn(Lt[[i]], 2))))
  }
  xinsym <- cbind(xin[,1], xin[,3], xin[,2])
  xintotal <- rbind(xin, xinsym)

  # yin: the raw Kendall's tau covariance for various assembled time pairs
  yin <- NULL
  for(i in sub_id){
    y_id <- which(lower.tri(Lk[[i]]) == T)
    yin <- c(yin, Lk[[i]][y_id])
  }
  yin <- c(yin, yin)

  # Bandwidth for the covariance surface (gridobs and gridreg are set above)
  if(is.numeric(bwK)){
    bw2 <- bwK
  }else if(bwK == "GCV"){
    bw2 <- GetGCVbw2D(xintotal[,2:3], yin, Lt, kernK, gridobs, gridreg, dataType = dataType)
  }else{
    stop('bwK must be a scalar or the character denoting GCV !')
  }

  # Local linear estimate of the Kendall's tau function
  Kend <- fdapace::Lwls2D(bw = bw2, xin = xintotal[,2:3], yin = yin,
                          xout1 = gridreg, xout2 = gridreg)

  # FPC_dis: the first nK eigenfunction estimates at regular grids
  # FPC_smooth: the functional data object for the first nK eigenfunction estimates
  #
  # The eigenvalues were previously computed here and DISCARDED, which left no
  # way to compute a fraction of variance explained on the robust branch and so
  # no way to choose nK from the data -- fdapace returns $lambda and $cumFVE,
  # KFPCA returned neither. They are kept below as $values.
  #
  # NOTE for anyone using them: these are eigenvalues of the DISCRETISED
  # covariance surface, so they carry a grid factor. That factor cancels in any
  # ratio, which is all FVE needs. A smoothed PACE-type surface is also not
  # guaranteed positive semi-definite, so NEGATIVE eigenvalues occur and must be
  # truncated at zero before any cumulative share is taken.
  .eig <- eigen(Kend)
  FPC_dis <- apply(as.matrix(.eig$vectors[,1:nK]), 2, unitize, grid = diff(gridreg)[1])
  FPC_smooth <- fda::smooth.basis(gridreg, FPC_dis, fdParobj)$fd

  # (the grid, the mean bandwidth and the mean itself are computed above, before
  # the scatter step, because the spatial-sign estimator needs a location)

  # The return object
  ret <- list()
  ret$ObsGrid <- gridobs
  ret$RegGrid <- gridreg
  ret$bwmean <- bw3
  ret$kernmean <- kernmean
  # THE MEAN-ESTIMATOR CHOICE TRAVELS ON THE OBJECT, and it has to. On the robust
  # branch `predict.KFPCA()` NEVER READS `object$mean` -- it re-estimates the mean
  # from `object$Lt` and `object$Ly` with its own `MeanEst()` call. An option
  # honoured here but not there would leave prediction centred on the
  # least-squares mean while every returned field claimed otherwise, and every
  # number would have looked plausible. Storing the option on the fit makes that
  # defect structurally impossible rather than merely remembered.
  ret$meanrobust <- meanrobust
  # Which scatter produced KendFun. Recorded for the same reason meanrobust is:
  # a fitted object should say how it was fitted, so nothing downstream has to
  # guess or be told separately.
  ret$scatter <- scatter
  ret$bwK <- bw2
  ret$kernK <- kernK
  ret$mean <- meanest
  ret$KendFun <- Kend
  ret$FPC_dis <- FPC_dis
  ret$values <- .eig$values      # full spectrum; see the note above
  ret$FPC_smooth <- FPC_smooth
  ret$Lt <- Lt
  ret$Ly <- Ly

  class(ret) <- 'KFPCA'

  if(more){
    # FPC scores
    score_ret <- FPCscoreLSE(Lt, Ly, kern = kernmean, bw3, FPC_dis, gridreg, more = TRUE)
    FPC_score <- score_ret$score

    # X_dis: the prediction of trajectories at regular grids
    # X_fd: the functional data object for the prediction of trajectories
    X_dis <- FPC_dis %*% t(FPC_score) + meanest
    X_fd <- fda::smooth.basis(gridreg, X_dis, fdParobj)$fd

    # X_dis_fine: the prediction of trajectories at all observation time points
    # Xest_ind: A list containing the prediction of each trajectory at their own
    # observation time points
    X_dis_fine <- score_ret$FPC_dis_fine %*% t(FPC_score) + score_ret$meanest_fine
    Xest_ind <- list()
    for(i in 1:n){
      Xest_ind[[i]] <- X_dis_fine[match(Lt[[i]], gridobs), i]
    }

    #Computation time
    endFPCA <- Sys.time()
    CompTime <- endFPCA - startFPCA

    ret$score <- FPC_score
    ret$X_fd <- X_fd
    ret$Xest_ind <- Xest_ind
    ret$CompTime <- CompTime

    return(ret)
  }else{
    #Computation time
    endFPCA <- Sys.time()
    CompTime <- endFPCA - startFPCA

    ret$CompTime <- CompTime

    return(ret)
  }

}


