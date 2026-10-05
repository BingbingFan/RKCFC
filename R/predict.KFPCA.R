#' @title Predict FPC scores
#'
#' @description Predict FPC scores using least square estimate (LSE) for a new sample.
#'
#' @param object A KFPCA object obtained from \code{\link[KFPCA]{KFPCA}}.
#' @param newLt A \code{list} of \emph{n} vectors, where \emph{n} is the new sample size. Each entry contains the observation time in ascending order for each new subject.
#' @param newLy A \code{list} of \emph{n} vectors, where \emph{n} is the new sample size. Each entry contains the measurements of each new subject at the observation time correspond to \code{newLt}.
#' @param nK An integer denoting the number of FPCs.
#' @param more Logical; If \code{FALSE}, only the predictions of FPC scores are returned. If \code{TRUE}, the mean function estimates and the eigenfunction estimates at the new observation time points are also returned.
#' @param ... Not used.
#'
#' @return If \code{more = FALSE}, a \emph{n} by \code{nK} \code{matrix} containing the predictions of the FPC scores is returned, where \emph{n} is the new sample size. If \code{more = TRUE}, a \code{list} containing the following components is returned:
#' \item{score_new}{a \emph{n} by \code{nK} \code{matrix} containing the predictions of the FPC scores.}
#' \item{meanest_new}{Mean function estimates at the new observation time points.}
#' \item{FPC_dis_new}{Eigenfunction estimates at the new observation time points.}
#'
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' interval <- c(0, 10)
#' eigfun <- list(function(x) cos(pi * x / 10) / sqrt(5),
#'                function(x) sin(pi * x / 10) / sqrt(5))
#' sim <- function(n) {
#'   score <- cbind(stats::rnorm(n, 0, 3), stats::rnorm(n, 0, sqrt(1.5)))
#'   Lt <- lapply(seq_len(n), function(i) sort(stats::runif(sample(6:8, 1), 0, 10)))
#'   Ly <- lapply(seq_len(n), function(i)
#'     score[i, 1] * eigfun[[1]](Lt[[i]]) + score[i, 2] * eigfun[[2]](Lt[[i]]) +
#'       stats::rnorm(length(Lt[[i]]), 0, sqrt(0.1)))
#'   list(Lt = Lt, Ly = Ly, score = score)
#' }
#' train <- sim(100)
#' test  <- sim(20)
#'
#' basis <- fda::create.bspline.basis(interval, nbasis = 13, norder = 4,
#'                                    breaks = seq(0, 10, length.out = 11))
#' Klist <- KFPCA(train$Lt, train$Ly, interval, nK = 2, bw = 1,
#'                nRegGrid = 51, fdParobj = basis)
#'
#' score_pre <- predict(Klist, test$Lt, test$Ly, nK = 2)
#' plot(test$score[, 1], score_pre[, 1])
#' }
predict.KFPCA <- function(object, newLt, newLy, nK, more = FALSE, ...){

  Lt <- object$Lt
  Ly <- object$Ly
  RegGrid <- object$RegGrid
  FPC_dis <- object$FPC_dis
  bwmean <- object$bwmean
  kernmean <- object$kernmean

  # ObsGridnew: observation grids of newLt
  # nObsGridnew: the number of new observation grids
  # FPC_dis_new: a nObsGridnew by nK matrix containing eigenfunction estimates at various
  # observation grids
  ObsGridnew <- sort(unique(unlist(newLt)))
  nObsGridnew <- length(ObsGridnew)
  bwFPC <- NULL
  FPC_dis_new <- matrix(0, ncol = nK, nrow = nObsGridnew)
  for(i in 1:nK){
    bwFPC[i] <- GetGCVbw1D(Lt = list(RegGrid), Ly = list(FPC_dis[,i]), kern = "epan",
                           dataType = "Dense")
    FPC_dis_new[,i] <- fdapace::Lwls1D(bwFPC[i], kernel_type = "epan", xin = RegGrid, yin = FPC_dis[,i],
                                       xout = ObsGridnew)
  }

  # meanest_new: mean function estimates at the new observation time points
  # Inherit the fit's own mean estimator. `%||%`-style fallback because objects
  # fitted before this field existed have no `$meanrobust`, and those must keep
  # behaving exactly as they did -- least squares.
  meanrobust <- isTRUE(object$meanrobust)
  meanest_new <- MeanEst(Lt, Ly, kern = kernmean, bw = bwmean,
                         gridout = ObsGridnew, robust = meanrobust)$mean

  # mean_ind: a list of vectors. The i-th vector contains the mean estimates at the observation
  # time of the i-th new subject.
  # phi: a list of matrices. The i-th matrix contains the eigenfunction estimates at the
  # observation time of the i-th new subject.
  # score_new: a n by nK matrix containing the estimates of the FPC scores. The (i, j)-th element is
  # the j-th FPC score estimate of the i-th new subject.
  n_new <- length(newLt)
  mean_ind <- list()
  phi <- list()
  score_new <- matrix(0, ncol = nK, nrow = n_new)
  for(i in 1:n_new){
    id <- sapply(newLt[[i]], function(x){which(ObsGridnew == x)})
    mean_ind[[i]] <- meanest_new[id]
    phi[[i]] <- FPC_dis_new[id,]
    response <- newLy[[i]] - mean_ind[[i]]
    score_new[i,] <- as.vector(.rkcfc_solve_covariance(
      crossprod(phi[[i]]),
      crossprod(phi[[i]], response)
    ))
  }

  if(more == FALSE){
    return(score_new)
  }else{
    ret <- list()
    ret$score_new <- score_new
    ret$meanest_new <- meanest_new
    ret$FPC_dis_new <- FPC_dis_new
    return(ret)
  }

}
