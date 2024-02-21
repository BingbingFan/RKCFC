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
  meanest_new <- MeanEst(Lt, Ly, kern = kernmean, bw = bwmean, gridout = ObsGridnew)$mean

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
    score_new[i,] <- MASS::ginv(t(phi[[i]]) %*% phi[[i]]) %*% t(phi[[i]]) %*% (newLy[[i]] - mean_ind[[i]])
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
