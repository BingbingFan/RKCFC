predicted_yhat <- function(FPCA_k, ith_yy,  target_t, robust, optns ){
  # suppressMessages(require(KFPCA))
  J <- length(FPCA_k)
  ## i th yy : curve that I want to predict
  if(robust){

    n_M <- lapply(1:length(ith_yy), function(x) dim(FPCA_k[[x]]$score)[2] )
    phi <- lapply(1:length(ith_yy), function(x) FPCA_k[[x]]$FPC_smooth )
    mu <- lapply(1:length(ith_yy), function(x) FPCA_k[[x]]$mean )
    workGrid <- lapply(1:length(ith_yy), function(x) FPCA_k[[x]]$RegGrid )

    prd <- lapply(1:J, function(x) predict.KFPCA(FPCA_k[[x]], newLt = list(target_t[[x]]), newLy = list(ith_yy[[x]]),
                                                 nK = optns$nK, more = TRUE))
    xi_Est_i <- lapply(1:J, function(x) prd[[x]]$score_new)
    pi_temp <- lapply(1:J, function(x) prd[[x]]$FPC_dis_new)
    # mfpca
    xiEst <- lapply(1:J, function(x) FPCA_k[[x]]$score)

    fittedX <- multi_eigen <- multi_score <- list()
    full_xiEst <- vector()
    for(j in 1:length(ith_yy) ){
      full_xiEst <- cbind(full_xiEst, xiEst[[j]])
    }

    Z <-  t(full_xiEst) %*% full_xiEst / (nrow(full_xiEst)-1)
    index <- 0
    for(j in 1: J   ){
      muxtemp <- spline(workGrid[[j]], mu[[j]], xout = target_t[[j]])$y
      index <- seq(from = (max(index)+1 ) , length.out = n_M[[j]])
      multi_eigen_j <- (svd(Z)$u)[index,index] %*% t(pi_temp[[j]])
      multi_score_j <- xi_Est_i[[j]]  %*% (svd(Z)$u)[index,index]
      multi_eigen[[j]]  <-  multi_eigen_j
      multi_score[[j]]  <-  multi_score_j
      fittedX[[j]]  <-  as.numeric(muxtemp + (t(multi_eigen_j) %*% t(multi_score_j)))
    }


  }else{
    n_M <- lapply(1:length(ith_yy), function(x) length(FPCA_k[[x]]$lambda) )
    lambda <- lapply(1:length(ith_yy), function(x) FPCA_k[[x]]$lambda )
    phi <- lapply(1:length(ith_yy), function(x) FPCA_k[[x]]$phi )
    mu <- lapply(1:length(ith_yy), function(x) FPCA_k[[x]]$mu )
    workGrid <- lapply(1:length(ith_yy), function(x) FPCA_k[[x]]$workGrid )
    sig2 <- lapply(1:length(ith_yy), function(x) FPCA_k[[x]]$sigma2 )

    xiEst <- fittedX <- xi_Est_i <-  list()
    multi_eigen <- multi_score <- list()

    for(j in 1:length(ith_yy) ){

      pi_temp <- matrix(nrow = length(ith_yy[[j]]), ncol = ncol( phi[[j]]) )
      for(l in 1:ncol(phi[[j]])){
        pi_temp[,l] <- spline(workGrid[[j]] ,phi[[j]][, l], xout = target_t[[j]])$y
      }
      mu_ij <- spline(workGrid[[j]], mu[[j]], xout = target_t[[j]])$y
      sigma_i <- matrix(nrow = length(ith_yy[[j]]), ncol = length(ith_yy[[j]])  )
      for (mi in 1:length(ith_yy[[j]]) ) {
        for (mj in 1:length(ith_yy[[j]])) {
          if(mi >= mj){
            sigma_ij <- (ith_yy[[j]][mi] - mu_ij[mi])*(ith_yy[[j]][mj] - mu_ij[mj])
            if(mi == mj){
              sigma_ij <- sigma_ij + sig2[[j]]
              sigma_i[mi, mj] <- sigma_ij
            }
            sigma_i[mi, mj] <- sigma_ij
            sigma_i[mj, mi] <- sigma_ij
          }


        }
      }
      xi_Est_i[[j]] <- unlist(lapply(1:n_M[[j]], function(x)
        lambda[[j]][x] %*% t(pi_temp)[x,] %*% MASS::ginv(sigma_i) %*% (ith_yy[[j]] - mu_ij)))
      xiEst[[j]]  <-  FPCA_k[[j]]$xiEst
    }

    full_xiEst <- vector()
    for(j in 1:length(ith_yy) ){
      full_xiEst <- cbind(full_xiEst, xiEst[[j]])
    }

    Z <-  t(full_xiEst) %*% full_xiEst / (nrow(full_xiEst)-1)

    index <- 0
    for(j in 1: length(ith_yy)   ){
      muxtemp <- spline(workGrid[[j]], mu[[j]], xout = target_t[[j]])$y
      index <- seq(from = (max(index)+1 ) , length.out = n_M[[j]])

      pi_temp <- matrix(nrow=length(muxtemp), ncol=ncol(phi[[j]]) )

      for(l in 1:ncol(phi[[j]])){
        pi_temp[,l] <- spline(workGrid[[j]] ,phi[[j]][, l], xout=target_t[[j]])$y
      }
      multi_eigen_j <- (svd(Z)$u)[index,index] %*% t(pi_temp)
      multi_score_j <- t(xi_Est_i[[j]] ) %*% (svd(Z)$u)[index,index]
      multi_eigen[[j]]  <-  multi_eigen_j
      multi_score[[j]]  <-  multi_score_j
      fittedX[[j]]  <-  as.numeric(muxtemp + (t(multi_eigen_j) %*% t(multi_score_j)))
    }

  }

  return(list(fittedX=fittedX, multi_score = multi_score, multi_eigen = multi_eigen ))
}
