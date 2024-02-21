#' @export
univar_predicted_yhat <- function(FPCA_k, ith_yy,  target_t, robust, optns ){

  ## i th yy : curve that I want to predict
  if(robust){
    prd <- predict.KFPCA(FPCA_k, newLt = list(target_t), newLy = list(ith_yy), nK = optns$nK, more = TRUE)
    xi_Est_i <- prd$score_new
    pi_temp <- prd$FPC_dis_new
    fittedX  <- as.numeric(prd$meanest_new + xi_Est_i %*% t(pi_temp) )

  }else{
    n_M <- length(FPCA_k$lambda)
    lambda <- FPCA_k$lambda
    phi <- FPCA_k$phi
    mu <- FPCA_k$mu
    workGrid <- FPCA_k$workGrid
    sig2 <- FPCA_k$sigma2

    pi_temp <- matrix(nrow = length(ith_yy), ncol = ncol( phi) )
    for(l in 1:ncol(phi)){
      pi_temp[,l] <- spline(workGrid ,phi[, l], xout = target_t)$y
    }
    mu_ij <- spline(workGrid, mu, xout = target_t)$y
    sigma_i <- matrix(nrow = length(ith_yy), ncol = length(ith_yy)  )
    for (mi in 1:length(ith_yy) ) {
      for (mj in 1:length(ith_yy)) {
        if(mi >= mj){
          sigma_ij <- (ith_yy[mi] - mu_ij[mi])*(ith_yy[mj] - mu_ij[mj])
          if(mi == mj){
            sigma_ij <- sigma_ij + sig2
            sigma_i[mi, mj] <- sigma_ij
          }
          sigma_i[mi, mj] <- sigma_ij
          sigma_i[mj, mi] <- sigma_ij
        }
      }
    }
    xi_Est_i <- unlist(lapply(1:n_M, function(x)
      lambda[x] %*% t(pi_temp)[x,] %*% MASS::ginv(sigma_i) %*% (ith_yy - mu_ij)))

    fittedX  <- as.numeric(mu_ij + t(xi_Est_i) %*% t(pi_temp))
  }

  return(list(fittedX=fittedX, xi_Est_i = xi_Est_i, pi_i = pi_temp ))
}
