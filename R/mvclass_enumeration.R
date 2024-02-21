#' @export

# CH ----
CH_mv_dt <- function(dtlist, mod_kcfc, mod_rkcfc){
  J <- length(dtlist)
  Ni <- length(unique(dtlist[[1]]$id))
  g1 <- data.frame(id = 1:Ni, g1 = mod_kcfc$res_clust$g)
  g2 <- data.frame(id = 1:Ni, g2 = mod_rkcfc$res_clust$g)
  for (ii in 1:J) {
    dtlist[[ii]] <- merge(dtlist[[ii]],g1, by="id", all.x = T )
    dtlist[[ii]] <- merge(dtlist[[ii]],g2, by="id", all.x = T )
  }
  return(dtlist)
}
CH_mv_dt1 <- function(dtlist, mod_kcfc, mod_rkcfc){
  J <- length(dtlist)
  Ni <- length(unique(dtlist[[1]]$id))
  g1 <- data.frame(id = 1:Ni, g1 = rep(1, Ni))
  g2 <- data.frame(id = 1:Ni, g2 = rep(1, Ni))
  for (ii in 1:J) {
    dtlist[[ii]] <- merge(dtlist[[ii]],g1, by="id", all.x = T )
    dtlist[[ii]] <- merge(dtlist[[ii]],g2, by="id", all.x = T )
  }
  return(dtlist)
}

CH_mv <- function(dtlist, class, id, y, time, option,interval, robust,  basis){

  CH_sparse <- function(data, class, id, y, time, basis, robust = TRUE,
                        option = NULL,
                        interval ){
    # require(KFPCA)
    # require(fdapace)
    # require(dplyr)
    data <- data[,c( id,class, time, y)]
    names(data) <- c( "id","g", "time", "y")
    K <- length(unique(data$g))
    n <- length(unique(data$id))
    # FPCA data
    dtk <- lapply(1:K, function(x){
      subset(data, g == x)
    })

    fdata <- lapply(1:K, function(x) fdapace::MakeFPCAInputs(IDs = dtk[[x]]$id, yVec = dtk[[x]]$y, tVec = dtk[[x]]$time ))

    # estimate mean function using FPCA or KFPCA
    mu_func <- function(dt, robust, option){
      if(robust){
        FPCA_obj <- KFPCA(Lt = dt$Lt,  Ly = dt$Ly, interval = option$interval, nK = option$nK,
                          bw = option$bw, bwK = option$bwK,
                          nRegGrid = option$nRegGrid, fdParobj = option$basis)
        mu <- FPCA_obj$mean
        Grid <- FPCA_obj$RegGrid
      }else{
        FPCA_obj <- fdapace::FPCA(Ly = dt$Ly, Lt = dt$Lt, optns = option)
        mu <- FPCA_obj$mu
        Grid <- FPCA_obj$workGrid
      }
      return(list(mu = mu, Grid = Grid))
    }
    # dt=fdata[[1]]
    FPCA_res <- lapply(1:K, function(x) mu_func(dt=fdata[[x]], robust, option = option) )

    mus <- lapply(1:K, function(x) FPCA_res[[x]]$mu )
    Grids <- lapply(1:K, function(x) FPCA_res[[x]]$Grid )
    # within-cluster mean sum of square
    Ws <- numeric(K)
    for (k in 1:K) {
      mean_fd <- fda::smooth.basis(Grids[[k]], mus[[k]], basis)$fd
      mean_i <- lapply(fdata[[k]]$Lt, function(x) fda:::eval.fd(x,mean_fd ))
      dtk[[k]]$mu <- unlist(mean_i)
      dtk[[k]]$d2 <- (dtk[[k]]$y - dtk[[k]]$mu)^2
      dtk[[k]] <- dtk[[k]] %>%
        group_by(id) %>%
        mutate(w_i = mean(d2, na.rm=T))
      w_i_k <- data.frame(unique(dtk[[k]][, c("id","w_i")]))
      Ws[k] <- mean(w_i_k$w_i, na.rm=T)
    }

    W_k<- mean(Ws)
    # bewteen-cluster distance
    dt <- fdapace::MakeFPCAInputs(IDs = data$id, yVec = data$y, tVec = data$time )
    if(robust){
      FPCA_obj <- KFPCA(Lt = dt$Lt,  Ly = dt$Ly, interval = option$interval, nK = option$nK,
                        bw = option$bw, bwK = option$bwK,
                        nRegGrid = option$nRegGrid, fdParobj = option$basis)


      mu <- FPCA_obj$mean
      Grid <- FPCA_obj$RegGrid
    }else{
      FPCA_obj <- fdapace::FPCA(Ly = dt$Ly, Lt = dt$Lt, optns = option)
      mu <- FPCA_obj$mu
      Grid <- FPCA_obj$workGrid
    }
    bss <- numeric(K)
    for (k in 1:K) {
      nk <- length(unique(dtk[[k]]$id))
      bss[k] <- nk * sum((mus[[k]]-mu)^2)
    }
    B_K <- sum(bss)/n
    # CH
    CH <- B_K / W_k
    return(list(CH = CH, W_K =W_k, Ws = Ws , B_K = B_K, bss_by_n = bss/n ))
  }
  nY <- length(dtlist)

  CH_Uni <- lapply(1:nY, function(x)  CH_sparse(data = dtlist[[x]], class = class ,id = id,y = y[x],
                                                time = time[x] , basis = basis, robust = robust,
                                                option = option ,interval = interval))

  multi_W <- multi_B <- NULL
  for (i in 1:nY) {
    multi_W <- c(multi_W, CH_Uni[[i]]$W_K)
    multi_B <- c(multi_B, CH_Uni[[i]]$B_K)
  }
  CH <- mean(multi_B)/mean(multi_W)

  return(list(CH = CH , W=mean(multi_W), B = mean(multi_B),
              multi_W = multi_W, multi_B = multi_B, CH_Uni = CH_Uni ))

}
# sil ----
sil_score_mv <- function(data, class, id, y, time,  robust = FALSE, option = NULL){
  silhouette_score <- function(data, class, id, y, time,  robust = FALSE, option = NULL ){
    # require(fdapace)
    # if(robust){
    #   require(KFPCA)
    # }
    # require(dplyr)
    data <- data[,c( id,class, time, y)]
    names(data) <- c( "id","g", "time", "y")
    K <- length(unique(data$g))
    n <- length(unique(data$id))
    # FPCA data
    dtk <- lapply(1:K, function(x){
      subset(data, g == x)
    })

    fdata_k <- lapply(1:K, function(x) MakeFPCAInputs(IDs = dtk[[x]]$id, yVec = dtk[[x]]$y, tVec = dtk[[x]]$time ))
    # estimate mean function using FPCA or KFPCA
    mu_func <- function(dt,robust, option){
      if(robust){
        FPCA_obj <- KFPCA(Lt = dt$Lt,  Ly = dt$Ly,  nK = option$nK, bw = option$bw, bwK =option$bwK,
                          nRegGrid = option$nRegGrid, fdParobj = option$basis , interval =  option$interval)
        mu <- FPCA_obj$mean
        Grid <- FPCA_obj$RegGrid
      }else{
        FPCA_obj <- fdapace::FPCA(Ly = dt$Ly, Lt = dt$Lt,
                                  option)
        mu <- FPCA_obj$mu
        Grid <- FPCA_obj$workGrid
      }
      return(list(mu = mu, Grid = Grid))
    }

    FPCA_res <- lapply(1:K, function(x) mu_func(dt = fdata_k[[x]] ,robust, option))
    mus <- lapply(1:K, function(x) FPCA_res[[x]]$mu )
    Grids <- lapply(1:K, function(x) FPCA_res[[x]]$Grid )
    # within-cluster mean sum of square
    # for subjects i , calculate L2 distance of i with u_k(t)
    fdata <- MakeFPCAInputs(IDs = data$id, yVec = data$y, tVec = data$time )

    L2_dis <- function(Lyi, Lti, mus, Grids){
      mean_i_k <-  lapply(1:K, function(x) spline(Grids[[x]] ,mus[[x]], xout = Lti)$y)
      d_L2 <- lapply(1:K, function(x) (Lyi - mean_i_k[[x]])^2 )
      d_L2 <- unlist(lapply(1:K, function(x) mean(d_L2[[x]])))
    }

    L2_i <- lapply(1:n, function(x) L2_dis(Lyi = fdata$Ly[[x]], Lti = fdata$Lt[[x]], mus, Grids) )
    #   cal calculate b_i and a_i
    uid <- data.frame(id = unlist(fdata$Lid))
    class_i <- data.frame(unique(data[, c("id","g")]))
    class_i <- merge(uid, class_i,  by="id", all.x = T )
    # bi ai fun
    a_i <- unlist(lapply(1:n, function(x) L2_i[[x]][ class_i$g[x] ]  ))
    b_i <- unlist(lapply(1:n, function(x) min(L2_i[[x]][ - class_i$g[x] ])  ))

    silhouette <- mean( unlist( lapply(1:n, function(x) (b_i[x] - a_i[x])/max(b_i[x], a_i[x]) ) ) )

    return(list(a_i = a_i, b_i = b_i , silhouette = silhouette, n = n))
  }

  J <- length(data)
  sil_mv <- lapply(1:J, function(x) silhouette_score(data = data[[x]], class , id , y=y[x] ,
                                                     time=time[x] ,  robust , option ))
  sil_score <- mean(unlist(lapply(1:J, function(x) sil_mv[[x]]$silhouette)))
  return(list(silhouette_score = sil_score, sil_mv = sil_mv ))
}

# ICs ----
pIC_mv <- function(data,id,time,y, clust, robust = F, option ){

  # data ----
  J <- length(y)
  K <- length(unique(clust))
  Ldt <- lapply(1:J, function(x) fdapace::MakeFPCAInputs(IDs =  data[[x]][[id]],
                                                         tVec =  data[[x]][[time[x]]],
                                                         yVec =  data[[x]][[y[x]]]))
  # hatY function----
  MRFPCA_hatY <- function(Ldt, clust, robust, option){
    # cluster number: K and number of Y :J
    K <-  length(unique(clust))
    J <- length(Ldt)
    clustids <- lapply(1:K, function(x) which(clust == x)  )
    Lid_list <- list()
    Lt_list <- list()
    Ly_list <- list()
    for (i in 1:K) {
      Lid_list[[i]] <- list()
      Lt_list[[i]] <- list()
      Ly_list[[i]] <- list()
      idk <- clustids[[i]]
      for (j in 1:J) {
        Lid_list[[i]][[j]] <- lapply( idk,  function(x) Ldt[[j]]$Lid[[x]])
        Lt_list[[i]][[j]] <- lapply( idk,  function(x) Ldt[[j]]$Lt[[x]])
        Ly_list[[i]][[j]] <- lapply( idk,  function(x) Ldt[[j]]$Ly[[x]])
      }
    }
    #  MFPCA fuctions ----
    RMFPCA <- function(Ly, Lt, robust, option){
      # library(fdapace)
      # library(KFPCA)
      J <- length(Ly)
      n <- length(Ly[[1]])
      if(robust){
        uFPCA <- lapply(1:J, function(x) KFPCA(Lt = Lt[[x]],  Ly = Ly[[x]], interval = option$interval,
                                               nK = option$nK, bw = option$bw,
                                               nRegGrid = option$nRegGrid, fdParobj = option$basis) )
        FPC_score <- lapply(1:J, function(x) uFPCA[[x]]$score)
        FPC <- lapply(1:J, function(x) uFPCA[[x]]$FPC_dis)
        mu <- lapply(1:J, function(x) uFPCA[[x]]$mean)
        workgrid <- lapply(1:J, function(x) uFPCA[[x]]$RegGrid)
      }else{
        uFPCA <- lapply(1:J, function(x) FPCA(Ly = Ly[[x]], Lt = Lt[[x]], optns = option))
        FPC_score <- lapply(1:J, function(x) uFPCA[[x]]$xiEst)
        FPC <- lapply(1:J, function(x) uFPCA[[x]]$phi)
        mu <- lapply(1:J, function(x) uFPCA[[x]]$mu)
        workgrid <- lapply(1:J, function(x) uFPCA[[x]]$workGrid)
      }
      n_M <- lapply(1:J, function(x) ncol(FPC_score[[x]]))
      # Z matrix
      full_xiEst <- vector()
      for(j in 1:J ){
        full_xiEst <- cbind(full_xiEst, FPC_score[[j]])
      }
      Z <-  t(full_xiEst) %*% full_xiEst / (nrow(full_xiEst)-1)

      index_M <- list()
      index <- 0
      for (j in 1:J) {
        index <- seq(from = (max(index)+1 ) , length.out = n_M[[j]])
        index_M[[j]] <- index
      }

      multi_score <- NULL
      multi_eigen <- fittedX <- list()
      for (j in 1:J) {
        index <- index_M[[j]]
        # multi score
        multi_score_j <- t(svd(Z)$u[index,index] %*% t(FPC_score[[j]]))
        multi_score <- cbind(multi_score,multi_score_j)
        # multi_eigen_j
        multi_eigen_j <- t(svd(Z)$u[index,index] %*% t(FPC[[j]]))
        multi_eigen[[j]] <- multi_eigen_j
        # hatY at workgrid
        fittedX[[j]]  <-  mu[[j]] + (multi_score_j %*% t(multi_eigen_j))
      }
      # hatY at tij
      hatY_fun <- function(j, i, Lt, mu, multi_eigen, multi_score, workgrid,  index_M){
        index <- index_M[[j]]
        tij <- Lt[[j]][[i]]
        mscore_ij <- multi_score[i,index]
        muj <- mu[[j]]
        multi_eigen_j <- multi_eigen[[j]]
        workgrid_j <- workgrid[[j]]
        muxtemp <- spline(workgrid_j, muj, xout = tij)$y
        phi_temp <- matrix(nrow=length(muxtemp), ncol=ncol(multi_eigen_j) )
        for(l in 1:ncol(multi_eigen_j)){
          phi_temp[,l] <- spline(workgrid_j ,multi_eigen_j[, l], xout = tij)$y
        }
        multi_eigen_j_i <- svd(Z)$u[index,index] %*% t(phi_temp)
        fittedX_j_i  <-  as.numeric(muxtemp + (mscore_ij %*% multi_eigen_j_i))
        return(fittedX_j_i)
      }

      hatY_ob <- list()
      for (j in 1:J) {

        hatY_ob[[j]] <- lapply(1:n, function(x)
          hatY_fun(j = j, i = x, Lt, mu, multi_eigen, multi_score, workgrid,  index_M))
      }


      return(list(mu = mu, uFPC_score = FPC_score, uFPC = FPC,
                  multi_score = multi_score, multi_FPC = multi_eigen,
                  multi_hatY_full = fittedX, multi_hatY_ob = hatY_ob , uFPCA = uFPCA,
                  data = list(Lid_list = Lid_list, Lt_list = Lt_list,Ly_list = Ly_list)))

    }
    # MFPCA  by class
    hatY_by_k <- lapply(1:K, function(x) RMFPCA(Ly = Ly_list[[x]] , Lt = Lt_list[[x]], robust = F, option = option))
    return(hatY_by_k)
  }

  hatY <- MRFPCA_hatY(Ldt, clust, robust = F, option)
  hatY_ob <- Ly_ob <- sigmas <- M <- list()
  for (j in 1:J) {
    hatY_ob[[j]] <- Ly_ob[[j]] <- list()
    sigmas[[j]] <- M[[j]] <- 0
    for (k in 1:K) {
      hatY_ob[[j]] <- c(hatY_ob[[j]], hatY[[k]]$multi_hatY_ob[[j]])
      Ly_ob[[j]] <- c(Ly_ob[[j]], hatY[[1]]$data$Ly_list[[k]][[j]])
      sigmas[[j]] <- c(sigmas[[j]], hatY[[k]]$uFPCA[[j]]$sigma2)
      M[[j]] <- c(M[[j]], dim(hatY[[k]]$uFPCA[[j]]$xiEst)[2] )
    }
    sigmas[[j]] <- sigmas[[j]][-1]
    M[[j]] <- M[[j]][-1]
  }
  # clustids <- lapply(1:K, function(x) which(clust == x)  )
  ns <- as.numeric(table(clust))
  class_new <- NULL
  for (k in 1:K) {
    class_new <- c(class_new, rep(k, ns[k]))
  }
  # pICs----
  pseudoIC<- function(Ly, hatY, sigma_k, class, M, K ){
    n <- length(Ly)
    N <- length(unlist(Ly))
    ni <- unlist(lapply(1:n ,function(x) length(Ly[[x]])))
    sigma_i <- sigma_k[class]
    ss <- unlist(lapply(1:n ,function(x)  sum( Ly[[x]] -  hatY[[x]] )^2 )  )
    hat_LL <- sum(-0.5*log(2*pi)*ni - 0.5*ni*log(sigma_i) -0.5/sigma_i* ss )
    nu <- sum(M)
    pAIC <- -2*hat_LL + 2*(K+nu)
    pBIC <- -2*hat_LL + log(N)*(K+nu)
    return(list(pAIC = pAIC, pBIC = pBIC, hat_LL = hat_LL))
  }
  #
  uICs <- lapply(1:J, function(x) pseudoIC(Ly = Ly_ob[[x]], hatY = hatY_ob[[x]],
                                           sigma_k = sigmas[[x]], class = class_new,
                                           M = M[[x]], K =K))
  #   calculate pAIC and pBIC
  hat_LL <- unlist(lapply(1:J, function(x) uICs[[x]]$hat_LL))
  pAIC <-  -2*sum(hat_LL) + 2*(K + sum(unlist(M)) )
  N <- length(unlist(Ly_ob))
  pBIC <-  -2*sum(hat_LL) + log(N)*(K + sum(unlist(M)))
  return(list(pAIC=pAIC, pBIC =pBIC,hat_LL = hat_LL, N_obs = N , M = M))


}
