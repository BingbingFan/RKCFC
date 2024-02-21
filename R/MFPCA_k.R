#' @export
MFPCA_k <- function(u, Ly, Lt,  clust_f, optns,robust ){
  xiEst <- multi_score <- NULL

  ClustIds <- lapply(unique(clust_f), function(x) which(clust_f == x)  )
  J <- length(Ly)
  ref_hatyc <- list()
  uu <- ClustIds[[u]]
  yy <- list()
  tt <- list()
  for(j in 1:J){
    yy[[j]] <-  lapply( uu,  function(x) Ly[[j]][[x]])
    tt[[j]] <- lapply( uu,  function(x) Lt[[j]][[x]])
  }

  if(robust){
    ref_hatyc <- lapply(1:J, function(x)
      KFPCA(Lt = tt[[x]],  Ly = yy[[x]], interval = optns$interval,
            nK = optns$nK, bw = optns$bw, bwK = optns$bwK,
            nRegGrid = optns$nRegGrid, fdParobj = optns$basis) )

    xiEst <- lapply(1:J, function(x) ref_hatyc[[x]]$score )
    n_M <- lapply(1:J, function(x) dim(ref_hatyc[[x]]$score)[2] )
    phi <- lapply(1:J, function(x) ref_hatyc[[x]]$FPC_dis )
    mux <- lapply(1:J, function(x) ref_hatyc[[x]]$mean )

  }else{

    ref_hatyc <- lapply(1:J, function(x)
      fdapace::FPCA(yy[[x]], tt[[x]], optns =  optns ))

    xiEst <- lapply(1:J, function(x) ref_hatyc[[x]]$xiEst )
    n_M <- lapply(1:J, function(x) ref_hatyc[[x]]$selectK )
    phi <- lapply(1:J, function(x) ref_hatyc[[x]]$phi )
    mux <- lapply(1:J, function(x) ref_hatyc[[x]]$mu )

  }
  full_xiEst <- vector()
  for(j in 1:J ){
    full_xiEst <- cbind(full_xiEst, xiEst[[j]])
  }

  Z <-  t(full_xiEst) %*% full_xiEst / (nrow(full_xiEst)-1)
  #
  eigen.Z  <-  eigen(Z)
  eigen.values  <-  eigen.Z$values
  pve  <-  cumsum(eigen.values)/sum(eigen.values)
  Cms <-  eigen.Z$vectors

  mfpca.score = function(predXi, Cms){
    rho = matrix(NA, nrow = nrow(predXi), ncol=dim(Cms)[2])
    for(i in 1:nrow(predXi)){
      for(m in 1:dim(Cms)[2]){
        rho[i,m] = predXi[i,] %*% Cms[,m]
      }
    }
    return(rho)
  }
  multi_score <- mfpca.score(full_xiEst, Cms)

  return(list(multi_score = multi_score,
              ref_hatyc = ref_hatyc  ))
}
