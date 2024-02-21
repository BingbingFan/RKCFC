univar_clust_update <- function(m,ClustIds,Ly,Lt,K,ref_hatyc,optns,robust,LOO){
  # m: index of curve (m-th curve )
  # ClustIds
  # packages and functions for parallel ----
  require(MASS)
  require(fdapace)
  require(KFPCA)
  # start of cluster updating---
  ith_yy <- Ly[[m]]
  ith_tt <- Lt[[m]]
  hatyc <-  list()
  for(u in 1:K){	## index for cluster group
    uu <- ClustIds[[u]]  ## index of curves in cluster u
    if(LOO){
      if( sum(m == uu)== 1  ){ ## If m-th curve is in the group u
        yy <- lapply( uu,  function(x) Ly[[x]])
        yy[[ which(m==uu)  ]] <- NULL
        tt <- lapply( uu,  function(x) Lt[[x]])
        tt[[ which(m==uu)  ]] <- NULL
        # FPCA without subject m
        if(robust){
          temp_fpca <- KFPCA(Lt = tt,  Ly = yy, interval = optns$interval,
                             nK = optns$nK, bw = optns$bw, bwK = optns$bwK,
                             nRegGrid = optns$nRegGrid, fdParobj = optns$basis)
        }else{
          temp_fpca <- FPCA(yy, tt, optns =  optns )

        }
        yhatc_pred <- univar_predicted_yhat(temp_fpca, ith_yy,  ith_tt, robust, optns)
      }

      if( sum(m == uu)== 0  ){## If m-th curve is not in the group u

        yhatc_pred <-  univar_predicted_yhat(ref_hatyc[[u]], ith_yy,  ith_tt, robust, optns)
      }
    }else{
      yhatc_pred <-  univar_predicted_yhat(ref_hatyc[[u]], ith_yy,  ith_tt, robust, optns)
    }
    hatyc[[u]] <- yhatc_pred$fittedX
  }
  # calculate predict error
  ss <- list()
  for(u in 1:K){
    uu <- ClustIds[[u]]
    if(sum(m == uu)==1){
      err <- hatyc[[u]]-  ( ith_yy )
    }
    ss[[u]]  <- sum((hatyc[[u]] - ( ith_yy ))^2 )
  }
  clust_i <- which.min(ss)
  return(list(clust_i = clust_i, err = err,
              hatyc = hatyc, ss = ss))
}
