univar_FPCA_k <- function( Ly, Lt,K,  clust_f , robust, optns){
  ClustIds <- lapply(unique(clust_f), function(x) which(clust_f == x)  )
  ref_hatyc <- list()
  for(u in 1:K){
    uu <- ClustIds[[u]]
    yy <-  lapply( uu,  function(u) Ly[[u]])
    tt <- lapply( uu,  function(u) Lt[[u]])

    if(robust){

      ref_hatyc[[u]] <- KFPCA(Lt = tt,  Ly = yy, interval = optns$interval,
                              nK = optns$nK, bw = optns$bw, bwK = optns$bwK,
                              nRegGrid = optns$nRegGrid, fdParobj = optns$basis)
    }else{
      ref_hatyc[[u]] <- FPCA(yy, tt, optns =  optns )
    }

  }

  return(ref_hatyc)
}
