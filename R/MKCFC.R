MKCFC <- function(Ly = NULL, Lt = NULL, data = NULL, id = NULL, ys= NULL, timevar= NULL,
                  K ,robust, optns, init.method =  "kmeans",  maxIter = 20,
                  parallel = TRUE, nmcores = NULL, trace = 1 , conv.prop = 0,
                  LOO = TRUE ){
  if(parallel && is.null(nmcores)){
    stop(' "nmcores" should be specified when "parallel = TRUE"!')
  }

  suppressMessages(require(MASS))
  suppressMessages(require(parallel))
  suppressMessages(require(doParallel))
  suppressMessages(require(fdapace))
  suppressMessages(require(KFPCA))
  if( !is.null(Ly) && !is.null(Lt) ){
    n <- check_input(Ly, Lt)
    J <- length(Ly)
    uid <- NULL
  }
  if(!is.null(data) && !is.null(ys) && !is.null(timevar) && !is.null(id) ){
    dt <- MFPCAInputs(id, ys,timevar ,data )
    Ly <-  dt$Ly
    Lt <-  dt$Lt
    uid <- unlist(dt$id)
    n <- length(uid)
    J <- length(Ly)
  }
  # initial cluster # Kmeans
  if(K == 1){
    if(!is.null(uid)){
      res_clust <- data.frame(id = uid, g = 1)
    }else{
      res_clust <- data.frame(g = rep(1,n))
    }
    it  <-  0
    MFPCA_C <- lapply(1:K, function(x) MFPCA_k(u = x, Ly, Lt,  res_clust$g, optns,robust ))
    conv <- T

  }else{
    clust0 <- MKCFC_init(Ly, Lt, K, init.method)

    clustlist <- list()
    it <- 1
    clustlist[[1]] <- clust0
    if (  min(table(clust0))  <  10 | length(unique(clustlist[[it]]) ) < K   ) {
      conv <- FALSE
      print('Stop at it = 1 & unbalabced result ! \n')
    }
    # iteration
    if (   min(table(clust0))  >=  10  && length(unique(clustlist[[it]]) ) == K  ) {
      for (it in 2:(maxIter)) {
        if(trace == 2){
          cat(paste(it,'th iteration \n'))
        }
        ClustIds <- lapply(unique(clustlist[[it-1]]), function(x) which(clustlist[[it-1]] == x)  )
        # ref_hatyc
        ref_hatyc <- list()
        for(u in 1:K){
          uu <- ClustIds[[u]]
          yy <- list()
          tt <- list()
          for(j in 1:J){
            yy[[j]] <-  lapply( uu,  function(u) Ly[[j]][[u]])
            tt[[j]] <- lapply( uu,  function(u) Lt[[j]][[u]])
          }

          if(robust){

            ref_hatyc[[u]] <- lapply(1:J, function(x)
              KFPCA(Lt = tt[[x]],  Ly = yy[[x]], interval = optns$interval,
                    nK = optns$nK, bw = optns$bw, bwK = optns$bwK,
                    nRegGrid = optns$nRegGrid, fdParobj = optns$basis) )

          }else{
            ref_hatyc[[u]] <- lapply(1:J, function(x)
              fdapace::FPCA(yy[[x]], tt[[x]], optns =  optns ))
          }

        }
        # clust
        if(parallel){
          cl <- makeCluster(nmcores)
          registerDoParallel(cl)
          clust_it <- foreach(x=1:n) %dopar%
            clust_update(m = x, ClustIds, Ly, Lt, K, ref_hatyc, optns,robust,LOO)
          stopCluster(cl)
        }else{
          clust_it <- lapply(1:n, function(x) clust_update(m = x, ClustIds, Ly, Lt, K, ref_hatyc, optns,robust,LOO))
        }

        clustlist[[it]] <- unlist(lapply(1:n, function(x) clust_it[[x]]$clust_i))

        ##break information
        conv <- FALSE
        if( length(unique(clustlist[[it]]) ) < K  ){
          if(trace >= 1){
            cat("Error: Number of group is less than K !\n")
          }
          break
        }

        # if ( (it > 1) && any(sapply(clustlist[1:(it - 1)], function(u) all(u == clustlist[[it]])))  ) {
        if ( (it > 1) && any( sapply(clustlist[1:(it - 1)], function(u) {
          mclust::classError(u,clustlist[[it]])$errorRate <= conv.prop   } ) )   ){
          if(trace >= 1){
            cat("Clustering Algorithm converged !\n")
            conv <- TRUE
          }
          break
        }
        if (   min(table(clustlist[[it]]))  <  10  ) {
          if(trace >= 1){
            cat('Error: Less than 10 curve in a group !\n')
          }
          break
        }
      }
    }
    if(min(table(clustlist[[it]]))  >=  10 & length(unique(clustlist[[it-1]]) ) == K){
      MFPCA_C <- lapply(1:K, function(x) MFPCA_k(u = x, Ly, Lt,  clustlist[[it]], optns,robust ))
      # multi_score <- lapply(1:K, function(x) MFPCA_C[[x]]$multi_score)
      # ref_hatyc <- lapply(1:K, function(x) MFPCA_C[[x]]$ref_hatyc)
    }else{
      MFPCA_C <- multi_score <- ref_hatyc <- NULL
    }

    if(!is.null(uid)){
      res_clust <- data.frame(id = uid, g = clustlist[[it]])
    }else{
      res_clust <- data.frame(g = clustlist[[it]])
    }

  }
  return(list(K = K, res_clust = res_clust, init.method = init.method, iter = it ,
              conv = conv, optns = optns, MFPCA_C = MFPCA_C )
  )

}
