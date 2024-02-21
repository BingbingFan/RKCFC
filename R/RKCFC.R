#' @export
RKCFC <- function(Ly = NULL, Lt = NULL, data = NULL, id = NULL, ys= NULL, timevar= NULL,
                  K ,robust = FALSE, optns, init.method =  "kmeans",  maxIter = 20,
                  parallel = TRUE, nmcores = NULL, trace = 1 ,
                  conv.prop= 0.01, LOO = TRUE){
  # require(KFPCA)
  # require(MASS)
  # library(parallel)
  # library(doParallel)
  # require(fdapace)
  # main ----
  if(parallel && is.null(nmcores)){
    stop(' "nmcores" should be specified when "parallel = TRUE"!')
  }

  if( !is.null(Ly) && !is.null(Lt)){
    n <- univar_check_input(Ly, Lt)
    J <- length(Ly)
  }
  if(!is.null(data) && !is.null(ys) && !is.null(timevar) && !is.null(id)){
    data <- data[,c(id, timevar, ys)]
    names(data) <- c("id","time","y")
    dt <- fdapace::MakeFPCAInputs(IDs = data$id, tVec = data$time, yVec = data$y)
    Ly <- dt$Ly
    Lt <- dt$Lt
    uid <- unlist(dt$Lid)
    n <- length(uid)
  }
  # initial cluster # Kmeans
  if(K == 1){
    if(!is.null(uid)){
      res_clust <- data.frame(id = uid, g = 1)
    }else{
      res_clust <- data.frame(g = 1)
    }

    FPCA_C <- univar_FPCA_k( Ly, Lt, K, res_clust$g, robust, optns )
    it  <-  0
    conv <- T
  }else{
    #
    clust0 <- KCFC_init(Ly, Lt, K, init.method)

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
        # clust
        if(parallel){
          cl <- parallel::makeCluster(nmcores)
          doParallel::registerDoParallel(cl)
          clust_it <- foreach::foreach(x=1:n) %dopar%
            univar_clust_update(m = x, ClustIds, Ly, Lt, K, ref_hatyc, optns,robust)
          parallel::stopCluster(cl)
        }else{
          clust_it <- lapply(1:n, function(x) univar_clust_update(m = x, ClustIds, Ly, Lt, K, ref_hatyc, optns,robust,LOO))
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

        if ( (it > 1) && any( sapply(clustlist[1:(it - 1)], function(u) {
          mclust::classError(u,clustlist[[it]])$errorRate <= conv.prop   } ) )   ) {
          conv <- TRUE
          if(trace >= 1){
            cat("Clustering algorithm converged with", it," times iteration !\n")
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

    #
    if(min(table(clustlist[[it]]))  >=  10 & length(unique(clustlist[[it-1]]) ) == K){
      FPCA_C <- univar_FPCA_k( Ly, Lt,  K, clustlist[[it-1]], robust, optns )
    }else{
      FPCA_C <- NULL
    }
    if(!is.null(uid)){
      res_clust <- data.frame(id = uid, g = clustlist[[it-1]])
    }else{
      res_clust <- data.frame(g = clustlist[[it-1]])
    }

  }



  return(list(K = K, res_clust = res_clust,clustlist = clustlist, init.method = init.method, iter = it ,
              conv = conv, optns = optns, FPCA_C = FPCA_C, data = list(Ly=Ly, Lt = Lt )))

}
