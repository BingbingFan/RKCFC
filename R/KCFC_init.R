#' @export
KCFC_init <- function(Ly, Lt, K, init.method){
  if(init.method == "kmeans"){
    fpcaObjY <- FPCA(Ly,Lt)
    mod <- kmeans(fpcaObjY$xiEst, centers = K)
    clust0 <- mod$cluster
  }
  if(init.method == "random"){
    clust0 <- sample(c(1:K),length(Ly),replace = T, prob = rep(1/K, K))
  }
  if(init.method == "two-staged"){
    fpcaObjY <- FPCA(Ly,Lt)
    mod <- kmeans(fpcaObjY$xiEst, centers = K)
    clust0 <- mod$cluster
    if (  min(table(clust0))  <  10 | length(unique(clust0) ) < K   ){
      clust0 <- sample(c(1:K),length(Ly),replace = T, prob = rep(1/K, K))
    }
  }
  return(clust0)
}
