#' @export
MKCFC_init <- function(Ly, Lt, K, init.method){
  if(init.method == "kmeans"){
    fpcaObjY <- FPCA(Ly[[1]], Lt[[1]])
    mod <- kmeans(fpcaObjY$xiEst, centers = K)
    c0 <- mod$cluster
  }
  if(init.method == "random"){
    c0 <- sample(c(1:K),length(Ly[[1]]),replace = T, prob = rep(1/K, K))
  }
  if(init.method == "two-staged"){
    fpcaObjY <- fdapace::FPCA(Ly[[1]], Lt[[1]])
    mod <- kmeans(fpcaObjY$xiEst, centers = K)
    c0 <- mod$cluster
    if (  min(table(c0))  <  10 | length(unique(c0) ) < K   ){
      c0 <- sample(c(1:K),length(Ly[[1]]),replace = T, prob = rep(1/K, K))
    }
  }
  return(c0)
}
