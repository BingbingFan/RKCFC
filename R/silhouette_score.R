silhouette_score <- function(data, class, id, y, time,  robust = FALSE, option = NULL ){
  require(fdapace)
  if(robust){
    require(KFPCA)
  }
  require(dplyr)
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
  mu_func <- function(dt, option){
    if(robust){
      FPCA_obj <- KFPCA(Lt = dt$Lt,  Ly = dt$Ly, interval = option$interval, nK = option$nK, bw = option$bw,
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

  FPCA_res <- lapply(1:K, function(x) mu_func(dt = fdata_k[[x]], option = option) )
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

  return(list(a_i = a_i, b_i = b_i , silhouette = silhouette))
}
