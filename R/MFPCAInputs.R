MFPCAInputs <- function(id, y, timevar, data){
  names(data)[which(names(data) == id)] <- "id"
  uid <- unique(data$id)
  Ly <- list()
  Lt <- list()
  for (j in 1:length(y)) {
    Ly[[j]] <- list()
    Lt[[j]] <- list()
    for (i in 1:length(uid)) {
      dti <- data.frame(na.omit(data[data$id == uid[i],c(y[j],timevar)]))
      Ly[[j]][[i]] <- dti[, y[j] ]
      Lt[[j]][[i]] <- dti[, timevar ]
    }

  }
  Ldt <- list(id = as.list(uid), Ly = Ly, Lt = Lt)
  return(Ldt)
}
