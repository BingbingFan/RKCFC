#' @export
check_input <- function(Ly, Lt){
  if(!is.list(Ly) ||  !is.list(Lt)){
    stop("Both Ly and Lt should be a list !!!\n")
  }
  if(length(Ly) != length(Lt)){
    stop("Ly should be the same lenth as Lt !!!\n")
  }
  J <- length(Ly)
  ny <- unlist(lapply(1:J, function(x) length(Ly[[x]]) ))
  nt <- unlist(lapply(1:J, function(x) length(Lt[[x]]) ))

  if(length(unique(ny)) != 1 ){
    stop("Each y in Ly should contains the same number of subjects!!!\n")
  }
  if(length(unique(nt)) != 1 ){
    stop("Each t in Lt should contains the same number of subjects!!!\n")
  }

  if(!all.equal(ny, nt) ){
    stop("Ly should contains the same number of subjects as Lt !!!\n")
  }else{
    n <- ny[[1]]
  }
  return(n)
}
