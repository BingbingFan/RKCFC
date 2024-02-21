#' @export
univar_check_input <- function(Ly, Lt){
  if(!is.list(Ly) ||  !is.list(Lt)){
    stop("Both Ly and Lt should be a list !!!\n")
  }
  if(length(Ly) != length(Lt)){
    stop("Ly should be the same lenth as Lt !!!\n")
  }
  n <- length(Ly)
  return(n)
}
