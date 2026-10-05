#' @title Validate paired measurement and time lists
#'
#' @description
#' Checks that \code{Ly} and \code{Lt} are non-empty lists of the same length and
#' that each subject's two vectors are non-empty and equally long.
#'
#' @param Ly A list of \emph{n} numeric vectors of measurements.
#' @param Lt A list of \emph{n} numeric vectors of observation times.
#'
#' @return The number of subjects, invisibly used by \code{\link{RKCFC}}. Throws
#'   an error describing the first violated condition otherwise.
#'
#' @examples
#' univar_check_input(list(c(1, 2, 3), c(2, 4)), list(c(0, 0.5, 1), c(0.2, 0.8)))
#'
#' @export
univar_check_input <- function(Ly, Lt) {
  .rkcfc_validate_univariate(Ly, Lt)
}
