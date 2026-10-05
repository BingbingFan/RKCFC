#' @title Modified silhouette coefficient for sparse functional clustering
#'
#' @description
#' The silhouette coefficient of Kaufman and Rousseeuw, adapted to sparse
#' irregularly observed curves. Equations (14) to (16) of the accompanying
#' manuscript. Used to choose the number of classes by fitting
#' \code{\link{RKCFC}} over a range of \code{K} and taking the maximiser.
#'
#' @details
#' Subjects observed at different times have no common coordinates, so the usual
#' subject-to-subject distances are not available. The within-class distance
#' \eqn{a_i} is therefore the mean squared distance of subject \emph{i}'s
#' observations to its own class mean function, and the between-class distance
#' \eqn{b_i} is the minimum of the same quantity taken over the other classes,
#' \deqn{b_i = \min_{l \neq k} n_i^{-1} \sum_j \{ Y_{ij} - \hat{\mu}^{(l)}(t_{ij}) \}^2.}
#'
#' Two limitations to state when reporting the result. The coefficient is
#' undefined at \code{K = 1}; a value of 0 there is a convention and should be
#' labelled as one rather than compared against. And because both distances are
#' measured to class *mean functions*, the criterion rewards classes separated in
#' level and penalises classes separated in shape or slope. A method producing
#' the latter can score lower than one producing the former on the same data
#' without being worse; interpret a small maximum against what the criterion
#' returns when the number of classes is known to be correct, not in isolation.
#' \code{\link{select_K}} supplies a score-space alternative for that reason.
#'
#' @note
#' Versions of this package before 0.2.0 omitted \code{bwK} in the robust branch,
#' and only there, so the class mean functions behind the coefficient were
#' smoothed at the \code{KFPCA} default rather than at the configured bandwidth.
#' Values computed for \code{robust = TRUE} with those versions are affected;
#' values for \code{robust = FALSE} are not.
#'
#' @param data A long-format data frame, one row per measurement.
#' @param class Name of the column holding the class label, coded \code{1:K}.
#' @param id Name of the subject identifier column.
#' @param y Name of the measurement column.
#' @param time Name of the observation time column.
#' @param robust Logical; \code{TRUE} estimates the class mean functions by
#'   \code{\link{KFPCA}}, \code{FALSE} by \code{fdapace::FPCA}. Use the same
#'   value that produced the partition being scored.
#' @param option Estimation options, as \code{optns} in \code{\link{RKCFC}}.
#'
#' @return A list with components
#'   \item{a_i}{the within-class distance of each subject}
#'   \item{b_i}{the between-class distance of each subject}
#'   \item{silhouette}{the mean silhouette coefficient over subjects}
#'
#' @references
#' Kaufman, L. and Rousseeuw, P. J. (2005). \emph{Finding Groups in Data: An
#' Introduction to Cluster Analysis}. Wiley, Hoboken.
#'
#' @seealso \code{\link{RKCFC}}, \code{\link{select_K}}
#'
#' @export
silhouette_score <- function(data, class, id, y, time, robust = FALSE, option = NULL) {
  data <- data[, c(id, class, time, y)]
  names(data) <- c("id", "g", "time", "y")
  K <- length(unique(data$g))
  n <- length(unique(data$id))
  # FPCA data
  # Indexed rather than subset()ed: subset() evaluates `g` non-standardly, which
  # R CMD check reads as an undefined global. Identical rows either way.
  dtk <- lapply(1:K, function(x) {
    data[which(data$g == x), , drop = FALSE]
  })

  fdata_k <- lapply(1:K, function(x)
    fdapace::MakeFPCAInputs(IDs = dtk[[x]]$id, yVec = dtk[[x]]$y, tVec = dtk[[x]]$time))

  # estimate mean function using FPCA or KFPCA
  mu_func <- function(dt, option) {
    if (robust) {
      # bwK must be passed here as it is at every other call site. Omitting it
      # smooths the covariance surface at the KFPCA default of "GCV" rather than
      # at the configured bandwidth, so the mean functions behind the
      # coefficient are estimated under different smoothing from the model being
      # scored. See the note in the help page.
      FPCA_obj <- KFPCA(Lt = dt$Lt, Ly = dt$Ly, interval = option$interval,
                        nK = option$nK, bw = option$bw, bwK = option$bwK,
                        nRegGrid = option$nRegGrid, fdParobj = option$basis)
      mu <- FPCA_obj$mean
      Grid <- FPCA_obj$RegGrid
    } else {
      FPCA_obj <- fdapace::FPCA(Ly = dt$Ly, Lt = dt$Lt, optns = option)
      mu <- FPCA_obj$mu
      Grid <- FPCA_obj$workGrid
    }
    list(mu = mu, Grid = Grid)
  }

  FPCA_res <- lapply(1:K, function(x) mu_func(dt = fdata_k[[x]], option = option))
  mus <- lapply(1:K, function(x) FPCA_res[[x]]$mu)
  Grids <- lapply(1:K, function(x) FPCA_res[[x]]$Grid)

  # for subject i, the mean squared L2 distance to each class mean function
  fdata <- fdapace::MakeFPCAInputs(IDs = data$id, yVec = data$y, tVec = data$time)

  L2_dis <- function(Lyi, Lti, mus, Grids) {
    mean_i_k <- lapply(1:K, function(x) stats::spline(Grids[[x]], mus[[x]], xout = Lti)$y)
    d_L2 <- lapply(1:K, function(x) (Lyi - mean_i_k[[x]])^2)
    unlist(lapply(1:K, function(x) mean(d_L2[[x]])))
  }

  L2_i <- lapply(1:n, function(x) L2_dis(Lyi = fdata$Ly[[x]], Lti = fdata$Lt[[x]], mus, Grids))

  uid <- data.frame(id = unlist(fdata$Lid))
  class_i <- data.frame(unique(data[, c("id", "g")]))
  class_i <- merge(uid, class_i, by = "id", all.x = TRUE)

  a_i <- unlist(lapply(1:n, function(x) L2_i[[x]][class_i$g[x]]))
  b_i <- unlist(lapply(1:n, function(x) min(L2_i[[x]][-class_i$g[x]])))

  silhouette <- mean(unlist(lapply(1:n, function(x)
    (b_i[x] - a_i[x]) / max(b_i[x], a_i[x]))))

  list(a_i = a_i, b_i = b_i, silhouette = silhouette)
}
