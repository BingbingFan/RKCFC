#' @title Robust K-Centers Functional Clustering for Sparse Longitudinal Data
#'
#' @description
#' Clustering of sparsely and irregularly observed univariate longitudinal
#' trajectories that may depart from normality. The subject's latent trajectory
#' is modelled by a Karhunen-Loeve expansion with class-specific mean function
#' and eigenfunctions, and subjects are assigned to the class that reconstructs
#' their observations with the smallest squared error.
#'
#' \code{\link{RKCFC}} fits both methods through one entry point. With
#' \code{robust = FALSE} it is the K-centers functional clustering (KCFC) model
#' of Chiou and Li (2007), whose eigenspace is estimated by the PACE covariance
#' estimator. With \code{robust = TRUE} it is the proposed RKCFC model, whose
#' eigenspace is estimated by the Kendall's tau operator of Zhong et al. (2022).
#'
#' @section What is and is not robust:
#' Only the covariance surface and the eigenfunctions are estimated robustly.
#' The mean function (\code{MeanEst}, local linear least squares), the principal
#' component scores (\code{FPCscoreLSE}, least squares), the subject-level
#' Nadaraya-Watson pre-smoothing, and the squared-error class-assignment rule are
#' not. RKCFC is therefore robust to non-Gaussian principal component scores and
#' to a covariance operator without finite fourth moments; it is not
#' distribution-free in general, and it is not protected against gross outliers
#' in the observed values. State this rather than inferring more from the name.
#'
#' @section Choosing the number of classes:
#' \code{\link{silhouette_score}} implements the modified silhouette coefficient
#' for sparse designs, in which the within-class distance is the mean squared
#' distance of a subject's observations to its own class mean function and the
#' between-class distance is the minimum of the same quantity over the other
#' classes. \code{\link{select_K}} supplies score-space alternatives, since a
#' single weak criterion carrying the whole decision is fragile.
#'
#' @section Reproducibility machinery:
#' \code{\link{rre}}, \code{\link{robust_normalize}} and \code{\link{gen_scores}}
#' are the reporting metric, the preprocessing transform and the simulation score
#' generator used by the accompanying manuscript. They ship here so that the
#' version producing the published numbers is the version that is installed.
#'
#' @references
#' Chiou, J.-M. and Li, P.-L. (2007). Functional clustering and identifying
#' substructures of longitudinal data. \emph{Journal of the Royal Statistical
#' Society Series B}, \bold{69}(4), 679-699.
#'
#' Zhong, R., Liu, S., Li, H. and Zhang, J. (2022). Robust functional principal
#' component analysis for non-Gaussian longitudinal data. \emph{Journal of
#' Multivariate Analysis}, \bold{189}, 104864.
#'
#' Yao, F., Mueller, H.-G. and Wang, J.-L. (2005). Functional data analysis for
#' sparse longitudinal data. \emph{Journal of the American Statistical
#' Association}, \bold{100}(470), 577-590.
#'
#' @keywords internal
"_PACKAGE"
