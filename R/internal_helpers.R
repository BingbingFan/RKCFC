# Internal helpers. Not exported; documented here rather than in man/ because
# they are implementation details of the fitting path.

.rkcfc_validate_univariate <- function(Ly, Lt) {
  if (!is.list(Ly) || !is.list(Lt) || length(Ly) == 0L || length(Lt) == 0L) {
    stop("Ly and Lt must be non-empty lists.", call. = FALSE)
  }
  if (length(Ly) != length(Lt)) {
    stop("Ly and Lt must contain the same number of subjects.", call. = FALSE)
  }
  valid <- vapply(seq_along(Ly), function(i) {
    length(Ly[[i]]) > 0L && length(Ly[[i]]) == length(Lt[[i]])
  }, logical(1))
  if (!all(valid)) {
    stop("Each subject's Ly and Lt vectors must be non-empty and have the same length.", call. = FALSE)
  }
  length(Ly)
}

.rkcfc_cluster_ids <- function(labels, K) {
  if (length(labels) == 0L || length(K) != 1L || K < 1L || K != as.integer(K)) {
    stop("K must be a positive integer and labels must be non-empty.", call. = FALSE)
  }
  lapply(seq_len(as.integer(K)), function(k) which(labels == k))
}

.rkcfc_require_cluster_sizes <- function(labels, K, minimum = 10L) {
  ids <- .rkcfc_cluster_ids(labels, K)
  if (any(lengths(ids) < minimum)) {
    stop(sprintf("Each cluster must contain at least %d curves.", minimum), call. = FALSE)
  }
  invisible(ids)
}

.rkcfc_one_cluster_labels <- function(n) {
  if (length(n) != 1L || n < 1L || n != as.integer(n)) {
    stop("n must be a positive integer.", call. = FALSE)
  }
  rep.int(1L, as.integer(n))
}

.rkcfc_validate_parallel_workers <- function(parallel, nmcores, n) {
  if (!is.logical(parallel) || length(parallel) != 1L || is.na(parallel)) {
    stop("parallel must be TRUE or FALSE.", call. = FALSE)
  }
  if (!parallel) {
    return(1L)
  }
  if (is.null(nmcores) || length(nmcores) != 1L || !is.finite(nmcores) ||
      nmcores < 1L || nmcores != as.integer(nmcores)) {
    stop("nmcores must be a positive integer when parallel = TRUE.", call. = FALSE)
  }
  min(as.integer(nmcores), as.integer(n))
}

# Make the workers able to run the update functions. Two modes must both work,
# and they fail differently, which is why this is not a plain clusterExport().
#
#   installed  The update functions live in the RKCFC namespace. When a worker
#              deserialises one, it resolves that namespace from the WORKER's
#              own .libPaths(). If a different build of RKCFC is installed
#              there, the worker silently runs THAT one. The failure surfaces
#              as "unused argument", not "object not found", because the worker
#              found a real function of the same name from an older version.
#              Pushing the master's library path fixes it, and loading the
#              namespace explicitly makes the dependency visible rather than
#              relying on deserialisation to trigger it.
#
#   sourced    R/*.R has been source()d, so the functions are in the global
#              environment and there may be no installed RKCFC at all, or only
#              a stale one. Here the closures themselves must be shipped, and
#              loading a namespace would be actively wrong.
#
# A hand-maintained export list is only correct for the sourced case, and
# forgetting a name there fails solely under parallel = TRUE.
.rkcfc_prepare_workers <- function(cluster, exports) {
  env <- topenv(environment(.rkcfc_prepare_workers))
  if (isNamespace(env)) {
    parallel::clusterCall(cluster, function(paths) .libPaths(paths), .libPaths())
    pkg <- environmentName(env)
    parallel::clusterCall(cluster, function(p) {
      suppressMessages(requireNamespace(p, quietly = TRUE))
      invisible(NULL)
    }, pkg)
    return(invisible("namespace"))
  }
  parallel::clusterExport(cluster, exports, envir = env)
  invisible("global")
}

# Solve S x = residual for a covariance matrix S that may be near-singular.
# Direct solve first, then a ridge scaled to the matrix, then a generalised
# inverse. Always forming a pseudo-inverse would be both slower and less
# accurate in the well-conditioned case, which is the common one.
.rkcfc_solve_covariance <- function(covariance, residual,
                                    tol = sqrt(.Machine$double.eps)) {
  if (!is.matrix(covariance) || nrow(covariance) != ncol(covariance)) {
    stop("covariance must be a square matrix.", call. = FALSE)
  }
  if (length(residual) != nrow(covariance)) {
    stop("residual must have one value per covariance row.", call. = FALSE)
  }
  covariance <- (covariance + t(covariance)) / 2
  direct <- tryCatch(solve(covariance, residual), error = function(e) NULL)
  if (!is.null(direct)) {
    return(direct)
  }
  scale <- max(1, max(abs(diag(covariance))))
  ridge <- sqrt(.Machine$double.eps) * scale
  regularized <- tryCatch(
    solve(covariance + diag(ridge, nrow(covariance)), residual),
    error = function(e) NULL
  )
  if (!is.null(regularized)) {
    return(regularized)
  }
  MASS::ginv(covariance, tol = tol) %*% residual
}

# Conditional expectation of the FPC scores given the observed values, Yao,
# Mueller and Wang (2005) equation (4). Least squares, and therefore NOT robust:
# see the package-level help.
.rkcfc_nonrobust_scores <- function(phi, lambda, sigma2, residual) {
  lambda <- pmax(as.numeric(lambda), 0)
  if (!is.matrix(phi) || ncol(phi) != length(lambda)) {
    stop("phi columns must match the number of eigenvalues.", call. = FALSE)
  }
  if (length(residual) != nrow(phi)) {
    stop("residual length must match the number of observation times.", call. = FALSE)
  }
  noise <- max(as.numeric(sigma2), .Machine$double.eps)
  weighted_phi <- sweep(phi, 2L, sqrt(lambda), "*")
  covariance <- tcrossprod(weighted_phi) + diag(noise, nrow(phi))
  as.vector(lambda * crossprod(phi, .rkcfc_solve_covariance(covariance, residual)))
}
