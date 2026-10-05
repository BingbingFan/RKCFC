#' @title Fit a K-centers functional clustering model
#'
#' @description
#' Clusters sparsely and irregularly observed univariate longitudinal
#' trajectories into \code{K} classes, alternating between class-specific
#' functional principal component analysis and reassignment by minimal squared
#' reconstruction error. With \code{robust = TRUE} the eigenspace is estimated by
#' Kendall's tau (RKCFC); with \code{robust = FALSE} by the PACE covariance
#' estimator (KCFC).
#'
#' @details
#' The algorithm targets the objective
#' \deqn{J(c) = \sum_i \sum_j \{ Y_{ij} - \hat{Y}_{ij,(c_i)} \}^2,}
#' the total squared reconstruction error at the current assignment. The
#' reassignment step minimises \eqn{J} over the labels with the class-specific
#' fits held fixed and is therefore monotone by construction. The re-estimation
#' step is not: the Kendall estimator does not minimise a sum of squared
#' reconstruction errors, and the local linear smoothing of the mean and
#' covariance surfaces carries a bandwidth-dependent bias that does not vanish in
#' \eqn{J}. Monotone descent of the whole iteration is therefore not guaranteed,
#' and neither is convergence of the partition. Use several starts and keep the
#' partition with the smallest \eqn{J}.
#'
#' Data are supplied either as the paired lists \code{Ly} and \code{Lt}, or as a
#' long-format \code{data} frame with the column names given by \code{id},
#' \code{ys} and \code{timevar}. Supply one form or the other, never both.
#'
#' A candidate partition placing fewer than ten curves in any class is rejected
#' and the previous partition retained, because the class-specific FPCA cannot be
#' estimated from fewer. This is reported through \code{conv}, which is
#' \code{FALSE} in that case; read it together with \code{iter}, since
#' \code{iter == maxIter} means the iteration cap bound the fit while
#' \code{iter < maxIter} with \code{conv == FALSE} means the small-class rule
#' fired. The two need opposite remedies and must not be pooled.
#'
#' @param Ly A list of \emph{n} numeric vectors, one per subject, holding the
#'   observed measurements. Supply with \code{Lt}.
#' @param Lt A list of \emph{n} numeric vectors, one per subject, holding the
#'   observation times in ascending order. Supply with \code{Ly}.
#' @param data A long-format data frame, one row per measurement. Supply with
#'   \code{id}, \code{ys} and \code{timevar} instead of \code{Ly} and \code{Lt}.
#' @param id Name of the subject identifier column in \code{data}.
#' @param ys Name of the measurement column in \code{data}.
#' @param timevar Name of the observation time column in \code{data}.
#' @param K Number of classes; a positive integer not exceeding the number of
#'   subjects. \code{K = 1} returns the trivial partition together with the
#'   pooled fit, which is useful as the reference point for
#'   \code{\link{silhouette_score}}.
#' @param robust Logical. \code{TRUE} estimates the eigenspace by Kendall's tau
#'   (RKCFC), \code{FALSE} by the PACE covariance estimator (KCFC).
#' @param optns A list of estimation options. When \code{robust = TRUE} it must
#'   supply \code{interval}, \code{nK}, \code{bw}, \code{bwK}, \code{nRegGrid}
#'   and \code{basis}; see \code{\link{KFPCA}} and the examples. Note that
#'   \code{bw}, the subject-level Nadaraya-Watson bandwidth, has no
#'   cross-validation path and must be given numerically. When
#'   \code{robust = FALSE} it is passed to \code{fdapace::FPCA}.
#' @param init.method Initialisation for the K-centers iteration. \code{"kmeans"}
#'   applies k-means to the pooled functional principal component scores;
#'   \code{"random"} assigns labels uniformly at random; \code{"two-staged"}
#'   applies k-means and falls back to random assignment when that returns a
#'   degenerate split. Under \code{robust = TRUE} the scores are Kendall scores;
#'   see \code{\link{KCFC_init}}.
#' @param init.labels Optional integer vector of length \emph{n} giving a
#'   starting partition directly. When supplied, \code{init.method} is ignored
#'   and \code{\link{KCFC_init}} is not called. Intended for studying
#'   initialisation: it is the only way to start from a partition the built-in
#'   schemes cannot produce, such as k-means on non-robust scores while fitting
#'   robustly, an external clustering, or the true partition as an oracle.
#'   Default \code{NULL}, which leaves the behaviour unchanged.
#' @param maxIter Maximum number of iterations, a positive integer.
#' @param parallel Logical; distribute the reassignment step over workers.
#' @param nmcores Number of workers when \code{parallel = TRUE}.
#' @param trace \code{0} silent, \code{1} reports convergence and refusals,
#'   \code{2} also reports each iteration.
#' @param conv.prop Convergence tolerance on the reclassification proportion. The
#'   iteration stops when the proportion of subjects changing label, relative to
#'   any previously visited partition, is at most this value; \code{0} therefore
#'   also detects cycles.
#' @param LOO Logical; when \code{TRUE}, the class-specific fit used to score
#'   subject \emph{i} against its own current class is re-estimated with subject
#'   \emph{i} removed, so that a subject's own data does not bias the decision to
#'   keep it where it already is. Methodologically preferable and substantially
#'   more expensive, since it refits the FPCA once per subject per iteration.
#'
#' @return A list with components
#'   \item{K}{the number of classes fitted}
#'   \item{res_clust}{a data frame of final class labels, carrying an \code{id}
#'     column when the model was fitted from \code{data}}
#'   \item{clustlist}{the partition at every iteration, in order}
#'   \item{init.method}{the initialisation used}
#'   \item{iter}{the iteration reached}
#'   \item{conv}{whether the convergence criterion was met}
#'   \item{stop_reason}{what ended the iteration: \code{"converged"} (the
#'     partition repeated the one before it), \code{"cycle"} (it repeated an
#'     earlier one, so the iteration is orbiting), \code{"small_cluster"} (a
#'     candidate put fewer than ten curves in a class and was refused),
#'     \code{"maxIter"} (the cap bound the fit), or \code{"trivial"} for
#'     \code{K = 1}. \code{conv} is \code{TRUE} for the first two}
#'   \item{cycle_length}{the length of the orbit when \code{stop_reason} is
#'     \code{"cycle"}, \code{NA} otherwise}
#'   \item{objective}{a data frame tracing the objective, one row per
#'     reassignment: \code{J_before} after the refit and before reassignment,
#'     \code{J_after} after it, \code{J_reestimation} the change attributable to
#'     the refit alone, and \code{reassignment_descends}, which must be
#'     \code{TRUE} in every row. See Details}
#'   \item{objective_final}{the objective at the final partition, computed
#'     without leave-one-out refitting whatever \code{LOO} was set to, so that
#'     fits are comparable across restarts}
#'   \item{optns}{the estimation options used}
#'   \item{FPCA_C}{the final class-specific FPCA fits, one per class}
#'   \item{data}{the \code{Ly} and \code{Lt} actually fitted}
#'
#' @seealso \code{\link{RKCFC_restarts}} for the multiple-start procedure,
#'   \code{\link{silhouette_score}} for choosing \code{K},
#'   \code{\link{KCFC_init}} for the initialisation, \code{\link{rre}} for
#'   scoring the reconstruction.
#'
#' @examples
#' \donttest{
#' data(pbc2)
#' # Serum bilirubin: sparse, irregular and strongly right-skewed.
#' pbc <- pbc2[!is.na(pbc2$serBilir), c("id", "year", "serBilir")]
#' keep <- names(which(table(pbc$id) >= 3))
#' pbc <- pbc[pbc$id %in% keep, ]
#'
#' interval <- range(pbc$year)
#' optns <- list(
#'   interval = interval, nK = 2, bw = 0.1 * diff(interval),
#'   # a fixed Kendall-surface bandwidth: GCV can choose one too small for the
#'   # class-specific fits on sparse data
#'   bwK = 0.2 * diff(interval),
#'   nRegGrid = 51,
#'   basis = fda::create.bspline.basis(interval, nbasis = 13, norder = 4,
#'             breaks = seq(interval[1], interval[2], length.out = 11))
#' )
#'
#' fit <- RKCFC(data = pbc, id = "id", ys = "serBilir", timevar = "year",
#'              K = 2, robust = TRUE, optns = optns,
#'              init.method = "two-staged", parallel = FALSE, LOO = FALSE)
#' table(fit$res_clust$g)
#' }
#'
#' @export
RKCFC <- function(Ly = NULL, Lt = NULL, data = NULL, id = NULL, ys = NULL,
                  timevar = NULL, K, robust = FALSE, optns,
                  init.method = "kmeans", init.labels = NULL, maxIter = 20,
                  parallel = TRUE, nmcores = NULL, trace = 1,
                  conv.prop = 0.01, LOO = TRUE) {
  list_input <- !is.null(Ly) || !is.null(Lt)
  data_input <- !is.null(data) || !is.null(id) || !is.null(ys) || !is.null(timevar)
  if (list_input == data_input) {
    stop("Supply either Ly and Lt, or data with id, ys, and timevar.", call. = FALSE)
  }
  if (length(K) != 1L || K < 1L || K != as.integer(K)) {
    stop("K must be a positive integer.", call. = FALSE)
  }
  if (length(maxIter) != 1L || maxIter < 1L || maxIter != as.integer(maxIter)) {
    stop("maxIter must be a positive integer.", call. = FALSE)
  }

  uid <- NULL
  if (list_input) {
    if (is.null(Ly) || is.null(Lt)) {
      stop("Both Ly and Lt must be supplied together.", call. = FALSE)
    }
    n <- univar_check_input(Ly, Lt)
  } else {
    if (is.null(data) || is.null(id) || is.null(ys) || is.null(timevar)) {
      stop("data, id, ys, and timevar must be supplied together.", call. = FALSE)
    }
    data <- data[, c(id, timevar, ys)]
    names(data) <- c("id", "time", "y")
    dt <- fdapace::MakeFPCAInputs(IDs = data$id, tVec = data$time, yVec = data$y)
    Ly <- dt$Ly
    Lt <- dt$Lt
    uid <- unlist(dt$Lid)
    n <- length(uid)
  }
  if (K > n) {
    stop("K cannot exceed the number of subjects.", call. = FALSE)
  }

  workers <- .rkcfc_validate_parallel_workers(parallel, nmcores, n)
  cluster <- NULL
  if (workers > 1L) {
    cluster <- parallel::makeCluster(workers)
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    .rkcfc_prepare_workers(
      cluster,
      c("univar_clust_update", "univar_predicted_yhat",
        ".rkcfc_solve_covariance", ".rkcfc_nonrobust_scores")
    )
  }

  if (K == 1L) {
    final_labels <- .rkcfc_one_cluster_labels(n)
    FPCA_C <- univar_FPCA_k(Ly, Lt, K, final_labels, robust, optns)
    res_clust <- if (is.null(uid)) {
      data.frame(g = final_labels)
    } else {
      data.frame(id = uid, g = final_labels)
    }
    return(list(
      K = K, res_clust = res_clust, clustlist = list(final_labels),
      init.method = init.method, iter = 0L, conv = TRUE,
      stop_reason = "trivial", cycle_length = NA_integer_,
      objective = NULL,
      objective_final = .rkcfc_objective(Ly, Lt, final_labels, FPCA_C, robust, optns),
      optns = optns, FPCA_C = FPCA_C, data = list(Ly = Ly, Lt = Lt)
    ))
  }

  # A supplied partition bypasses KCFC_init() entirely. This exists so that
  # initialisation can be STUDIED rather than only chosen from the three built-in
  # schemes, and three of the schemes the revision needs are not expressible any
  # other way:
  #
  #   - k-means on NON-ROBUST PACE scores while fitting robustly. That is exactly
  #     what the released 0.1.x package did regardless of the `robust` flag
  #     (defect B.1), and running it deliberately turns the correction from an
  #     assertion into a measurement. It cannot be reached through init.method
  #     because RKCFC() passes its own `robust` down to KCFC_init().
  #   - Ward, or any other external clustering.
  #   - the ORACLE start, at the true partition, which bounds what any
  #     initialiser could buy.
  #
  # Default NULL, so nothing about the existing behaviour changes.
  clust0 <- if (is.null(init.labels)) {
    KCFC_init(Ly, Lt, K, init.method, robust, optns)
  } else {
    if (length(init.labels) != length(Ly)) {
      stop("init.labels must have one entry per subject.", call. = FALSE)
    }
    lab <- as.integer(as.factor(init.labels))
    if (length(unique(lab)) != K) {
      stop("init.labels must contain exactly K distinct classes.", call. = FALSE)
    }
    lab
  }
  # Report what actually happened, not what was asked for. A results table that
  # records init.method = "kmeans" for a fit that started from a supplied oracle
  # partition is worse than one that records nothing.
  if (!is.null(init.labels)) init.method <- "supplied"
  .rkcfc_require_cluster_sizes(clust0, K)
  clustlist <- list(clust0)
  final_labels <- clust0
  iteration <- 1L
  conv <- FALSE
  stop_reason <- if (maxIter < 2L) "maxIter" else NA_character_
  cycle_length <- NA_integer_
  trace_rows <- list()

  if (maxIter >= 2L) {
    for (candidate_iteration in seq.int(2L, as.integer(maxIter))) {
      if (trace == 2) {
        cat(candidate_iteration, "th iteration\n")
      }
      ClustIds <- .rkcfc_cluster_ids(final_labels, K)
      ref_hatyc <- univar_FPCA_k(Ly, Lt, K, final_labels, robust, optns)
      if (is.null(cluster)) {
        clust_it <- lapply(seq_len(n), univar_clust_update,
                           ClustIds = ClustIds, Ly = Ly, Lt = Lt, K = K,
                           ref_hatyc = ref_hatyc, optns = optns,
                           robust = robust, LOO = LOO)
      } else {
        clust_it <- parallel::parLapply(
          cluster, seq_len(n), univar_clust_update,
          ClustIds = ClustIds, Ly = Ly, Lt = Lt, K = K,
          ref_hatyc = ref_hatyc, optns = optns, robust = robust, LOO = LOO
        )
      }
      candidate <- vapply(clust_it, function(x) x$clust_i, integer(1))

      # Recorded before the small-class break, so a refused iteration still
      # reports what the objective was doing when it was refused.
      trace_rows[[length(trace_rows) + 1L]] <- c(
        iteration = candidate_iteration,
        .rkcfc_objective_row(clust_it, final_labels, candidate)
      )

      candidate_ids <- .rkcfc_cluster_ids(candidate, K)
      if (any(lengths(candidate_ids) < 10L)) {
        if (trace >= 1) cat("Error: Less than 10 curves in a group.\n")
        stop_reason <- "small_cluster"
        break
      }
      clustlist[[candidate_iteration]] <- candidate
      final_labels <- candidate
      iteration <- candidate_iteration

      # The convergence test is unchanged: the candidate matching ANY previously
      # visited partition to within conv.prop. What is new is telling apart the
      # two ways that happens. Matching the partition immediately before is
      # convergence; matching only an earlier one is a cycle, and the two need
      # opposite remedies. classError() is invariant to label permutation, so
      # this compares partitions rather than labellings.
      history <- clustlist[-length(clustlist)]
      rates <- vapply(history, function(previous) {
        mclust::classError(previous, candidate)$errorRate
      }, numeric(1))
      hit <- which(rates <= conv.prop)
      if (length(hit) > 0L) {
        last_hit <- max(hit)
        if (last_hit == length(history)) {
          stop_reason <- "converged"
        } else {
          stop_reason <- "cycle"
          cycle_length <- as.integer(length(history) - last_hit + 1L)
        }
        if (trace >= 1) cat("Clustering algorithm converged.\n")
        conv <- TRUE
        break
      }
      if (candidate_iteration == as.integer(maxIter)) stop_reason <- "maxIter"
    }
  }

  FPCA_C <- univar_FPCA_k(Ly, Lt, K, final_labels, robust, optns)
  res_clust <- if (is.null(uid)) {
    data.frame(g = final_labels)
  } else {
    data.frame(id = uid, g = final_labels)
  }

  objective <- if (length(trace_rows)) {
    out <- as.data.frame(do.call(rbind, trace_rows))
    # J_after(b) -> J_before(b+1) is the re-estimation step, isolated. NA on the
    # first row because nothing precedes it.
    out$J_reestimation <- out$J_before - c(NA_real_, utils::head(out$J_after, -1L))
    out$reassignment_descends <- out$J_after <= out$J_before + 1e-8
    out
  } else {
    NULL
  }

  list(
    K = K, res_clust = res_clust, clustlist = clustlist,
    init.method = init.method, iter = iteration, conv = conv,
    stop_reason = stop_reason, cycle_length = cycle_length,
    objective = objective,
    objective_final = .rkcfc_objective(Ly, Lt, final_labels, FPCA_C, robust, optns),
    optns = optns, FPCA_C = FPCA_C, data = list(Ly = Ly, Lt = Lt)
  )
}
