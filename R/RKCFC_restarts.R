#' @title Fit from several starts and keep the best
#'
#' @description
#' Runs \code{\link{RKCFC}} from \code{nstart} initialisations and returns the
#' fit attaining the smallest objective, together with a summary of what every
#' start produced. This is the recommended way to use the method.
#'
#' @details
#' The K-centers iteration is a Lloyd-type scheme converging to a local optimum,
#' so a single fit reports whichever basin its start landed in. A median near the
#' accuracy ceiling with a first quartile near chance is the signature of that,
#' and it is a property of the algorithm rather than of any particular dataset.
#'
#' Selecting by the objective makes the procedure deterministic given the seed:
#' rather than reporting one start and hoping, it reports the best partition
#' found under a stated rule, and reports the objective so the reader can see how
#' much better it was. The comparison uses \code{objective_final}, which is
#' computed without leave-one-out refitting whatever \code{LOO} is set to,
#' precisely so that starts are comparable.
#'
#' The \code{starts} data frame is what makes stability reportable. Read three
#' things from it: the spread of \code{J}, which says how much the start matters;
#' \code{prop_best}, the proportion of starts reaching the best partition, which
#' is the honest stability statistic; and \code{stop_reason}, split by cause,
#' since a fit stopped by the iteration cap and one stopped by a class falling
#' below ten curves need opposite remedies and must not be pooled.
#'
#' Failed starts are recorded, never fatal. A start whose initialisation is
#' degenerate throws inside \code{\link{RKCFC}}; here it becomes a row with
#' \code{error} set and \code{J} missing, so a failure rate can be reported
#' rather than silently changing the denominator.
#'
#' @param ... Arguments passed to \code{\link{RKCFC}}. Supply the data, \code{K},
#'   \code{robust} and \code{optns} exactly as you would there.
#' @param nstart Number of starts.
#' @param init.method Initialisation for each start. \code{"random"} is the usual
#'   choice here, since repeating a deterministic k-means start gives the same
#'   partition every time and measures nothing.
#' @param seed Optional integer. Set for a reproducible set of starts; the
#'   \emph{k}th start is seeded from \code{seed + k}, so an individual start can
#'   be re-run in isolation.
#' @param agreement Logical; also compute the pairwise adjusted Rand index among
#'   the successful partitions. \eqn{O(\mathrm{nstart}^2)} comparisons, cheap
#'   relative to fitting.
#'
#' @return A list with components
#'   \item{fit}{the \code{\link{RKCFC}} fit attaining the smallest objective}
#'   \item{starts}{a data frame with one row per start: \code{start}, \code{seed},
#'     \code{J}, \code{iter}, \code{conv}, \code{stop_reason},
#'     \code{aRI_to_best} and \code{error}}
#'   \item{prop_best}{the proportion of successful starts reaching the best
#'     partition, taking agreement above 0.999 as the same partition}
#'   \item{agreement}{the pairwise adjusted Rand index matrix among successful
#'     starts, or \code{NULL}}
#'   \item{n_failed}{the number of starts that threw}
#'
#' @seealso \code{\link{RKCFC}}
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 60
#' g <- rep(1:2, length.out = n)
#' grid <- seq(0, 1, length.out = 101)
#' Lt <- lapply(seq_len(n), function(i) sort(sample(grid, max(2, stats::rpois(1, 7)))))
#' Ly <- lapply(seq_len(n), function(i)
#'   2 * sin(4 * Lt[[i]]) + 5 * (g[i] == 2) + stats::rnorm(length(Lt[[i]]), 0, 0.5))
#'
#' optns <- list(dataType = "Sparse", methodMuCovEst = "smooth", nRegGrid = 21,
#'               methodSelectK = 2, useBinnedData = "OFF")
#' out <- RKCFC_restarts(Ly = Ly, Lt = Lt, K = 2, robust = FALSE, optns = optns,
#'                       nstart = 5, seed = 42, parallel = FALSE, trace = 0,
#'                       LOO = FALSE)
#' out$starts
#' out$prop_best
#' }
#'
#' @export
RKCFC_restarts <- function(..., nstart = 20, init.method = "random",
                           seed = NULL, agreement = TRUE) {
  if (length(nstart) != 1L || nstart < 1L || nstart != as.integer(nstart)) {
    stop("nstart must be a positive integer.", call. = FALSE)
  }
  nstart <- as.integer(nstart)

  seeds <- if (is.null(seed)) rep(NA_integer_, nstart) else as.integer(seed) + seq_len(nstart)

  fits <- vector("list", nstart)
  rows <- vector("list", nstart)
  for (b in seq_len(nstart)) {
    if (!is.na(seeds[b])) set.seed(seeds[b])
    fit <- tryCatch(RKCFC(..., init.method = init.method),
                    error = function(e) conditionMessage(e))
    if (is.character(fit)) {
      fits[[b]] <- NULL
      rows[[b]] <- data.frame(start = b, seed = seeds[b], J = NA_real_,
                              iter = NA_integer_, conv = NA,
                              stop_reason = NA_character_, error = fit,
                              stringsAsFactors = FALSE)
    } else {
      fits[[b]] <- fit
      rows[[b]] <- data.frame(start = b, seed = seeds[b],
                              J = fit$objective_final, iter = fit$iter,
                              conv = fit$conv, stop_reason = fit$stop_reason,
                              error = NA_character_, stringsAsFactors = FALSE)
    }
  }
  starts <- do.call(rbind, rows)

  ok <- which(!vapply(fits, is.null, logical(1)))
  if (length(ok) == 0L) {
    stop("Every start failed. First message: ", starts$error[1], call. = FALSE)
  }

  best <- ok[which.min(starts$J[ok])]
  best_labels <- fits[[best]]$res_clust$g

  starts$aRI_to_best <- NA_real_
  starts$aRI_to_best[ok] <- vapply(ok, function(b) {
    mclust::adjustedRandIndex(fits[[b]]$res_clust$g, best_labels)
  }, numeric(1))

  # "Reached the best partition" is agreement, not label equality: the labels
  # themselves are arbitrary up to permutation.
  prop_best <- mean(starts$aRI_to_best[ok] > 0.999)

  agree <- NULL
  if (isTRUE(agreement) && length(ok) > 1L) {
    agree <- matrix(1, length(ok), length(ok),
                    dimnames = list(starts$start[ok], starts$start[ok]))
    for (i in seq_along(ok)) for (j in seq_len(i - 1L)) {
      a <- mclust::adjustedRandIndex(fits[[ok[i]]]$res_clust$g,
                                     fits[[ok[j]]]$res_clust$g)
      agree[i, j] <- agree[j, i] <- a
    }
  }

  list(fit = fits[[best]], starts = starts, prop_best = prop_best,
       agreement = agree, n_failed = nstart - length(ok))
}
