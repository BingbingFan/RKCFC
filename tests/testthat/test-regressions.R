# Regression tests for the two defects present in RKCFC 0.1.x.
#
# Both were silent: they produced plausible numbers rather than errors, and both
# reached published results. Structural tests are used deliberately. A
# behavioural test would have to construct data on which the robust and
# non-robust starts provably disagree, which is exactly the kind of test that
# passes for the wrong reason; asserting that the wiring exists is narrower and
# cannot pass while the wiring is absent.

test_that("KCFC_init receives robust and does not hard-wire the non-robust engine", {
  # DEFECT (0.1.x): KCFC_init(Ly, Lt, K, init.method) called fdapace::FPCA()
  # unconditionally, so RKCFC(robust = TRUE) was started from a non-robust
  # estimator. K-centers converges to a local optimum, so the start selects the
  # basin and the basin selects the answer.
  expect_true(all(c("robust", "optns") %in% names(formals(KCFC_init))))

  body_text <- paste(deparse(body(KCFC_init)), collapse = "\n")
  expect_false(grepl("fdapace::FPCA", body_text, fixed = TRUE))
  expect_true(grepl("univar_FPCA_k", body_text, fixed = TRUE))

  # and RKCFC must actually forward both arguments
  call_text <- paste(deparse(body(RKCFC)), collapse = "\n")
  expect_match(call_text, "KCFC_init\\(Ly, Lt, K, init\\.method, robust, optns\\)")
})

test_that("silhouette_score passes bwK in the robust branch", {
  # DEFECT (0.1.x): bwK was omitted here and ONLY here, so the class mean
  # functions behind the coefficient were smoothed at the KFPCA default of
  # "GCV" rather than at the configured bandwidth -- different smoothing from
  # the model being scored. Affected robust = TRUE only.
  body_text <- paste(deparse(body(silhouette_score)), collapse = "\n")
  robust_call <- regmatches(
    body_text,
    regexpr("KFPCA\\(Lt = dt\\$Lt[^)]*\\)", body_text)
  )
  expect_length(robust_call, 1L)
  expect_match(robust_call, "bwK")
})

test_that("the KFPCA runtime is vendored rather than depended on", {
  # Asserted against the INSTALLED package rather than the source tree, because
  # under R CMD check the tests run from the check directory and the source tree
  # is not reachable from test_path().
  deps <- unlist(utils::packageDescription("RKCFC")[c("Depends", "Imports", "Suggests")])
  expect_false(any(grepl("\\bKFPCA\\b", deps)))

  expect_true(is.function(RKCFC::KFPCA))
  expect_true(is.function(utils::getS3method("predict", "KFPCA")))

  # The vendored runtime must be reachable internally, not merely re-exported.
  for (fn in c("MeanEst", "FPCscoreLSE", "GetMatList", "kernfun")) {
    expect_true(is.function(get(fn, envir = asNamespace("RKCFC"))), info = fn)
  }
})

test_that("no source file reaches for KFPCA as a package", {
  # Source-tree assertion: meaningful during development, skipped when the
  # package is checked from a tarball.
  source_dir <- file.path(testthat::test_path(), "..", "..", "R")
  skip_if_not(dir.exists(source_dir), "not running from the source tree")
  source_text <- unlist(lapply(
    list.files(source_dir, pattern = "\\.R$", full.names = TRUE),
    readLines, warn = FALSE
  ))
  expect_false(any(grepl("KFPCA::", source_text, fixed = TRUE)))
})

test_that("no multivariate entry point leaked into the univariate package", {
  exported <- getNamespaceExports("RKCFC")
  multivariate <- c("MKCFC", "MKCFC_init", "MFPCA_k", "MFPCAInputs",
                    "predicted_yhat", "clust_update", "CH_mv", "CH_mv_dt",
                    "sil_score_mv", "check_input")
  expect_length(intersect(exported, multivariate), 0L)
})

# --- the robust mean option ---------------------------------------------------
#
# These guard a REPRODUCIBILITY promise, not a defect. 537.9 core-h of simulation
# results were produced by this package before the option existed, and every one
# of those numbers must survive it. "The default is unchanged" is a claim about
# intent; the first test below is the measurement.

test_that("MeanEst(robust = FALSE) is the original least-squares call, untouched", {
  set.seed(404)
  Lt <- lapply(1:60, function(i) sort(stats::runif(sample(4:9, 1))))
  Ly <- lapply(Lt, function(t) 2 + 3 * t + stats::rnorm(length(t), 0, 0.3))
  g  <- seq(0, 1, length.out = 41)

  # The reference is fdapace::Lwls1D called EXACTLY as the pre-change body called
  # it -- no `win`, so its own integer default applies. If the default branch is
  # ever rerouted through the weighted call with unit weights, this fails.
  xin <- sort(unlist(Lt)); yin <- unlist(Ly)[order(unlist(Lt))]
  ref <- fdapace::Lwls1D(0.15, "epan", xin = xin, yin = yin, xout = g)

  expect_identical(RKCFC:::MeanEst(Lt, Ly, "epan", 0.15, g)$mean, ref)
  expect_identical(RKCFC:::MeanEst(Lt, Ly, "epan", 0.15, g, robust = FALSE)$mean, ref)
})

test_that("the mean-estimator choice travels on the object and reaches predict()", {
  # THE TRAP, and it is the reason this test exists. predict.KFPCA() NEVER READS
  # object$mean on the robust branch -- it re-estimates the mean from the fit's
  # own Lt/Ly. An option honoured in KFPCA() but not there would leave prediction
  # centred on the least-squares mean while every returned field claimed
  # otherwise, and every number would still have looked plausible.
  skip_on_cran()
  set.seed(505)
  n <- 120
  Lt <- lapply(seq_len(n), function(i) sort(stats::runif(sample(4:9, 1))))
  Ly <- lapply(Lt, function(t) 2 + 3 * t + stats::rnorm(length(t), 0, 0.3))
  for (i in sample(n, round(0.08 * n))) Ly[[i]] <- Ly[[i]] + 12   # contaminate

  fp <- fda::create.bspline.basis(c(0, 1), nbasis = 8, norder = 4)
  mk <- function(rb) KFPCA(Lt = Lt, Ly = Ly, interval = c(0, 1), nK = 2,
                           bw = 0.15, bwK = 0.15, bwmean = 0.15, nRegGrid = 51,
                           fdParobj = fp, meanrobust = rb)
  k0 <- mk(FALSE); k1 <- mk(TRUE)

  expect_false(k0$meanrobust)
  expect_true(k1$meanrobust)
  expect_false(isTRUE(all.equal(k0$mean, k1$mean)))

  p0 <- predict(k0, newLt = Lt[1:25], newLy = Ly[1:25], nK = 2, more = TRUE)
  p1 <- predict(k1, newLt = Lt[1:25], newLy = Ly[1:25], nK = 2, more = TRUE)
  expect_false(isTRUE(all.equal(p0$meanest_new, p1$meanest_new)))

  # An object fitted BEFORE this field existed has no $meanrobust, and must keep
  # behaving exactly as it did: least squares.
  kold <- k1; kold$meanrobust <- NULL
  pold <- predict(kold, newLt = Lt[1:25], newLy = Ly[1:25], nK = 2, more = TRUE)
  expect_equal(pold$meanest_new, p0$meanest_new)
})

test_that("the robust mean resists contamination the least-squares mean absorbs", {
  skip_on_cran()
  set.seed(606)
  n <- 150; truth <- function(t) 2 + 3 * t
  Lt <- lapply(seq_len(n), function(i) sort(stats::runif(sample(4:9, 1))))
  Ly <- lapply(Lt, function(t) truth(t) + stats::rnorm(length(t), 0, 0.3))
  for (i in sample(n, round(0.08 * n))) Ly[[i]] <- Ly[[i]] + 12

  g <- seq(0, 1, length.out = 51)
  e <- function(rb) sqrt(mean((RKCFC:::MeanEst(Lt, Ly, "epan", 0.15, g,
                                               robust = rb)$mean - truth(g))^2))
  # Generous factor: the point is that robustness works at all, not to pin a
  # ratio that a bandwidth change could move. Measured 2026-09-07: 0.043 against
  # 0.973, a factor of 23.
  expect_lt(e(TRUE), e(FALSE) / 3)
})

test_that("optns$meanrobust reaches KFPCA through univar_FPCA_k", {
  # The study passes settings as an `optns` list, so the option is only usable if
  # this call site forwards it. Structural: an optns list that predates the
  # option must still select least squares, via isTRUE(NULL).
  body_text <- paste(deparse(body(univar_FPCA_k)), collapse = "\n")
  # fixed = TRUE, not a regex: the string being matched is full of characters a
  # regex would treat as metacharacters, and escaping them is how this test
  # would quietly start matching nothing.
  expect_true(grepl("meanrobust = isTRUE(optns$meanrobust)", body_text,
                    fixed = TRUE))
})

# --- alternative robust estimators --------------------------------------------

test_that("scatter = 'sign' changes the covariance but not the location", {
  skip_on_cran()
  set.seed(707)
  n <- 100
  Lt <- lapply(seq_len(n), function(i) sort(stats::runif(sample(4:9, 1))))
  Ly <- lapply(Lt, function(t) 2 + 3 * t + stats::rt(1, 1.5) * sqrt(2) * sin(pi * t) +
                 stats::rnorm(length(t), 0, 0.3))
  fp <- fda::create.bspline.basis(c(0, 1), nbasis = 8, norder = 4)
  mk <- function(sc) KFPCA(Lt = Lt, Ly = Ly, interval = c(0, 1), nK = 2, bw = 0.15,
                           bwK = 0.15, bwmean = 0.15, nRegGrid = 51, fdParobj = fp,
                           scatter = sc)
  k1 <- mk("kendall"); k2 <- mk("sign")

  expect_identical(k1$scatter, "kendall")
  expect_identical(k2$scatter, "sign")
  expect_false(isTRUE(all.equal(k1$KendFun, k2$KendFun)))

  # THE POINT OF THE CONTRAST. Both estimators use the same mean, so a difference
  # between them is the scatter and nothing else. If a future edit gives the sign
  # branch its own location, this fails and it should.
  expect_equal(k1$mean, k2$mean)
})

test_that("the default scatter is kendall, and it is the pre-existing estimator", {
  expect_identical(eval(formals(KFPCA)$scatter)[1], "kendall")
  body_text <- paste(deparse(body(KFPCA)), collapse = "\n")
  # The kendall branch must still average signed PAIRWISE DIFFERENCES; the sign
  # branch signs a deviation from the mean. Both idioms must be present.
  expect_true(grepl("Xlist[[i]][1, ] - Xlist[[i]][j, ]", body_text, fixed = TRUE))
  expect_true(grepl("Xlist[[i]][1, ] - mu_i", body_text, fixed = TRUE))
})

test_that("normal_scores preserves shape and removes marginal tails", {
  set.seed(808)
  Ly <- lapply(1:40, function(i) stats::rt(sample(4:8, 1), df = 1.5))
  Lz <- normal_scores(Ly)

  expect_identical(lengths(Ly), lengths(Lz))
  expect_length(Lz, length(Ly))

  kurt <- function(x) mean((x - mean(x))^4) / stats::var(x)^2
  expect_gt(kurt(unlist(Ly)), 5)          # heavy by construction
  expect_lt(abs(kurt(unlist(Lz)) - 3), 1) # normal by construction

  # Strictly monotone in the pooled data: the transform may not reorder anything,
  # which is what makes a partition computed on it comparable with the raw scale.
  expect_false(is.unsorted(unlist(Lz)[order(unlist(Ly))]))

  # Nothing maps to an infinite score -- the (N + 1) denominator exists for this.
  expect_true(all(is.finite(unlist(Lz))))
})

test_that("normal_scores keeps NA in place and rejects bad input", {
  Ly <- list(c(1, 2, NA), c(4, 5))
  Lz <- normal_scores(Ly)
  expect_identical(lengths(Ly), lengths(Lz))
  expect_true(is.na(Lz[[1]][3]))
  expect_true(all(is.finite(unlist(Lz)[-3])))

  expect_error(normal_scores(1:10), "must be a list")
  expect_error(normal_scores(list(1)), "at least two finite")
})
