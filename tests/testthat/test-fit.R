# End-to-end fits on a small, well-separated design. These are wiring tests, not
# accuracy claims: n is far too small and the separation far too easy to say
# anything about the method's performance.

sim_two_class <- function(n = 60, d = 6, arm = 1, seed = 1) {
  set.seed(seed)
  mu <- list(function(t) 2 * sin(4 * t),
             function(t) 2 * sin(4 * t) + d)
  phi <- list(
    list(function(t) sqrt(2) * sin(2 * pi * t), function(t) sqrt(2) * sin(3 * pi * t)),
    list(function(t) sqrt(2) * cos(2 * pi * t), function(t) sqrt(2) * cos(3 * pi * t))
  )
  g <- rep(1:2, length.out = n)
  xi <- list(gen_scores(n, arm, diag(c(9, 4))), gen_scores(n, arm, diag(c(8.1, 3.6))))
  grid <- seq(0, 1, length.out = 101)
  Lt <- Ly <- vector("list", n)
  for (i in seq_len(n)) {
    ni <- max(2L, stats::rpois(1, 7))
    Lt[[i]] <- sort(sample(grid, ni))
    k <- g[i]
    Ly[[i]] <- mu[[k]](Lt[[i]]) +
      xi[[k]][i, 1] * phi[[k]][[1]](Lt[[i]]) +
      xi[[k]][i, 2] * phi[[k]][[2]](Lt[[i]]) +
      stats::rnorm(ni, 0, 0.1)
  }
  list(Ly = Ly, Lt = Lt, g = g)
}

robust_optns <- list(
  interval = c(0, 1), nK = 2, bw = 0.1, bwK = 0.15, bwmean = 0.1, nRegGrid = 21,
  basis = fda::create.bspline.basis(c(0, 1), nbasis = 7, norder = 4,
                                    breaks = seq(0, 1, length.out = 5))
)
plain_optns <- list(dataType = "Sparse", methodMuCovEst = "smooth",
                    nRegGrid = 21, methodSelectK = 2, useBinnedData = "OFF")

test_that("RKCFC recovers a well-separated two-class partition, robust branch", {
  skip_on_cran()
  d <- sim_two_class()
  fit <- RKCFC(Ly = d$Ly, Lt = d$Lt, K = 2, robust = TRUE, optns = robust_optns,
               init.method = "two-staged", maxIter = 10, parallel = FALSE,
               trace = 0, LOO = FALSE)

  expect_s3_class(fit$res_clust, "data.frame")
  expect_length(fit$res_clust$g, length(d$Ly))
  expect_length(fit$FPCA_C, 2L)
  expect_true(fit$iter >= 1L)
  expect_gt(mclust::adjustedRandIndex(fit$res_clust$g, d$g), 0.5)
})

test_that("RKCFC recovers a well-separated two-class partition, non-robust branch", {
  skip_on_cran()
  d <- sim_two_class()
  fit <- RKCFC(Ly = d$Ly, Lt = d$Lt, K = 2, robust = FALSE, optns = plain_optns,
               init.method = "two-staged", maxIter = 10, parallel = FALSE,
               trace = 0, LOO = FALSE)
  expect_gt(mclust::adjustedRandIndex(fit$res_clust$g, d$g), 0.5)
})

test_that("K = 1 short-circuits to the pooled fit", {
  skip_on_cran()
  d <- sim_two_class(n = 30)
  fit <- RKCFC(Ly = d$Ly, Lt = d$Lt, K = 1, robust = TRUE, optns = robust_optns,
               parallel = FALSE, trace = 0, LOO = FALSE)
  expect_equal(fit$iter, 0L)
  expect_true(fit$conv)
  expect_length(fit$FPCA_C, 1L)
  expect_equal(unique(fit$res_clust$g), 1L)
})

test_that("the data-frame interface matches the list interface and keeps ids", {
  skip_on_cran()
  d <- sim_two_class(n = 40)
  long <- data.frame(
    subject = rep(seq_along(d$Ly), lengths(d$Ly)),
    t = unlist(d$Lt), value = unlist(d$Ly)
  )
  fit <- RKCFC(data = long, id = "subject", ys = "value", timevar = "t",
               K = 1, robust = FALSE, optns = plain_optns,
               parallel = FALSE, trace = 0, LOO = FALSE)
  expect_true("id" %in% names(fit$res_clust))
  expect_equal(nrow(fit$res_clust), 40L)
})

test_that("RKCFC validates its arguments", {
  d <- sim_two_class(n = 20)
  expect_error(RKCFC(K = 2, robust = FALSE, optns = plain_optns),
               "Supply either")
  expect_error(RKCFC(Ly = d$Ly, Lt = d$Lt, data = data.frame(), K = 2,
                     robust = FALSE, optns = plain_optns), "Supply either")
  expect_error(RKCFC(Ly = d$Ly, Lt = d$Lt, K = 2.5, robust = FALSE,
                     optns = plain_optns), "positive integer")
  expect_error(RKCFC(Ly = d$Ly, Lt = d$Lt, K = 100, robust = FALSE,
                     optns = plain_optns), "cannot exceed")
})

test_that("silhouette_score returns one value per subject and a scalar summary", {
  skip_on_cran()
  d <- sim_two_class(n = 40)
  long <- data.frame(
    subject = rep(seq_along(d$Ly), lengths(d$Ly)),
    t = unlist(d$Lt), value = unlist(d$Ly),
    g = rep(d$g, lengths(d$Ly))
  )
  s <- silhouette_score(long, class = "g", id = "subject", y = "value",
                        time = "t", robust = TRUE, option = robust_optns)
  expect_length(s$a_i, 40L)
  expect_length(s$b_i, 40L)
  expect_length(s$silhouette, 1L)
  expect_true(is.finite(s$silhouette))
})
