# The objective trace and the multiple-start procedure. They test whether the
# two steps monotonically optimise a common objective, so what they assert is
# exactly what the documentation claims: the reassignment step descends, the
# re-estimation step is not guaranteed to, and the difference is measured rather
# than assumed.

sim_two_class <- function(n = 60, d = 6, seed = 1) {
  set.seed(seed)
  g <- rep(1:2, length.out = n)
  grid <- seq(0, 1, length.out = 101)
  Lt <- lapply(seq_len(n), function(i) sort(sample(grid, max(2L, stats::rpois(1, 7)))))
  Ly <- lapply(seq_len(n), function(i)
    2 * sin(4 * Lt[[i]]) + d * (g[i] == 2) +
      3 * stats::rnorm(1) * sqrt(2) * sin(2 * pi * Lt[[i]]) +
      stats::rnorm(length(Lt[[i]]), 0, 0.1))
  list(Ly = Ly, Lt = Lt, g = g)
}

plain_optns <- list(dataType = "Sparse", methodMuCovEst = "smooth",
                    nRegGrid = 21, methodSelectK = 2, useBinnedData = "OFF")

test_that("the reassignment step descends in every recorded iteration", {
  skip_on_cran()
  d <- sim_two_class()
  fit <- RKCFC(Ly = d$Ly, Lt = d$Lt, K = 2, robust = FALSE, optns = plain_optns,
               init.method = "two-staged", maxIter = 15, parallel = FALSE,
               trace = 0, LOO = FALSE)

  expect_s3_class(fit$objective, "data.frame")
  expect_true(nrow(fit$objective) >= 1L)
  expect_true(all(fit$objective$reassignment_descends))
  expect_true(all(fit$objective$J_after <= fit$objective$J_before + 1e-8))
  # first row has no predecessor, so the re-estimation change is undefined there
  expect_true(is.na(fit$objective$J_reestimation[1]))
  expect_true(all(is.finite(fit$objective$J_reestimation[-1])))
})

test_that("the final objective is positive, finite and LOO-invariant", {
  skip_on_cran()
  d <- sim_two_class()
  args <- list(Ly = d$Ly, Lt = d$Lt, K = 2, robust = FALSE, optns = plain_optns,
               init.method = "two-staged", maxIter = 15, parallel = FALSE,
               trace = 0)
  set.seed(11); a <- do.call(RKCFC, c(args, LOO = FALSE))
  set.seed(11); b <- do.call(RKCFC, c(args, LOO = TRUE))

  expect_true(is.finite(a$objective_final) && a$objective_final > 0)
  # objective_final is defined from ONE set of class fits with no leave-one-out
  # refitting, so two fits reaching the same partition must score identically
  # whatever LOO was set to. That is what makes restarts comparable.
  if (mclust::adjustedRandIndex(a$res_clust$g, b$res_clust$g) > 0.999) {
    expect_equal(a$objective_final, b$objective_final)
  }
})

test_that("stop_reason distinguishes convergence from a cycle and from the cap", {
  skip_on_cran()
  d <- sim_two_class()
  common <- list(Ly = d$Ly, Lt = d$Lt, K = 2, robust = FALSE, optns = plain_optns,
                 init.method = "two-staged", parallel = FALSE, trace = 0,
                 LOO = FALSE)

  set.seed(3)
  free <- do.call(RKCFC, c(common, maxIter = 15))
  expect_true(free$stop_reason %in% c("converged", "cycle", "maxIter", "small_cluster"))
  expect_equal(free$conv, free$stop_reason %in% c("converged", "cycle"))

  # A cap of 2 permits exactly one reassignment, so unless it converges
  # immediately the cap must be what stopped it.
  set.seed(3)
  capped <- do.call(RKCFC, c(common, maxIter = 2))
  expect_true(capped$stop_reason %in% c("converged", "maxIter"))
  if (capped$stop_reason == "maxIter") expect_false(capped$conv)

  expect_true(is.na(free$cycle_length) || free$stop_reason == "cycle")
})

test_that("K = 1 reports the trivial stop and still scores an objective", {
  skip_on_cran()
  d <- sim_two_class(n = 30)
  fit <- RKCFC(Ly = d$Ly, Lt = d$Lt, K = 1, robust = FALSE, optns = plain_optns,
               parallel = FALSE, trace = 0, LOO = FALSE)
  expect_equal(fit$stop_reason, "trivial")
  expect_null(fit$objective)
  expect_true(is.finite(fit$objective_final))
})

test_that("RKCFC_restarts returns the smallest-objective fit and summarises the rest", {
  skip_on_cran()
  d <- sim_two_class()
  out <- RKCFC_restarts(Ly = d$Ly, Lt = d$Lt, K = 2, robust = FALSE,
                        optns = plain_optns, nstart = 5, seed = 42,
                        maxIter = 10, parallel = FALSE, trace = 0, LOO = FALSE)

  expect_equal(nrow(out$starts), 5L)
  expect_equal(out$starts$seed, 42L + 1:5)
  ok <- !is.na(out$starts$J)
  expect_equal(out$fit$objective_final, min(out$starts$J[ok]))
  expect_equal(max(out$starts$aRI_to_best[ok]), 1)
  expect_true(out$prop_best > 0 && out$prop_best <= 1)
  expect_equal(dim(out$agreement), c(sum(ok), sum(ok)))
  expect_equal(diag(out$agreement), rep(1, sum(ok)), ignore_attr = TRUE)
})

test_that("RKCFC_restarts is reproducible under a seed and records failures", {
  skip_on_cran()
  d <- sim_two_class()
  args <- list(Ly = d$Ly, Lt = d$Lt, K = 2, robust = FALSE, optns = plain_optns,
               nstart = 3, seed = 7, maxIter = 8, parallel = FALSE, trace = 0,
               LOO = FALSE)
  a <- do.call(RKCFC_restarts, args)
  b <- do.call(RKCFC_restarts, args)
  expect_equal(a$starts$J, b$starts$J)
  expect_equal(a$fit$res_clust$g, b$fit$res_clust$g)
  expect_equal(a$n_failed, 0L)
  expect_true(all(is.na(a$starts$error)))
})

test_that("RKCFC_restarts validates nstart", {
  expect_error(RKCFC_restarts(nstart = 0), "positive integer")
  expect_error(RKCFC_restarts(nstart = 2.5), "positive integer")
})
