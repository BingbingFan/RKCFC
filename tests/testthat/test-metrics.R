test_that("rre is 0 for a perfect fit and 1 for the mean", {
  set.seed(1)
  y <- stats::rnorm(200)
  expect_equal(rre(y, y), 0)
  expect_equal(rre(y, rep(mean(y), length(y))), 1)
  expect_gt(rre(y, rep(0, length(y)) + 5), 1)   # worse than the mean
})

test_that("rre trims on the observation, not on the error", {
  set.seed(2)
  y <- c(stats::rnorm(199), 50)                 # one extreme observation
  yhat <- y
  yhat[200] <- 0                                # all of the error is there
  # Untrimmed the single point dominates; trimmed it is excluded because the
  # OBSERVATION is extreme, which is the behaviour the trimmed value exists for.
  expect_gt(rre(y, yhat), rre(y, yhat, trim = 0.02))
  expect_lt(rre(y, yhat, trim = 0.02), 1e-8)
})

test_that("robust_normalize centres and scales, and round-trips", {
  set.seed(3)
  y <- c(stats::rlnorm(200), 1e5)               # one extreme value
  z <- robust_normalize(y)
  expect_equal(stats::median(z), 0)
  expect_equal(stats::IQR(z), 1)
  expect_equal(back_transform(z, attr(z, "centre"), attr(z, "scale")),
               y, ignore_attr = TRUE)
})

test_that("robust_normalize falls back only on a genuinely zero IQR", {
  # Half the observations at the floor is NOT enough: with type-7 quantiles Q3
  # is still above it, so the IQR stays positive and no fallback fires.
  half <- c(rep(0, 100), seq_len(100))
  expect_silent(robust_normalize(half))

  floored <- c(rep(0, 90), seq_len(10))         # more than three quarters
  expect_warning(z <- robust_normalize(floored), "IQR was zero")
  expect_true(is.finite(attr(z, "scale")) && attr(z, "scale") > 0)

  # A constant column: the fallback fires, still finds no spread, and stops.
  expect_error(suppressWarnings(robust_normalize(rep(1, 50))), "no usable spread")
})

test_that("the simulation arms are IQR-matched to the Gaussian arm", {
  set.seed(4)
  Z <- diag(c(9, 4))
  iqr_of <- function(arm) mean(apply(gen_scores(20000, arm, Z), 2, stats::IQR) /
                                 sqrt(diag(Z)))
  reference <- iqr_of(1)
  for (arm in 2:6) {
    expect_equal(iqr_of(arm), reference, tolerance = 0.06,
                 info = paste("arm", arm))
  }
})

test_that("the Fleishman arms hit their target moments", {
  skip_on_cran()
  set.seed(5)
  moments <- function(arm) {
    x <- as.vector(gen_scores(2e5, arm, diag(1)))
    x <- x / stats::sd(x)
    c(skew = mean(x^3), exkurt = mean(x^4) - 3)
  }
  expect_equal(unname(moments(2)), c(0, 4), tolerance = 0.15)
  expect_equal(unname(moments(3)), c(1.5, 4), tolerance = 0.2)
})

test_that("gen_scores rejects an unknown arm", {
  expect_error(gen_scores(10, 7, diag(2)), "arm must be 1..6")
})
