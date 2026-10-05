test_that("robust_scatter orders components by the recovered scale", {
  # The sign covariance matrix has trace 1, so its eigenvalues are not the
  # score spreads. Ordering by them silently returns the wrong subspace when
  # truncating; ordering by the recovered scale does not.
  set.seed(1)
  Xi <- matrix(stats::rt(3000, df = 2), ncol = 3) %*% diag(c(5, 2, 1))
  values <- robust_scatter(Xi)$values
  expect_equal(values, sort(values, decreasing = TRUE))
  expect_gt(values[1], values[3])
})

test_that("robust_scatter falls back to the moment estimator when asked", {
  set.seed(2)
  Xi <- matrix(stats::rnorm(3000), ncol = 3) %*% diag(c(5, 2, 1))
  moment <- robust_scatter(Xi, robust = FALSE)
  expect_equal(moment$values, sort(moment$values, decreasing = TRUE))
  # On Gaussian data the two estimators should recover the same subspace.
  expect_lt(max(principal_angles(moment$vectors, robust_scatter(Xi)$vectors)), 10)
})

test_that("whitening turns Mahalanobis distance into Euclidean distance", {
  set.seed(3)
  X <- matrix(stats::rnorm(2000), ncol = 2) %*% matrix(c(3, 1, 0, 2), 2)
  W <- whiten_scores(X, robust = FALSE)
  # covariance of the whitened coordinates is approximately the identity
  expect_equal(stats::cov(W), diag(2), tolerance = 0.1)
})

test_that("ch_index keeps its degrees-of-freedom divisors", {
  # Without (K-1) and (n-K) the bare B/W ratio is monotone in K and its argmax
  # is always the largest K offered, which would make the criterion useless.
  # Here K = 2 is the truth and must win.
  set.seed(4)
  X <- rbind(matrix(stats::rnorm(400, -3), ncol = 2),
             matrix(stats::rnorm(400,  3), ncol = 2))
  parts <- lapply(2:6, function(k) stats::kmeans(X, centers = k, nstart = 5)$cluster)
  ch <- vapply(parts, ch_index, numeric(1), scores = X, robust = FALSE)
  expect_equal(which.max(ch), 1L)
})

test_that("silhouette_index is undefined below two clusters and agrees broadly across methods", {
  set.seed(5)
  X <- rbind(matrix(stats::rnorm(300, -3), ncol = 2),
             matrix(stats::rnorm(300,  3), ncol = 2))
  g <- rep(1:2, each = 150)
  expect_true(is.na(silhouette_index(X, rep(1L, 300))))
  pairwise <- silhouette_index(X, g, robust = FALSE, method = "pairwise")
  centroid <- silhouette_index(X, g, robust = FALSE, method = "centroid")
  expect_gt(pairwise, 0.5)
  expect_gt(centroid, 0.5)
})

test_that("select_K returns one row per partition with both criteria", {
  set.seed(6)
  X <- rbind(matrix(stats::rnorm(300, -3), ncol = 2),
             matrix(stats::rnorm(300,  3), ncol = 2))
  parts <- lapply(2:4, function(k) stats::kmeans(X, centers = k)$cluster)
  out <- select_K(X, parts, robust = FALSE)
  expect_equal(out$K, 2:4)
  expect_true(all(is.finite(out$CH)))
  expect_true(all(is.finite(out$silhouette)))
})
