test_that("the covariance solver returns a direct solution", {
  result <- RKCFC:::.rkcfc_solve_covariance(diag(2), c(3, -1))
  expect_equal(as.vector(result), c(3, -1))
})

test_that("the covariance solver handles a singular matrix", {
  result <- RKCFC:::.rkcfc_solve_covariance(matrix(1, 2, 2), c(2, 2))
  expect_equal(as.vector(result), c(1, 1), tolerance = 1e-7)
})

test_that("non-robust scores use the model covariance", {
  score <- RKCFC:::.rkcfc_nonrobust_scores(diag(2), c(4, 1), 1, c(2, 2))
  expect_equal(score, c(8 / 5, 1))
})

test_that("input validation rejects mismatched lists", {
  expect_equal(univar_check_input(list(1:3, 1:2), list(1:3, 1:2)), 2L)
  expect_error(univar_check_input(list(1:3), list(1:2)), "same length")
  expect_error(univar_check_input(list(1:3), list(1:3, 1:3)), "same number")
  expect_error(univar_check_input(list(), list()), "non-empty")
})

test_that("cluster helpers enforce the minimum class size", {
  labels <- c(rep(1L, 12), rep(2L, 3))
  expect_length(RKCFC:::.rkcfc_cluster_ids(labels, 2L), 2L)
  expect_error(RKCFC:::.rkcfc_require_cluster_sizes(labels, 2L), "at least 10")
  expect_silent(RKCFC:::.rkcfc_require_cluster_sizes(rep(1:2, each = 10), 2L))
})

test_that("parallel worker validation is strict about nmcores", {
  expect_equal(RKCFC:::.rkcfc_validate_parallel_workers(FALSE, NULL, 100), 1L)
  expect_equal(RKCFC:::.rkcfc_validate_parallel_workers(TRUE, 4, 100), 4L)
  expect_equal(RKCFC:::.rkcfc_validate_parallel_workers(TRUE, 40, 10), 10L)
  expect_error(RKCFC:::.rkcfc_validate_parallel_workers(TRUE, NULL, 10), "nmcores")
})
