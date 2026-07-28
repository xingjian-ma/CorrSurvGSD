# test_calendar_cutoff.R — unit tests for calendar_cutoff.R
#
# Run from the package root with: Rscript tests/test_calendar_cutoff.R

source("R/utils.R")
source("R/calendar_cutoff.R")
library(testthat)

# -----------------------------------------------------------------
# pfs_mean
# -----------------------------------------------------------------
test_that("pfs_mean at A = 0 is zero", {
  expect_equal(pfs_mean(0, 12, c(6, 12), c(1/12, 1/12), 0.15), 0)
})

test_that("pfs_mean increases with A", {
  A_seq <- seq(0, 24, length.out = 10)
  vals  <- vapply(A_seq, function(a) pfs_mean(a, 12, c(6, 12), c(1/12, 1/12), 0.15), 0)
  expect_true(all(diff(vals) >= 0))
})

test_that("pfs_mean approaches 1 for large A", {
  expect_gt(pfs_mean(120, 12, c(6, 12), c(1/12, 1/12), 0.15), 0.99)
})

test_that("pfs_mean hand-checked value", {
  # A = R = 12, B = 12, Lambda = 0.15
  # I1 = 1, I4 = (1 - exp(-1.8)) / 1.8
  # m_delta = 1 - I4 ≈ 0.5363
  actual   <- pfs_mean(12, 12, c(6, 12), c(1/12, 1/12), 0.15)
  expected <- 1 - (1 - exp(-1.8)) / 1.8
  expect_equal(actual, expected, tolerance = 1e-6)
})

# -----------------------------------------------------------------
# pfs_event_diff
# -----------------------------------------------------------------
test_that("pfs_event_diff is negative at A = 0", {
  diff0 <- pfs_event_diff(0, 100, 100, 0.15, 0.15, 12, c(6, 12), c(1/12, 1/12), 50)
  expect_lt(diff0, 0)
})

test_that("pfs_event_diff is positive for large A", {
  diff_large <- pfs_event_diff(120, 100, 100, 0.15, 0.15, 12, c(6, 12), c(1/12, 1/12), 50)
  expect_gt(diff_large, 0)
})

test_that("pfs_event_diff is strictly increasing", {
  A_seq <- seq(0, 48, length.out = 20)
  diffs <- vapply(A_seq, function(a) {
    pfs_event_diff(a, 100, 100, 0.15, 0.15, 12, c(6, 12), c(1/12, 1/12), 50)
  }, 0)
  expect_true(all(diff(diffs) > 0))
})

# -----------------------------------------------------------------
# calendar_cutoff_single
# -----------------------------------------------------------------
test_that("calendar_cutoff_single returns a valid root", {
  A_root <- calendar_cutoff_single(100, 100, 0.15, 0.15,
                                   12, c(6, 12), c(1/12, 1/12), 50, 0)
  residual <- pfs_event_diff(A_root, 100, 100, 0.15, 0.15, 12, c(6, 12), c(1/12, 1/12), 50)
  expect_equal(residual, 0, tolerance = 1e-6)
})

test_that("calendar_cutoff_single respects A_prev", {
  A_root <- calendar_cutoff_single(100, 100, 0.15, 0.15,
                                   12, c(6, 12), c(1/12, 1/12), 120, A_prev = 10)
  expect_gt(A_root, 10)
})

test_that("calendar_cutoff_single extends bracket when needed", {
  # small R, large target → bracket needs to grow past 2R
  A_root <- calendar_cutoff_single(200, 200, 0.07, 0.07,
                                   6, c(3, 6), c(1/6, 1/6), 300, 0)
  expect_gt(A_root, 12)
  residual <- pfs_event_diff(A_root, 200, 200, 0.07, 0.07, 6, c(3, 6), c(1/6, 1/6), 300)
  expect_equal(residual, 0, tolerance = 1e-6)
})

# -----------------------------------------------------------------
# calendar_cutoff — full pipeline
# -----------------------------------------------------------------
test_that("calendar_cutoff returns non-decreasing A_vec", {
  state <- validate_trial_input(
    n_T = 100, n_C = 100, d_PFS_vec = c(30, 60, 90),
    t_vec = c(6), v_vec = c(200 / 12, 200 / 12),
    median_PFS_C = log(2) / 0.15,
    median_OS_C = log(2) / 0.05,
    median_PFS_T = log(2) / 0.15,
    median_OS_T = log(2) / 0.05,
    alpha_spending_PFS = "OF",
    alpha_spending_OS = "Pocock"
  )
  state <- calendar_cutoff(state)
  design <- state$design
  expect_true(all(diff(design$A_vec) >= 0))
  expect_length(design$A_vec, 3)
})

test_that("calendar_cutoff roots satisfy equation (7)", {
  state <- validate_trial_input(
    n_T = 150, n_C = 150, d_PFS_vec = c(50, 100, 150),
    t_vec = c(12), v_vec = c(300 / 24, 300 / 24),
    median_PFS_C = log(2) / 0.25,
    median_OS_C = log(2) / 0.10,
    median_PFS_T = log(2) / 0.15,
    median_OS_T = log(2) / 0.05,
    alpha_spending_PFS = "OF",
    alpha_spending_OS = "Pocock"
  )
  state <- calendar_cutoff(state)
  design <- state$design
  for (ell in seq_along(design$A_vec)) {
    residual <- pfs_event_diff(
      design$A_vec[ell], 150, 150, 0.15, 0.25,
      24, c(12, 24), c(1/24, 1/24), design$d_PFS_vec[ell]
    )
    expect_equal(residual, 0, tolerance = 1e-6)
  }
})
