# test_utils.R — unit tests for utils.R
#
# Run from the package root with: Rscript tests/test_utils.R

source("R/utils.R")
library(testthat)

test_that("uniform accrual (s=1) validates correctly", {
  state <- validate_trial_input(
    n_T = 100, n_C = 100, d_PFS_vec = c(50),
    v_vec = c(200/12),
    median_PFS_C = log(2) / 0.15,
    median_OS_C = log(2) / 0.05,
    median_PFS_T = log(2) / 0.15,
    median_OS_T = log(2) / 0.05,
    alpha_spending_PFS = "OF",
    alpha_spending_OS = "Pocock"
  )
  design <- state$design
  expect_s3_class(state, "trial_state")
  expect_equal(design$n, 200)
  expect_equal(design$L, 1)
  expect_equal(design$R, 12)
  expect_equal(design$Lambda_T, 0.15)
  expect_identical(design$futility_analysis, c(PFS = FALSE, OS = FALSE))
  expect_true(all(is.na(design$futility_HR)))
  expect_identical(design$efficacy_start, c(PFS = 1L, OS = 1L))
  expect_identical(
    design$hierarchy_order,
    c(primary = "PFS", secondary = "OS")
  )
  expect_equal(state$options$tol, 1e-8)
  expect_equal(state$options$seed, 777)
})

test_that("HR input logic derives treatment hazards correctly", {
  state <- validate_trial_input(
    n_T = 100, n_C = 100, d_PFS_vec = c(50),
    v_vec = c(200/12),
    median_PFS_C = log(2) / 0.25,
    median_OS_C = log(2) / 0.10,
    HR_PFS = 0.6,
    HR_OS = 0.5,
    alpha_spending_PFS = "OF",
    alpha_spending_OS = "Pocock"
  )
  design <- state$design
  expect_equal(design$Lambda_T, 0.15)
  expect_equal(design$lambda_2_T, 0.05)
  expect_equal(design$lambda_1_T, 0.10)
  expect_equal(design$HR_PFS, 0.6)
  expect_equal(design$HR_OS, 0.5)
})

test_that("t_vec input derives R_vec correctly", {
  state <- validate_trial_input(
    n_T = 100, n_C = 100, d_PFS_vec = c(50),
    t_vec = c(6), v_vec = c(100/6, 100/6),
    median_PFS_C = log(2) / 0.15,
    median_OS_C = log(2) / 0.05,
    median_PFS_T = log(2) / 0.15,
    median_OS_T = log(2) / 0.05,
    alpha_spending_PFS = "OF",
    alpha_spending_OS = "Pocock"
  )
  design <- state$design
  expect_equal(design$R_vec, c(6, 6))
  expect_equal(design$R, 12)
})

test_that("t_vec non-NULL for s=1 throws error", {
  expect_error(
    validate_trial_input(
      n_T = 100, n_C = 100, d_PFS_vec = c(50),
      t_vec = c(6), v_vec = c(200/12),
      median_PFS_C = log(2) / 0.15,
      median_OS_C = log(2) / 0.05,
      median_PFS_T = log(2) / 0.15,
      median_OS_T = log(2) / 0.05,
      alpha_spending_PFS = "OF",
      alpha_spending_OS = "Pocock"
    ),
    "t_vec must be NULL for single-segment"
  )
})

test_that("v_vec length mismatch with t_vec throws error", {
  expect_error(
    validate_trial_input(
      n_T = 100, n_C = 100, d_PFS_vec = c(50),
      t_vec = c(6), v_vec = c(100/6, 100/6, 100/6),
      median_PFS_C = log(2) / 0.15,
      median_OS_C = log(2) / 0.05,
      median_PFS_T = log(2) / 0.15,
      median_OS_T = log(2) / 0.05,
      alpha_spending_PFS = "OF",
      alpha_spending_OS = "Pocock"
    ),
    "t_vec must have length"
  )
})

test_that("non-strictly-increasing t_vec throws error", {
  expect_error(
    validate_trial_input(
      n_T = 100, n_C = 100, d_PFS_vec = c(50),
      t_vec = c(6, 6), v_vec = c(100/6, 100/6, 100/6),
      median_PFS_C = log(2) / 0.15,
      median_OS_C = log(2) / 0.05,
      median_PFS_T = log(2) / 0.15,
      median_OS_T = log(2) / 0.05,
      alpha_spending_PFS = "OF",
      alpha_spending_OS = "Pocock"
    ),
    "t_vec must be strictly increasing"
  )
})

test_that("piecewise s>1 fully derives all quantities", {
  state <- validate_trial_input(
    n_T = 80, n_C = 120, d_PFS_vec = c(50, 100),
    t_vec = c(6), v_vec = c(80/6, 120/6),
    median_PFS_C = log(2) / 0.25,
    median_OS_C = log(2) / 0.10,
    median_PFS_T = log(2) / 0.15,
    median_OS_T = log(2) / 0.05,
    alpha_spending_PFS = "OF",
    alpha_spending_OS = "Pocock"
  )
  design <- state$design
  expect_equal(design$n, 200)
  expect_equal(design$n_T, 80)
  expect_equal(design$n_C, 120)
  expect_equal(design$R, 12)
  expect_equal(design$R_vec, c(6, 6))
  expect_equal(design$t_vec, c(6, 12))
  expect_equal(design$v_vec, c(80/6, 120/6))
  expect_equal(design$p_vec, c(80/6/200, 120/6/200))
  expect_equal(design$lambda_1_T, 0.10)
  expect_equal(design$lambda_2_T, 0.05)
  expect_equal(design$lambda_1_C, 0.15)
  expect_equal(design$lambda_2_C, 0.10)
  expect_equal(design$Lambda_T, 0.15)
  expect_equal(design$Lambda_C, 0.25)
  expect_equal(design$HR_PFS, 0.15 / 0.25)
  expect_equal(design$HR_OS, 0.05 / 0.10)
  expect_equal(design$L, 2)
})

test_that("d_PFS_ell <= n is enforced", {
  expect_error(
    validate_trial_input(
      n_T = 100, n_C = 100, d_PFS_vec = c(250),
      v_vec = c(200/12),
      median_PFS_C = log(2) / 0.15,
      median_OS_C = log(2) / 0.05,
      median_PFS_T = log(2) / 0.15,
      median_OS_T = log(2) / 0.05,
      alpha_spending_PFS = "OF",
      alpha_spending_OS = "Pocock"
    ),
    "must not exceed total sample size"
  )
})

test_that("R_s <= 0 throws error", {
  expect_error(
    validate_trial_input(
      n_T = 100, n_C = 100, d_PFS_vec = c(50),
      t_vec = c(20), v_vec = c(200/12, 10),
      median_PFS_C = log(2) / 0.15,
      median_OS_C = log(2) / 0.05,
      median_PFS_T = log(2) / 0.15,
      median_OS_T = log(2) / 0.05,
      alpha_spending_PFS = "OF",
      alpha_spending_OS = "Pocock"
    ),
    "R_s must be positive"
  )
})

test_that("treatment input logic must not mix medians and HRs", {
  expect_error(
    validate_trial_input(
      n_T = 100, n_C = 100, d_PFS_vec = c(50),
      v_vec = c(200/12),
      median_PFS_C = log(2) / 0.15,
      median_OS_C = log(2) / 0.05,
      median_PFS_T = log(2) / 0.15,
      median_OS_T = log(2) / 0.05,
      HR_PFS = 0.8,
      HR_OS = 0.8,
      alpha_spending_PFS = "OF",
      alpha_spending_OS = "Pocock"
    ),
    "not both"
  )
})

test_that("HR input logic validates derived lambda_1_T", {
  expect_error(
    validate_trial_input(
      n_T = 100, n_C = 100, d_PFS_vec = c(50),
      v_vec = c(200/12),
      median_PFS_C = log(2) / 0.15,
      median_OS_C = log(2) / 0.05,
      HR_PFS = 0.2,
      HR_OS = 0.9,
      alpha_spending_PFS = "OF",
      alpha_spending_OS = "Pocock"
    ),
    "lambda_1_T must be positive"
  )
})

test_that("efficacy_start requires futility when it exceeds one", {
  expect_error(
    validate_trial_input(
      n_T = 100, n_C = 100, d_PFS_vec = c(50, 100),
      v_vec = c(200 / 12),
      median_PFS_C = log(2) / 0.25,
      median_OS_C = log(2) / 0.10,
      HR_PFS = 0.8,
      HR_OS = 0.8,
      alpha_spending_PFS = "OF",
      alpha_spending_OS = "Pocock",
      efficacy_start = c(PFS = 2, OS = 1)
    ),
    "efficacy_start must be 1 when futility is disabled"
  )
})
