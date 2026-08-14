# test-run_pipeline.R — public pipeline tests.

test_that("run_pipeline returns a complete theoretical result", {
  result <- run_pipeline(
    n_T = 100,
    n_C = 100,
    d_PFS_vec = c(50, 100),
    v_vec = c(200 / 12),
    median_PFS_C = log(2) / 0.25,
    median_OS_C = log(2) / 0.10,
    median_PFS_T = log(2) / 0.15,
    median_OS_T = log(2) / 0.05,
    alpha_spending_PFS = "OF",
    alpha_spending_OS = "Pocock",
    efficacy_looks = list(PFS = c(1, 2), OS = c(1, 2)),
    futility_looks = list(PFS = integer(0), OS = integer(0)),
    integration_seed = 777
  )

  expect_s3_class(result, "trial_state")
  expect_false(result$options$simulation)
  expect_equal(result$design$L, 2)
  expect_equal(result$design$efficacy_looks,
               list(PFS = c(1L, 2L), OS = c(1L, 2L)))
  expect_named(
    result$theoretical_results,
    c(
      "joint_correlation_matrix",
      "joint_mean_vector",
      "marginal_power",
      "joint_power_matrix",
      "joint_power"
    ),
    ignore.order = TRUE
  )
})
