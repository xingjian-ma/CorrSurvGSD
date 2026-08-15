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
    futility_looks = list(PFS = 1, OS = integer(0)),
    integration_seed = 777
  )

  expect_s3_class(result, "trial_state")
  expect_false(result$options$simulation)
  expect_equal(result$design$L, 2)
  expect_true(is.finite(result$design$boundary["PFS_1", "futility"]))
  expect_true(is.infinite(result$design$boundary["PFS_2", "futility"]))
  expect_true(all(is.infinite(result$design$boundary[c("OS_1", "OS_2"), "futility"])))
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

test_that("run_pipeline defaults to final efficacy without futility", {
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
    alpha_spending_OS = "Pocock"
  )

  expect_identical(
    result$design$efficacy_looks,
    list(PFS = 2L, OS = 2L)
  )
  expect_identical(
    result$design$futility_looks,
    list(PFS = integer(0), OS = integer(0))
  )
  expect_identical(
    result$design$futility_HR,
    list(PFS = numeric(0), OS = numeric(0))
  )
  expect_true(all(is.infinite(result$design$boundary[, "futility"])))
})

test_that("run_pipeline supports HR inputs and OS-primary simulation", {
  result <- run_pipeline(
    n_T = 100,
    n_C = 100,
    d_PFS_vec = 40,
    v_vec = c(200 / 12),
    median_PFS_C = log(2) / 0.25,
    median_OS_C = log(2) / 0.10,
    HR_PFS = 0.8,
    HR_OS = 0.85,
    alpha_spending_PFS = "OF",
    alpha_spending_OS = "OF",
    efficacy_looks = list(PFS = 1, OS = 1),
    futility_looks = list(PFS = integer(0), OS = integer(0)),
    hierarchy_order = c(primary = "OS", secondary = "PFS"),
    simulation = TRUE,
    simulation_seed = 456,
    n_sim = 5
  )

  expect_equal(result$design$HR_PFS, 0.8)
  expect_equal(result$design$HR_OS, 0.85)
  expect_equal(
    rownames(result$theoretical_results$joint_power_matrix),
    "OS_1"
  )
  expect_equal(
    colnames(result$theoretical_results$joint_power_matrix),
    "PFS_1"
  )
  expect_true(!is.null(result$empirical_results))
  expect_named(
    result$empirical_results,
    c(
      "joint_mean_vector",
      "joint_correlation_matrix",
      "marginal_power",
      "joint_power_matrix",
      "joint_power"
    )
  )
})

test_that("run_pipeline supports HSD alpha spending", {
  result <- run_pipeline(
    n_T = 100,
    n_C = 100,
    d_PFS_vec = c(50, 100),
    v_vec = c(200 / 12),
    median_PFS_C = log(2) / 0.25,
    median_OS_C = log(2) / 0.10,
    median_PFS_T = log(2) / 0.15,
    median_OS_T = log(2) / 0.05,
    alpha_spending_PFS = "HSD",
    alpha_spending_OS = "HSD",
    alpha_spending_gamma_PFS = -4,
    alpha_spending_gamma_OS = -4,
    efficacy_looks = list(PFS = c(1, 2), OS = c(1, 2)),
    futility_looks = list(PFS = integer(0), OS = integer(0))
  )

  expect_equal(result$design$alpha_spending_gamma_PFS, -4)
  expect_equal(result$design$alpha_spending_gamma_OS, -4)
  expect_true(all(is.finite(result$design$boundary[, "efficacy"])))
})
