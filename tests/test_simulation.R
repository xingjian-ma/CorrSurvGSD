# test_simulation.R — unit tests for simulation.R
#
# Run from the package root with: Rscript tests/test_simulation.R

source("R/simulation.R")
library(testthat)

simulation_design <- list(
  n_T = 100,
  n_C = 100,
  n = 200,
  L = 2,
  R_vec = 12,
  t_vec = 12,
  v_vec = 200 / 12,
  d_PFS_vec = c(40, 80),
  lambda_1_T = 0.10,
  lambda_2_T = 0.05,
  lambda_1_C = 0.15,
  lambda_2_C = 0.10
)

simulation_boundary <- rbind(
  c(futility = -Inf, efficacy = 1),
  c(futility = -Inf, efficacy = 1),
  c(futility = -Inf, efficacy = 1),
  c(futility = -Inf, efficacy = 1)
)
rownames(simulation_boundary) <- c("PFS_1", "PFS_2", "OS_1", "OS_2")

# =================================================================
# 1. Trial-level simulation statistics
# =================================================================

test_that("simulate_trial_statistics returns reproducible valid Z-scores", {
  z_first <- simulate_trial_statistics(
    design = simulation_design,
    n_sim = 6,
    seed = 123
  )
  z_second <- simulate_trial_statistics(
    design = simulation_design,
    n_sim = 6,
    seed = 123
  )

  expect_equal(dim(z_first), c(6, 4))
  expect_equal(colnames(z_first), c("PFS_1", "PFS_2", "OS_1", "OS_2"))
  expect_true(all(is.finite(z_first)))
  expect_equal(z_first, z_second)
})

test_that("simulate_trial_statistics stops after exhausted invalid attempts", {
  no_valid_design <- simulation_design
  no_valid_design$n_T <- 1
  no_valid_design$n_C <- 1
  no_valid_design$n <- 2
  no_valid_design$L <- 1
  no_valid_design$v_vec <- 2 / 12
  no_valid_design$d_PFS_vec <- 1

  expect_error(
    simulate_trial_statistics(
      design = no_valid_design,
      n_sim = 1,
      seed = 123,
      max_attempts = 1
    ),
    "Unable to generate n_sim valid replicates"
  )
})

# =================================================================
# 2. First efficacy crossing
# =================================================================

test_that("first_efficacy_crossing evaluates futility before efficacy", {
  boundary <- rbind(
    c(futility = -1, efficacy = Inf),
    c(futility = -1, efficacy = 2)
  )

  expect_identical(
    first_efficacy_crossing(c(-1.1, 3), boundary, gate_start = 2),
    NA_integer_
  )
})

test_that("first_efficacy_crossing respects gate_start", {
  boundary <- rbind(
    c(futility = -Inf, efficacy = 2),
    c(futility = -Inf, efficacy = 2)
  )

  expect_identical(
    first_efficacy_crossing(c(3, 2.1), boundary, gate_start = 2),
    2L
  )
})

# =================================================================
# 3. Marginal and gatekeeping joint power
# =================================================================

test_that("summarize_simulation_power computes PFS-primary power", {
  z <- rbind(
    c(1.5, 2.0, 1.5, 2.0),
    c(0.0, 1.5, 1.5, 1.5),
    c(0.0, 0.0, 1.5, 1.5),
    c(1.5, 0.0, 0.0, 1.5)
  )
  result <- summarize_simulation_power(
    z = z,
    boundary = simulation_boundary,
    hierarchy_order = c(primary = "PFS", secondary = "OS")
  )
  expected_marginal <- rbind(
    c(incremental = 0.50, cumulative = 0.50),
    c(incremental = 0.25, cumulative = 0.75),
    c(incremental = 0.75, cumulative = 0.75),
    c(incremental = 0.25, cumulative = 1.00)
  )
  rownames(expected_marginal) <- c("PFS_1", "PFS_2", "OS_1", "OS_2")
  expected_joint <- matrix(
    c(0.25, 0.25, NA_real_, 0.25),
    nrow = 2,
    byrow = TRUE
  )
  rownames(expected_joint) <- c("PFS_1", "PFS_2")
  colnames(expected_joint) <- c("OS_1", "OS_2")

  expect_equal(result$marginal_power, expected_marginal)
  expect_equal(result$joint_power_matrix, expected_joint)
  expect_equal(result$joint_power, 0.75)
})

test_that("summarize_simulation_power supports an OS-primary hierarchy", {
  z <- rbind(
    c(1.5, 2.0, 1.5, 2.0),
    c(0.0, 1.5, 1.5, 1.5),
    c(0.0, 0.0, 1.5, 1.5),
    c(1.5, 0.0, 0.0, 1.5)
  )
  result <- summarize_simulation_power(
    z = z,
    boundary = simulation_boundary,
    hierarchy_order = c(primary = "OS", secondary = "PFS")
  )

  expect_equal(rownames(result$joint_power_matrix), c("OS_1", "OS_2"))
  expect_equal(colnames(result$joint_power_matrix), c("PFS_1", "PFS_2"))
  expect_equal(result$joint_power_matrix[1, 1], 0.25)
  expect_equal(result$joint_power_matrix[1, 2], 0.25)
  expect_equal(result$joint_power_matrix[2, 2], 0)
  expect_true(is.na(result$joint_power_matrix[2, 1]))
  expect_equal(result$joint_power, 0.50)
})

# =================================================================
# 4. Module 5 state handler
# =================================================================

test_that("run_simulation stores empirical results without mutation", {
  state <- list(
    design = c(
      simulation_design,
      list(
        boundary = simulation_boundary,
        hierarchy_order = c(primary = "PFS", secondary = "OS")
      )
    ),
    options = list(n_sim = 6, seed = 456),
    theoretical_results = list(marker = "unchanged")
  )
  result <- run_simulation(state)

  expect_identical(result$design, state$design)
  expect_identical(result$theoretical_results, state$theoretical_results)
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
  expect_equal(length(result$empirical_results$joint_mean_vector), 4)
  expect_equal(dim(result$empirical_results$joint_correlation_matrix), c(4, 4))
  expect_equal(dim(result$empirical_results$marginal_power), c(4, 2))
  expect_equal(dim(result$empirical_results$joint_power_matrix), c(2, 2))
  expect_equal(
    result$empirical_results$joint_power,
    sum(result$empirical_results$joint_power_matrix, na.rm = TRUE)
  )
})
