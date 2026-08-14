# test-run_simulation.R — simulation tests.
#
# Run with \`devtools::test()\` from the package root.

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

test_that("trial-level simulation is reproducible and finite", {
  first <- simulate_trial_statistics(
    design = simulation_design,
    n_sim = 6,
    seed = 123
  )
  second <- simulate_trial_statistics(
    design = simulation_design,
    n_sim = 6,
    seed = 123
  )

  expect_equal(dim(first), c(6, 4))
  expect_equal(colnames(first), c("PFS_1", "PFS_2", "OS_1", "OS_2"))
  expect_true(all(is.finite(first)))
  expect_equal(first, second)
})

test_that("simulation reports exhausted invalid attempts", {
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

test_that("first efficacy crossing applies futility and gate start", {
  futility_boundary <- rbind(
    c(futility = -1, efficacy = Inf),
    c(futility = -1, efficacy = 2)
  )
  gated_boundary <- rbind(
    c(futility = -Inf, efficacy = 2),
    c(futility = -Inf, efficacy = 2)
  )

  expect_identical(
    first_efficacy_crossing(
      c(-1.1, 3),
      futility_boundary,
      gate_start = 2
    ),
    NA_integer_
  )
  expect_identical(
    first_efficacy_crossing(
      c(3, 2.1),
      gated_boundary,
      gate_start = 2
    ),
    2L
  )
})

test_that("simulation power supports both hierarchy orders", {
  z <- rbind(
    c(1.5, 2.0, 1.5, 2.0),
    c(0.0, 1.5, 1.5, 1.5),
    c(0.0, 0.0, 1.5, 1.5),
    c(1.5, 0.0, 0.0, 1.5)
  )

  pfs_primary <- summarize_simulation_power(
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
  expected_pfs_joint <- matrix(
    c(0.25, 0.25, NA_real_, 0.25),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("PFS_1", "PFS_2"), c("OS_1", "OS_2"))
  )

  expect_equal(pfs_primary$marginal_power, expected_marginal)
  expect_equal(pfs_primary$joint_power_matrix, expected_pfs_joint)
  expect_equal(pfs_primary$joint_power, 0.75)

  os_primary <- summarize_simulation_power(
    z = z,
    boundary = simulation_boundary,
    hierarchy_order = c(primary = "OS", secondary = "PFS")
  )
  expected_os_joint <- matrix(
    c(0.25, 0.25, NA_real_, 0),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("OS_1", "OS_2"), c("PFS_1", "PFS_2"))
  )

  expect_equal(os_primary$joint_power_matrix, expected_os_joint)
  expect_equal(os_primary$joint_power, 0.50)
})

test_that("run_simulation preserves theoretical state and stores empirical results", {
  state <- list(
    design = c(
      simulation_design,
      list(
        boundary = simulation_boundary,
        hierarchy_order = c(primary = "PFS", secondary = "OS")
      )
    ),
    options = list(n_sim = 6, simulation_seed = 456),
    theoretical_results = list(marker = "unchanged")
  )
  result <- run_simulation(state)
  repeated <- run_simulation(state)

  expect_identical(result$design, state$design)
  expect_identical(result$theoretical_results, state$theoretical_results)
  expect_identical(result$empirical_results, repeated$empirical_results)
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
  expect_equal(
    dim(result$empirical_results$joint_correlation_matrix),
    c(4, 4)
  )
  expect_equal(dim(result$empirical_results$marginal_power), c(4, 2))
  expect_equal(dim(result$empirical_results$joint_power_matrix), c(2, 2))
  expect_equal(
    result$empirical_results$joint_power,
    sum(result$empirical_results$joint_power_matrix, na.rm = TRUE)
  )
})
