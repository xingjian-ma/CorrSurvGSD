# test-calculate_closed_form_results.R — closed-form result tests.
#
# Run with `devtools::test()` from the package root.

state <- make_test_state(
  stage = "closed",
  n_T = 100,
  n_C = 100,
  d_PFS_vec = c(50, 100),
  v_vec = c(200 / 12),
  median_PFS_C = log(2) / 0.25,
  median_OS_C = log(2) / 0.10,
  HR_PFS = 0.8,
  HR_OS = 0.85,
  alpha_spending_PFS = "OF",
  alpha_spending_OS = "OF",
  efficacy_looks = list(PFS = c(1, 2), OS = c(1, 2)),
  futility_looks = list(PFS = 1, OS = integer(0)),
  futility_HR = list(PFS = 1.2, OS = numeric(0))
)

result <- state

# =================================================================
# 1. Boundary and drift helpers
# =================================================================

test_that("build_boundary handles a single active efficacy look", {
  actual <- build_boundary(
    d_vec = 100,
    L = 1,
    alpha = 0.025,
    alpha_spending = "OF",
    alpha_spending_gamma = NA_real_,
    efficacy_looks = 1,
    futility_looks = integer(0),
    futility_HR = numeric(0),
    r = 1
  )
  expect_equal(
    unname(actual[1, "efficacy"]),
    unname(stats::qnorm(0.975))
  )
  expect_true(is.infinite(actual[1, "futility"]))
})

test_that("build_boundary retains futility before efficacy starts", {
  actual <- build_boundary(
    d_vec = c(50, 100),
    L = 2,
    alpha = 0.025,
    alpha_spending = "OF",
    alpha_spending_gamma = NA_real_,
    efficacy_looks = 2,
    futility_looks = 1,
    futility_HR = 1.2,
    r = 1
  )
  expect_true(is.infinite(actual[1, "efficacy"]))
  expect_true(is.finite(actual[1, "futility"]))
  expect_true(is.finite(actual[2, "efficacy"]))
})

test_that("mean_drift matches the closed-form expression", {
  actual <- mean_drift(0.8, 1, c(50, 100))
  expected <- -log(0.8) * sqrt(c(50, 100) / 4)
  expect_equal(actual, expected)
})

# =================================================================
# 2. First-crossing and marginal power
# =================================================================

test_that("first_crossing_bounds uses efficacy at the current look", {
  boundary <- cbind(
    futility = c(-0.5, -0.4),
    efficacy = c(2.5, 2.0)
  )
  bounds <- first_crossing_bounds(boundary, look = 2)
  expect_equal(bounds$lower, c(-0.5, 2.0))
  expect_equal(bounds$upper, c(2.5, Inf))
})

test_that("first_crossing_bounds delays efficacy to gate_start", {
  boundary <- cbind(
    futility = c(-1, -1, -1),
    efficacy = c(2, 2, 2)
  )
  bounds <- first_crossing_bounds(boundary, look = 3, gate_start = 2)

  expect_equal(bounds$lower, c(-1, -1, 2))
  expect_equal(bounds$upper, c(Inf, 2, Inf))
})

test_that("single-look marginal power matches the normal tail", {
  power <- calculate_marginal_power(
    mean_vector = 0,
    correlation_matrix = matrix(1, 1, 1),
    boundary = cbind(futility = -Inf, efficacy = 0),
    L = 1,
    seed = 777
  )
  expect_equal(unname(power[1, "incremental"]), 0.5)
  expect_equal(unname(power[1, "cumulative"]), 0.5)
})

test_that("marginal power is zero before efficacy_start", {
  power <- calculate_marginal_power(
    mean_vector = c(0, 0),
    correlation_matrix = diag(2),
    boundary = cbind(
      futility = c(-Inf, -Inf),
      efficacy = c(Inf, 0)
    ),
    L = 2,
    seed = 777
  )

  expect_equal(unname(power[1, "incremental"]), 0)
  expect_equal(unname(power[1, "cumulative"]), 0)
  expect_equal(unname(power[2, "incremental"]), 0.5)
  expect_equal(unname(power[2, "cumulative"]), 0.5)
})

# =================================================================
# 3. Joint power
# =================================================================

test_that("joint power supports OS-primary hierarchy", {
  boundary <- rbind(
    c(futility = -Inf, efficacy = 0),
    c(futility = -Inf, efficacy = 0),
    c(futility = -Inf, efficacy = 0),
    c(futility = -Inf, efficacy = 0)
  )
  joint_power_matrix <- calculate_joint_power(
    joint_mean_vector = rep(0, 4),
    joint_correlation_matrix = diag(4),
    boundary = boundary,
    hierarchy_order = c(primary = "OS", secondary = "PFS"),
    L = 2,
    seed = 777
  )

  expect_equal(rownames(joint_power_matrix), c("OS_1", "OS_2"))
  expect_equal(colnames(joint_power_matrix), c("PFS_1", "PFS_2"))
  expect_equal(unname(joint_power_matrix[1, 1]), 0.25)
  expect_equal(unname(joint_power_matrix[1, 2]), 0.125)
  expect_equal(unname(joint_power_matrix[2, 2]), 0.125)
  expect_true(is.na(joint_power_matrix[2, 1]))
  expect_equal(sum(joint_power_matrix, na.rm = TRUE), 0.5)
})

# =================================================================
# 4. Module 4 state handler
# =================================================================

test_that("calculate_closed_form_results appends all required results", {
  expected_mean <- c(
    mean_drift(0.8, 1, result$design$d_PFS_vec),
    mean_drift(0.85, 1, result$design$d_OS_vec)
  )
  expected_futility <- -log(1.2) *
    sqrt(result$design$d_PFS_vec * result$design$r / (1 + result$design$r)^2)

  expect_equal(dim(result$design$boundary), c(4, 2))
  expect_true(all(is.finite(result$design$boundary[, "efficacy"])))
  expect_equal(result$theoretical_results$joint_mean_vector, expected_mean)
  expect_equal(
    unname(result$design$boundary[1:2, "futility"]),
    c(expected_futility[1], -Inf)
  )
  expect_true(all(is.infinite(result$design$boundary[3:4, "futility"])))
  expect_equal(dim(result$theoretical_results$marginal_power), c(4, 2))
  expect_equal(rownames(result$theoretical_results$marginal_power),
               c("PFS_1", "PFS_2", "OS_1", "OS_2"))
  expect_equal(dim(result$theoretical_results$joint_power_matrix), c(2, 2))
  expect_true(is.na(result$theoretical_results$joint_power_matrix[2, 1]))
  expect_equal(
    result$theoretical_results$joint_power,
    sum(result$theoretical_results$joint_power_matrix, na.rm = TRUE)
  )
})
