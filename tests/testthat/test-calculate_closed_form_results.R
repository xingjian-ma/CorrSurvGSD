# test-calculate_closed_form_results.R — closed-form result tests.
#
# Run with \`devtools::test()\` from the package root.

test_that("boundary helpers and mean drift preserve reference values", {
  single <- build_boundary(
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
  staged <- build_boundary(
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

  expect_equal(
    unname(single[1, "efficacy"]),
    stats::qnorm(0.975)
  )
  expect_true(is.infinite(single[1, "futility"]))
  expect_true(is.infinite(staged[1, "efficacy"]))
  expect_true(is.finite(staged[1, "futility"]))
  expect_true(is.finite(staged[2, "efficacy"]))
  expect_equal(
    mean_drift(0.8, 1, c(50, 100)),
    -log(0.8) * sqrt(c(50, 100) / 4)
  )
})

test_that("first-crossing bounds respect current look and gate start", {
  boundary <- cbind(
    futility = c(-1, -1, -1),
    efficacy = c(2, 2, 2)
  )
  current <- first_crossing_bounds(boundary, look = 3)
  gated <- first_crossing_bounds(boundary, look = 3, gate_start = 2)

  expect_equal(current$lower, c(-1, -1, 2))
  expect_equal(current$upper, c(2, 2, Inf))
  expect_equal(gated$lower, c(-1, -1, 2))
  expect_equal(gated$upper, c(Inf, 2, Inf))
})

test_that("marginal power handles single look and delayed efficacy", {
  single <- calculate_marginal_power(
    mean_vector = 0,
    correlation_matrix = matrix(1, 1, 1),
    boundary = cbind(futility = -Inf, efficacy = 0),
    L = 1,
    seed = 777
  )
  delayed <- calculate_marginal_power(
    mean_vector = c(0, 0),
    correlation_matrix = diag(2),
    boundary = cbind(
      futility = c(-Inf, -Inf),
      efficacy = c(Inf, 0)
    ),
    L = 2,
    seed = 777
  )

  expect_equal(unname(single[1, ]), c(0.5, 0.5))
  expect_equal(unname(delayed[1, ]), c(0, 0))
  expect_equal(unname(delayed[2, ]), c(0.5, 0.5))
})

test_that("joint power respects OS-primary hierarchy", {
  boundary <- rbind(
    c(futility = -Inf, efficacy = 0),
    c(futility = -Inf, efficacy = 0),
    c(futility = -Inf, efficacy = 0),
    c(futility = -Inf, efficacy = 0)
  )
  result <- calculate_joint_power(
    joint_mean_vector = rep(0, 4),
    joint_correlation_matrix = diag(4),
    boundary = boundary,
    hierarchy_order = c(primary = "OS", secondary = "PFS"),
    L = 2,
    seed = 777
  )

  expect_equal(rownames(result), c("OS_1", "OS_2"))
  expect_equal(colnames(result), c("PFS_1", "PFS_2"))
  expect_equal(unname(result[1, 1]), 0.25)
  expect_equal(unname(result[1, 2]), 0.125)
  expect_equal(unname(result[2, 2]), 0.125)
  expect_true(is.na(result[2, 1]))
  expect_equal(sum(result, na.rm = TRUE), 0.5)
})

test_that("closed-form state contains boundaries and power results", {
  result <- make_test_state(
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
  design <- result$design
  theoretical <- result$theoretical_results
  finite_boundary <- is.finite(design$boundary)

  expect_equal(dim(design$boundary), c(4, 2))
  expect_true(all(is.finite(design$boundary[, "efficacy"])))
  expect_equal(dim(design$boundary_HR), c(4, 2))
  expect_equal(dim(design$boundary_p), c(4, 2))
  expect_equal(dimnames(design$boundary_HR), dimnames(design$boundary))
  expect_equal(dimnames(design$boundary_p), dimnames(design$boundary))
  expect_true(all(is.finite(design$boundary_HR[finite_boundary])))
  expect_true(all(is.finite(design$boundary_p[finite_boundary])))
  expect_true(all(is.na(design$boundary_HR[!finite_boundary])))
  expect_true(all(is.na(design$boundary_p[!finite_boundary])))

  efficacy_scale <- sqrt(
    design$d_PFS_vec[1] * design$r / (1 + design$r)^2
  )
  expect_equal(
    design$boundary_HR[1, "efficacy"],
    exp(-design$boundary[1, "efficacy"] / efficacy_scale)
  )
  expect_equal(
    design$boundary_p[1, "efficacy"],
    stats::pnorm(design$boundary[1, "efficacy"], lower.tail = FALSE)
  )

  expected_mean <- c(
    mean_drift(0.8, 1, design$d_PFS_vec),
    mean_drift(0.85, 1, design$d_OS_vec)
  )
  expect_equal(theoretical$joint_mean_vector, expected_mean)
  expected_futility <- -log(1.2) * sqrt(
    design$d_PFS_vec * design$r / (1 + design$r)^2
  )
  expect_equal(
    unname(design$boundary[1:2, "futility"]),
    c(expected_futility[1], -Inf)
  )
  expect_true(all(is.infinite(design$boundary[3:4, "futility"])))
  expect_equal(dim(theoretical$marginal_power), c(4, 2))
  expect_equal(
    rownames(theoretical$marginal_power),
    c("PFS_1", "PFS_2", "OS_1", "OS_2")
  )
  expect_equal(dim(theoretical$joint_power_matrix), c(2, 2))
  expect_true(is.na(theoretical$joint_power_matrix[2, 1]))
  expect_equal(
    theoretical$joint_power,
    sum(theoretical$joint_power_matrix, na.rm = TRUE)
  )
})
