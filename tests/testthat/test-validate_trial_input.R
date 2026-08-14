# test-validate_trial_input.R — input validation tests.
#
# Run with \`devtools::test()\` from the package root.

test_that("uniform accrual builds the expected default design", {
  state <- validate_test_input(
    n_T = 100,
    n_C = 100,
    d_PFS_vec = c(50),
    v_vec = c(200 / 12),
    median_PFS_C = log(2) / 0.15,
    median_OS_C = log(2) / 0.05,
    median_PFS_T = log(2) / 0.15,
    median_OS_T = log(2) / 0.05,
    alpha_spending_PFS = "OF",
    alpha_spending_OS = "Pocock"
  )
  design <- state$design

  expect_s3_class(state, "trial_state")
  expect_equal(c(design$n, design$L, design$R), c(200, 1, 12))
  expect_equal(design$Lambda_T, 0.15)
  expect_identical(
    design$efficacy_looks,
    list(PFS = 1L, OS = 1L)
  )
  expect_identical(
    design$futility_looks,
    list(PFS = integer(0), OS = integer(0))
  )
  expect_identical(
    design$futility_HR,
    list(PFS = numeric(0), OS = numeric(0))
  )
  expect_identical(
    design$hierarchy_order,
    c(primary = "PFS", secondary = "OS")
  )
  expect_true(is.na(design$alpha_spending_gamma_PFS))
  expect_equal(state$options[c("tol", "integration_seed", "simulation_seed")],
               list(tol = 1e-8, integration_seed = 1, simulation_seed = 1))
})

test_that("HR inputs derive treatment hazards and ratios", {
  state <- validate_test_input(
    n_T = 100,
    n_C = 100,
    d_PFS_vec = c(50),
    v_vec = c(200 / 12),
    median_PFS_C = log(2) / 0.25,
    median_OS_C = log(2) / 0.10,
    HR_PFS = 0.6,
    HR_OS = 0.5,
    alpha_spending_PFS = "OF",
    alpha_spending_OS = "Pocock"
  )
  design <- state$design

  expect_equal(
    c(design$Lambda_T, design$lambda_2_T, design$lambda_1_T),
    c(0.15, 0.05, 0.10)
  )
  expect_equal(c(design$HR_PFS, design$HR_OS), c(0.6, 0.5))
})

test_that("piecewise accrual derives segment quantities", {
  state <- validate_test_input(
    n_T = 80,
    n_C = 120,
    d_PFS_vec = c(50, 100),
    t_vec = c(6),
    v_vec = c(80 / 6, 120 / 6),
    median_PFS_C = log(2) / 0.25,
    median_OS_C = log(2) / 0.10,
    median_PFS_T = log(2) / 0.15,
    median_OS_T = log(2) / 0.05,
    alpha_spending_PFS = "OF",
    alpha_spending_OS = "Pocock"
  )
  design <- state$design

  expect_equal(design$n, 200)
  expect_equal(design$R_vec, c(6, 6))
  expect_equal(design$t_vec, c(6, 12))
  expect_equal(design$v_vec, c(80 / 6, 120 / 6))
  expect_equal(design$p_vec, c(80 / 6 / 200, 120 / 6 / 200))
  expect_equal(
    c(design$lambda_1_T, design$lambda_2_T, design$Lambda_C),
    c(0.10, 0.05, 0.25)
  )
  expect_equal(design$L, 2)
})

test_that("three-segment accrual derives all segment durations", {
  state <- validate_test_input(
    n_T = 100,
    n_C = 100,
    d_PFS_vec = c(50, 100),
    t_vec = c(4, 8),
    v_vec = c(10, 15, 25),
    median_PFS_C = log(2) / 0.25,
    median_OS_C = log(2) / 0.10,
    median_PFS_T = log(2) / 0.15,
    median_OS_T = log(2) / 0.05,
    alpha_spending_PFS = "OF",
    alpha_spending_OS = "Pocock"
  )
  design <- state$design

  expect_equal(design$R_vec, c(4, 4, 4))
  expect_equal(design$t_vec, c(4, 8, 12))
  expect_equal(design$p_vec, c(0.05, 0.075, 0.125))
  expect_equal(design$R, 12)
})

test_that("valid futility HR and HSD specifications are normalized", {
  common <- list(
    n_T = 100,
    n_C = 100,
    d_PFS_vec = c(50, 100),
    v_vec = c(200 / 12),
    median_PFS_C = log(2) / 0.25,
    median_OS_C = log(2) / 0.10,
    HR_PFS = 0.8,
    HR_OS = 0.85,
    alpha_spending_PFS = "HSD",
    alpha_spending_OS = "HSD",
    alpha_spending_gamma_PFS = -40,
    alpha_spending_gamma_OS = 39.9,
    efficacy_looks = list(PFS = c(1, 2), OS = c(1, 2)),
    futility_looks = list(PFS = 1, OS = 1)
  )
  cases <- list(
    list(
      specification = 1.2,
      expected = list(PFS = 1.2, OS = 1.2)
    ),
    list(
      specification = c(PFS = 1.2, OS = 1.3),
      expected = list(PFS = 1.2, OS = 1.3)
    ),
    list(
      specification = list(PFS = 1.2, OS = 1.3),
      expected = list(PFS = 1.2, OS = 1.3)
    )
  )

  for (case in cases) {
    args <- c(common, list(futility_HR = case$specification))
    state <- do.call(validate_test_input, args)
    expect_equal(state$design$futility_HR, case$expected)
  }

  state <- do.call(validate_test_input, common)
  expect_equal(state$design$alpha_spending_gamma_PFS, -40)
  expect_equal(state$design$alpha_spending_gamma_OS, 39.9)
})

test_that("treatment and accrual input relationships are enforced", {
  expect_error(
    validate_test_input(
      median_PFS_T = log(2) / 0.15,
      median_OS_T = log(2) / 0.05,
      HR_PFS = 0.8,
      HR_OS = 0.8
    ),
    "not both"
  )
  expect_error(
    validate_test_input(
      median_PFS_T = log(2) / 0.05,
      median_OS_T = log(2) / 0.15
    ),
    "median_PFS_T must be less than median_OS_T"
  )
  expect_error(
    validate_test_input(
      HR_PFS = 0.2,
      HR_OS = 0.9
    ),
    "lambda_1_T must be positive"
  )
  expect_error(
    validate_test_input(
      median_PFS_T = NULL,
      median_OS_T = NULL,
      HR_PFS = NULL,
      HR_OS = NULL
    ),
    "Provide either treatment"
  )
  expect_error(
    validate_test_input(
      median_PFS_C = log(2) / 0.05,
      median_OS_C = log(2) / 0.10,
      median_PFS_T = log(2) / 0.15,
      median_OS_T = log(2) / 0.05
    ),
    "median_PFS_C must be less than median_OS_C"
  )
})

test_that("accrual and event-count constraints are enforced", {
  expect_error(
    validate_test_input(
      t_vec = c(6),
      v_vec = c(200 / 12),
      median_PFS_T = log(2) / 0.15,
      median_OS_T = log(2) / 0.05
    ),
    "t_vec must be NULL"
  )
  expect_error(
    validate_test_input(
      t_vec = c(6),
      v_vec = c(100 / 6, 100 / 6, 100 / 6),
      median_PFS_T = log(2) / 0.15,
      median_OS_T = log(2) / 0.05
    ),
    "t_vec must have length"
  )
  expect_error(
    validate_test_input(
      t_vec = c(6, 6),
      v_vec = c(100 / 6, 100 / 6, 100 / 6),
      median_PFS_T = log(2) / 0.15,
      median_OS_T = log(2) / 0.05
    ),
    "t_vec must be strictly increasing"
  )
  expect_error(
    validate_test_input(
      t_vec = c(20),
      v_vec = c(200 / 12, 10),
      median_PFS_T = log(2) / 0.15,
      median_OS_T = log(2) / 0.05
    ),
    "R_s must be positive"
  )
  expect_error(
    validate_test_input(
      d_PFS_vec = c(250),
      median_PFS_T = log(2) / 0.15,
      median_OS_T = log(2) / 0.05
    ),
    "must not exceed total sample size"
  )
  expect_error(
    validate_test_input(
      d_PFS_vec = c(50, 50),
      median_PFS_T = log(2) / 0.15,
      median_OS_T = log(2) / 0.05
    ),
    "strictly increasing"
  )
})

test_that("look schedules require valid endpoint structures", {
  expect_error(
    validate_test_input(
      d_PFS_vec = c(50, 100),
      efficacy_looks = list(PFS = 1, OS = c(1, 2)),
      median_PFS_T = log(2) / 0.15,
      median_OS_T = log(2) / 0.05
    ),
    "efficacy_looks must contain"
  )
  expect_error(
    validate_test_input(
      d_PFS_vec = c(50, 100),
      efficacy_looks = list(PFS = c(1, 2), OS = c(1, 2)),
      futility_looks = list(PFS = 2, OS = integer(0)),
      median_PFS_T = log(2) / 0.15,
      median_OS_T = log(2) / 0.05
    ),
    "futility_looks must contain"
  )
  expect_error(
    validate_test_input(
      d_PFS_vec = c(50, 100),
      efficacy_looks = list(PFS = 1, OS = c(1, 2), extra = 1),
      median_PFS_T = log(2) / 0.15,
      median_OS_T = log(2) / 0.05
    ),
    "named list"
  )
})

test_that("futility HR specifications reject malformed values", {
  expect_error(
    validate_test_input(
      d_PFS_vec = c(50, 100),
      futility_looks = list(PFS = 1, OS = 1),
      futility_HR = 1,
      median_PFS_T = log(2) / 0.15,
      median_OS_T = log(2) / 0.05
    ),
    "greater than 1"
  )
  expect_error(
    validate_test_input(
      d_PFS_vec = c(50, 100),
      futility_looks = list(PFS = 1, OS = 1),
      futility_HR = c(1.2, 1.3),
      median_PFS_T = log(2) / 0.15,
      median_OS_T = log(2) / 0.05
    ),
    "named PFS and OS"
  )
  expect_error(
    validate_test_input(
      d_PFS_vec = c(50, 100),
      futility_looks = list(PFS = 1, OS = integer(0)),
      futility_HR = list(PFS = c(1.2, 1.3), OS = numeric(0)),
      median_PFS_T = log(2) / 0.15,
      median_OS_T = log(2) / 0.05
    ),
    "match its futility_looks length"
  )
})

test_that("hierarchy, spending, and runtime options reject invalid values", {
  expect_error(
    validate_test_input(
      hierarchy_order = c(primary = "PFS", secondary = "PFS"),
      median_PFS_T = log(2) / 0.15,
      median_OS_T = log(2) / 0.05
    ),
    "hierarchy_order"
  )
  expect_error(
    validate_test_input(
      alpha_spending_PFS = "HSD",
      alpha_spending_gamma_PFS = 40,
      median_PFS_T = log(2) / 0.15,
      median_OS_T = log(2) / 0.05
    ),
    "interval \\[-40, 40\\)"
  )
  expect_error(
    validate_test_input(
      alpha = 0,
      median_PFS_T = log(2) / 0.15,
      median_OS_T = log(2) / 0.05
    ),
    "strictly between 0 and 1"
  )
  expect_error(
    validate_test_input(
      simulation = TRUE,
      n_sim = 1,
      median_PFS_T = log(2) / 0.15,
      median_OS_T = log(2) / 0.05
    ),
    "at least 2"
  )
  expect_error(
    validate_test_input(
      display_digits = 11,
      median_PFS_T = log(2) / 0.15,
      median_OS_T = log(2) / 0.05
    ),
    "at most 10"
  )
})
