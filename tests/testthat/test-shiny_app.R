# test-shiny_app.R — Shiny helper and result table tests.

load_shiny_app_environment <- function() {
  app_path <- system.file("shiny", "app.R", package = "CorrSurvGSD")
  expect_true(nzchar(app_path))
  expect_true(file.exists(app_path))

  app_environment <- new.env(parent = globalenv())
  sys.source(app_path, envir = app_environment)
  app_environment
}

test_that("Shiny parsing and formatting helpers handle edge cases", {
  skip_if_not_installed("shiny")
  app <- load_shiny_app_environment()

  expect_equal(
    app$parse_numeric_vector("1, 2, 3", "values"),
    c(1, 2, 3)
  )
  expect_null(app$parse_numeric_vector("", "values", allow_empty = TRUE))
  expect_error(
    app$parse_numeric_vector("1, invalid", "values"),
    "values must be a comma-separated numeric vector"
  )
  expect_equal(app$format_result_value(1.2345, digits = 2), "1.23")
  expect_equal(app$format_result_value(NA_real_), "-")
  expect_equal(app$format_boundary_value(Inf), "-")
})

test_that("Shiny result tables expose the expected structures", {
  skip_if_not_installed("shiny")
  app <- load_shiny_app_environment()
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

  trial_design <- app$trial_design_table(state)
  accrual <- app$accrual_design_table(state)
  schedule <- app$analysis_schedule_table(state)
  testing <- app$endpoint_testing_table(state)
  marginal <- app$marginal_power_results_table(state)
  joint <- app$joint_power_results_table(
    state$theoretical_results$joint_power_matrix
  )

  expect_equal(nrow(trial_design), 9)
  expect_equal(nrow(accrual), 1)
  expect_equal(nrow(schedule), 2)
  expect_equal(nrow(testing), 4)
  expect_equal(nrow(marginal), 4)
  expect_equal(nrow(joint), 2)
  expect_named(trial_design, c("Item", "Value"))
  expect_named(schedule, c("Look", "Target PFS events", "Derived OS events",
                           "Calendar cutoff"))
  expect_true(all(vapply(joint, is.character, logical(1))))
})
