# test-shiny_app.R — Shiny helper and result table tests.

load_shiny_app_environment <- function() {
  app_path <- system.file("shiny", "app.R", package = "CorrSurvGSD")
  expect_true(nzchar(app_path))
  expect_true(file.exists(app_path))

  app_environment <- new.env(parent = globalenv())
  sys.source(app_path, envir = app_environment)
  app_environment
}

complete_shiny_input <- function() {
  list(
    n_T = 100,
    n_C = 100,
    median_PFS_C = 4,
    median_OS_C = 8,
    effect_mode = "hr",
    HR_PFS = 0.8,
    HR_OS = 0.85,
    accrual_segments = 1,
    accrual_rate_1 = 20,
    d_PFS_vec = "50, 100",
    alpha_spending_PFS = "OF",
    alpha_spending_OS = "Pocock"
  )
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

test_that("Shiny required-input guidance checks completeness only", {
  skip_if_not_installed("shiny")
  app <- load_shiny_app_environment()
  input <- complete_shiny_input()

  expect_false(app$is_missing_input(1))
  expect_true(app$is_missing_input(NULL))
  expect_true(app$is_missing_input("   "))
  expect_length(app$required_input_labels(input), 0)

  input$n_T <- NULL
  expect_identical(
    app$required_input_labels(input),
    list(n_T = "Treatment group sample size")
  )

  input <- complete_shiny_input()
  input$effect_mode <- "median"
  input$median_PFS_T <- NULL
  input$median_OS_T <- NULL
  expect_identical(
    app$required_input_labels(input),
    list(
      median_PFS_T = "Treatment group median PFS",
      median_OS_T = "Treatment group median OS"
    )
  )

  input <- complete_shiny_input()
  input$alpha_spending_PFS <- "HSD"
  input$alpha_spending_gamma_PFS <- NULL
  expect_identical(
    app$required_input_labels(input),
    list(alpha_spending_gamma_PFS = "PFS HSD gamma")
  )
})

test_that("Shiny translates pipeline errors for users", {
  skip_if_not_installed("shiny")
  app <- load_shiny_app_environment()

  model_error <- app$translate_pipeline_error(
    "lambda_1_T must be positive after applying HR_PFS and HR_OS."
  )
  expect_match(model_error, "incompatible with the illness-to-death model")
  expect_false(grepl("lambda_1_T", model_error, fixed = TRUE))
  expect_identical(
    app$pipeline_error_input_ids("lambda_1_T must be positive"),
    c("HR_PFS", "HR_OS", "median_PFS_T", "median_OS_T")
  )

  vector_error <- app$translate_pipeline_error(
    "d_PFS_vec must be a comma-separated numeric vector."
  )
  expect_match(vector_error, "comma-separated numbers")
  expect_false(grepl("d_PFS_vec", vector_error, fixed = TRUE))
  expect_identical(
    app$pipeline_error_input_ids("d_PFS_vec must be invalid"),
    "d_PFS_vec"
  )

  expect_identical(
    app$translate_pipeline_error("unmapped internal failure"),
    paste0(
      "The trial design could not be evaluated. Check the entered inputs ",
      "and try again."
    )
  )
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
  boundary_z <- app$boundary_results_table(state, scale = "z")
  boundary_hr <- app$boundary_results_table(state, scale = "hr")
  boundary_p <- app$boundary_results_table(state, scale = "p")
  framework <- app$testing_framework_details(state)
  marginal <- app$marginal_power_results_table(state)
  joint <- app$joint_power_results_table(
    state$theoretical_results$joint_power_matrix
  )
  context <- app$result_context_text(state)

  expect_equal(nrow(trial_design), 3)
  expect_equal(nrow(accrual), 1)
  expect_equal(nrow(schedule), 2)
  expect_equal(nrow(testing), 4)
  expect_equal(nrow(boundary_z), 4)
  expect_equal(nrow(framework), 4)
  expect_equal(nrow(marginal), 4)
  expect_equal(nrow(joint), 2)
  expect_named(
    trial_design,
    c("Characteristic", "Treatment", "Control", "Comparative measure")
  )
  expect_named(schedule, c("Look", "Target PFS events", "Derived OS events",
                           "Calendar cutoff"))
  expect_named(
    boundary_z,
    c(
      "Endpoint",
      "Look",
      "Futility HR",
      "Futility boundary",
      "Efficacy boundary"
    )
  )
  expect_false(identical(boundary_z, boundary_hr))
  expect_false(identical(boundary_hr, boundary_p))
  expect_match(trial_design$`Comparative measure`[2], "PFS HR (T/C)",
               fixed = TRUE)
  expect_identical(framework$Label[2], "Gatekeeping order")
  expect_match(context, "200 participants")
  expect_identical(names(joint)[1], "Endpoint, look")
  expect_identical(joint[[1]][1], "PFS, look 1")
  expect_identical(names(joint)[2], "OS, look 1")
  expect_true(all(vapply(joint, is.character, logical(1))))
})

test_that("Shiny power tables place simulation values in parentheses", {
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
    alpha_spending_OS = "OF"
  )
  state$options$simulation <- TRUE
  state$empirical_results <- list(
    marginal_power = state$theoretical_results$marginal_power,
    joint_power_matrix = state$theoretical_results$joint_power_matrix
  )

  marginal <- app$marginal_power_results_table(state)
  joint <- app$joint_power_results_table(
    state$theoretical_results$joint_power_matrix,
    simulation_power = state$empirical_results$joint_power_matrix
  )

  expect_named(
    marginal,
    c("Endpoint", "Look", "Incremental power", "Cumulative power")
  )
  expect_match(marginal$`Incremental power`[1], "\\(")
  expect_match(joint[[2]][1], "\\(")
})

test_that("Shiny renders the Results view after a successful pipeline run", {
  skip_if_not_installed("shiny")
  app <- load_shiny_app_environment()

  shiny::testServer(app$server, {
    session$setInputs(
      n_T = 100,
      n_C = 100,
      median_PFS_C = 4,
      median_OS_C = 8,
      effect_mode = "hr",
      HR_PFS = 0.8,
      HR_OS = 0.85,
      accrual_segments = 1,
      accrual_rate_1 = 20,
      d_PFS_vec = "50, 100",
      alpha_spending_PFS = "OF",
      alpha_spending_OS = "Pocock"
    )
    session$setInputs(run = 1)

    expect_match(output$status$html, "Pipeline completed")
    expect_match(output$results_panel$html, "Design summary")
    expect_match(output$results_panel$html, "Joint power")

    session$setInputs(n_T = 120)

    expect_match(output$status$html, "Pipeline completed")
    expect_match(output$result_context$html, "200 participants")
  })
})
