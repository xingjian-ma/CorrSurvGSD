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
  expect_match(model_error, "implied treatment median OS must be longer")
  expect_match(model_error, "Adjust the treatment hazard ratios")
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
    app$translate_pipeline_error(
      "median_PFS_C must be less than median_OS_C"
    ),
    "Control median OS must be longer than control median PFS."
  )
  expect_identical(
    app$translate_pipeline_error(
      "median_PFS_T must be less than median_OS_T"
    ),
    "Treatment median OS must be longer than treatment median PFS."
  )

  expect_identical(
    app$translate_pipeline_error("unmapped internal failure"),
    paste0(
      "The trial design could not be evaluated. Check the entered inputs ",
      "and try again."
    )
  )
})

test_that("Shiny explains treatment median order errors for HR inputs", {
  skip_if_not_installed("shiny")
  app <- load_shiny_app_environment()

  shiny::testServer(app$server, {
    session$setInputs(
      n_T = 100,
      n_C = 100,
      median_PFS_C = 4,
      median_OS_C = 8,
      effect_mode = "hr",
      HR_PFS = 0.1,
      HR_OS = 0.9,
      accrual_segments = 1,
      accrual_rate_1 = 20,
      d_PFS_vec = "50, 100",
      alpha_spending_PFS = "OF",
      alpha_spending_OS = "Pocock"
    )
    session$setInputs(run = 1)

    expect_match(
      output$status$html,
      "implied treatment median OS must be longer"
    )
  })
})

test_that("Shiny introduction provides a fixed-sequence demo design", {
  skip_if_not_installed("shiny")
  app <- load_shiny_app_environment()
  demo <- app$demo_input_values()
  ui_html <- htmltools::renderTags(app$ui)$html

  expect_identical(demo$n_T, 100)
  expect_identical(demo$n_C, 100)
  expect_identical(demo$accrual_rate_1, 20)
  expect_identical(demo$hierarchy_primary, "PFS")
  expect_identical(demo$d_PFS_vec, "50, 100")
  expect_match(ui_html, "Load demo values", fixed = TRUE)
  expect_match(ui_html, "fixed testing sequence", fixed = TRUE)
  expect_match(ui_html, "Testing order", fixed = TRUE)
  expect_match(
    ui_html,
    "Closed-form evaluation engine for sequential PFS and OS testing",
    fixed = TRUE
  )
  expect_match(ui_html, "Time to progression (TTP)", fixed = TRUE)
  expect_match(ui_html, "overall survival (OS)", fixed = TRUE)
  expect_match(
    ui_html,
    "Time to progression (TTP) and OS are modeled",
    fixed = TRUE
  )
  expect_false(grepl("OS (overall survival)", ui_html, fixed = TRUE))
  expect_match(ui_html, "independent exponential", fixed = TRUE)
  expect_match(ui_html, "smaller of TTP and OS", fixed = TRUE)
  expect_false(grepl("Fleischer", ui_html, fixed = TRUE))
  expect_false(grepl("illness-death", ui_html, fixed = TRUE))
  expect_match(ui_html, "2. Calculate", fixed = TRUE)
  expect_match(ui_html, "3. Evaluate", fixed = TRUE)
  expect_match(ui_html, "Set the main design inputs", fixed = TRUE)
  expect_false(grepl("Use the Configure step", ui_html, fixed = TRUE))
  expect_false(grepl("Fill the inputs", ui_html, fixed = TRUE))
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
  expect_equal(nrow(joint), 3)
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
  expect_identical(framework$Label[2], "Testing order")
  expect_match(context, "200 participants")
  expect_identical(names(joint)[1], "Endpoint: look")
  expect_identical(joint[[1]][1], "PFS: look 1")
  expect_identical(
    names(joint),
    c("Endpoint: look", "OS: look 1", "OS: look 2", "OS: all looks")
  )
  expect_identical(joint[[1]][3], "PFS: all looks")
  expect_true(all(vapply(joint, is.character, logical(1))))
  expect_identical(
    joint$`OS: all looks`[3],
    app$format_result_value(
      sum(state$theoretical_results$joint_power_matrix, na.rm = TRUE),
      state$options$display_digits
    )
  )

  reverse_joint <- app$joint_power_results_table(matrix(
    c(0.01, 0.02, NA_real_, 0.03),
    nrow = 2,
    dimnames = list(c("OS_1", "OS_2"), c("PFS_1", "PFS_2"))
  ))
  expect_identical(
    names(reverse_joint),
    c("Endpoint: look", "PFS: look 1", "PFS: look 2", "PFS: all looks")
  )
  expect_identical(reverse_joint[[1]][3], "OS: all looks")
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
  expect_match(joint$`OS: all looks`[1], "\\(")
  expect_match(joint$`OS: all looks`[3], "\\(")
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
    expect_match(output$result_context$html, "download_control")
    expect_match(
      output$download_control$html,
      "Download tables (CSV)",
      fixed = TRUE
    )

    session$setInputs(power_view = "joint")
    expect_match(output$joint_power_results, "OS: all looks")
    expect_match(output$joint_power_results, "PFS: all looks")

    session$setInputs(n_T = 120)

    expect_match(output$status$html, "Pipeline completed")
    expect_match(output$result_context$html, "200 participants")
  })
})
