# app.R — minimal interface for the complete analysis pipeline.

library(shiny)

parse_numeric_vector <- function(value, name, allow_empty = FALSE) {
  value <- trimws(value)

  if (allow_empty && !nzchar(value)) {
    return(NULL)
  }

  values <- suppressWarnings(
    as.numeric(trimws(strsplit(value, ",", fixed = TRUE)[[1]]))
  )

  if (length(values) == 0 || anyNA(values)) {
    stop(name, " must be a comma-separated numeric vector.")
  }

  values
}

format_result_value <- function(value, digits = 4) {
  if (is.na(value)) {
    return("—")
  }

  formatC(value, format = "f", digits = digits)
}

format_boundary_value <- function(value, digits = 4) {
  if (is.na(value)) {
    return("—")
  }
  if (is.infinite(value)) {
    return("-")
  }

  format_result_value(value, digits = digits)
}

marginal_power_results_table <- function(state) {
  theoretical <- state$theoretical_results$marginal_power
  digits <- state$options$display_digits
  simulation <- if (isTRUE(state$options$simulation)) {
    state$empirical_results$marginal_power
  } else {
    NULL
  }
  labels <- strsplit(rownames(theoretical), "_", fixed = TRUE)

  result <- data.frame(
    Endpoint = vapply(labels, `[[`, character(1), 1),
    Look = as.integer(vapply(labels, `[[`, character(1), 2)),
    "Theoretical incremental power" = vapply(
      theoretical[, "incremental"], format_result_value, character(1),
      digits = digits
    ),
    check.names = FALSE
  )

  if (!is.null(simulation)) {
    result[["Simulation incremental power"]] <- vapply(
      simulation[, "incremental"], format_result_value, character(1),
      digits = digits
    )
  }
  result[["Theoretical cumulative power"]] <- vapply(
    theoretical[, "cumulative"], format_result_value, character(1),
    digits = digits
  )

  if (!is.null(simulation)) {
    result[["Simulation cumulative power"]] <- vapply(
      simulation[, "cumulative"], format_result_value, character(1),
      digits = digits
    )
  }

  result
}

joint_power_results_table <- function(power, digits = 4) {
  formatted_power <- apply(
    power, c(1, 2), format_result_value, digits = digits
  )

  result <- data.frame(
    row_label = rownames(power),
    formatted_power,
    check.names = FALSE
  )
  names(result)[1] <- ""
  names(result)[-1] <- colnames(power)
  result
}

trial_design_table <- function(state) {
  design <- state$design
  digits <- state$options$display_digits

  data.frame(
    Item = c(
      "Treatment sample size",
      "Control sample size",
      "Allocation ratio",
      "Control median PFS",
      "Control median OS",
      "Assumed PFS HR",
      "Assumed OS HR",
      "Derived treatment median PFS",
      "Derived treatment median OS"
    ),
    Value = c(
      as.character(design$n_T),
      as.character(design$n_C),
      format_result_value(design$r, digits),
      format_result_value(log(2) / design$Lambda_C, digits),
      format_result_value(log(2) / design$lambda_2_C, digits),
      format_result_value(design$HR_PFS, digits),
      format_result_value(design$HR_OS, digits),
      format_result_value(log(2) / design$Lambda_T, digits),
      format_result_value(log(2) / design$lambda_2_T, digits)
    ),
    check.names = FALSE
  )
}

accrual_design_table <- function(state) {
  design <- state$design
  digits <- state$options$display_digits
  segment_count <- length(design$v_vec)
  start_time <- c(0, head(design$t_vec, -1))

  data.frame(
    Segment = seq_len(segment_count),
    "Start time" = vapply(start_time, format_result_value, character(1),
                           digits = digits),
    Duration = vapply(design$R_vec, format_result_value, character(1),
                      digits = digits),
    "Accrual rate" = vapply(design$v_vec, format_result_value, character(1),
                             digits = digits),
    check.names = FALSE
  )
}

analysis_schedule_table <- function(state) {
  design <- state$design
  digits <- state$options$display_digits

  data.frame(
    Look = seq_len(design$L),
    "Target PFS events" = as.integer(design$d_PFS_vec),
    "Derived OS events" = vapply(
      design$d_OS_vec, format_result_value, character(1), digits = digits
    ),
    "Calendar cutoff" = vapply(
      design$A_vec, format_result_value, character(1), digits = digits
    ),
    check.names = FALSE
  )
}

endpoint_testing_table <- function(state) {
  design <- state$design
  digits <- state$options$display_digits
  endpoints <- c("PFS", "OS")

  do.call(rbind, lapply(endpoints, function(endpoint) {
    looks <- seq_len(design$L)
    futility_looks <- design$futility_looks[[endpoint]]
    futility_positions <- match(looks, futility_looks)
    boundary <- design$boundary[paste0(endpoint, "_", looks), , drop = FALSE]

    data.frame(
      Endpoint = endpoint,
      Look = looks,
      "Futility HR" = vapply(futility_positions, function(position) {
        if (is.na(position)) {
          return("-")
        }
        format_result_value(design$futility_HR[[endpoint]][position], digits)
      }, character(1)),
      "Futility boundary" = vapply(
        boundary[, "futility"], format_boundary_value,
        character(1), digits = digits
      ),
      "Efficacy boundary" = vapply(
        boundary[, "efficacy"], format_boundary_value,
        character(1), digits = digits
      ),
      check.names = FALSE
    )
  }))
}

testing_summary_table <- function(state) {
  design <- state$design
  digits <- state$options$display_digits

  data.frame(
    Item = c("One-sided alpha", "Gatekeeping order"),
    Value = c(
      format_result_value(design$alpha, digits),
      paste(design$hierarchy_order, collapse = " -> ")
    ),
    check.names = FALSE
  )
}

write_csv_section <- function(connection, title, table) {
  writeLines(title, connection)
  utils::write.table(
    table,
    file = connection,
    sep = ",",
    row.names = FALSE,
    col.names = TRUE,
    quote = TRUE
  )
  writeLines("", connection)
}

parse_single_number <- function(value, name) {
  value <- trimws(value)
  number <- suppressWarnings(as.numeric(value))

  if (!nzchar(value) || length(number) != 1 || is.na(number) ||
      !is.finite(number)) {
    stop(name, " must be a finite number.")
  }

  number
}

optional_number <- function(value, name) {
  if (is.null(value) || length(value) != 1 || is.na(value)) {
    return(NULL)
  }
  if (!is.numeric(value) || !is.finite(value)) {
    stop(name, " must be a finite number.")
  }

  value
}

look_input_id <- function(type, endpoint, look) {
  paste(type, endpoint, look, sep = "_")
}

accrual_input_id <- function(type, segment) {
  paste(type, segment, sep = "_")
}

parse_accrual_segment_count <- function(value) {
  if (is.null(value) || length(value) != 1 || is.na(value) ||
      !is.finite(value) || value != as.integer(value) || value < 1) {
    stop("Number of accrual segments must be a positive integer.")
  }

  as.integer(value)
}

collect_accrual_configuration <- function(input) {
  segment_count <- parse_accrual_segment_count(input$accrual_segments)
  accrual_rate <- vapply(seq_len(segment_count), function(segment) {
    value <- input[[accrual_input_id("accrual_rate", segment)]]
    if (is.null(value) || is.na(value)) {
      stop("Each accrual segment must have an accrual rate.")
    }
    value
  }, numeric(1))
  start_time <- 0

  if (segment_count > 1) {
    start_time <- c(start_time, vapply(
      2:segment_count,
      function(segment) {
        value <- input[[accrual_input_id("accrual_start", segment)]]
        if (is.null(value) || is.na(value)) {
          stop("Each accrual segment after the first must have a start time.")
        }
        value
      },
      numeric(1)
    ))
  }

  list(
    t_vec = if (segment_count == 1) NULL else start_time[-1],
    v_vec = accrual_rate
  )
}

accrual_table <- function(segment_count) {
  rows <- list()

  for (segment in seq_len(segment_count)) {
    start_time_cell <- if (segment == 1) {
      tags$span("0")
    } else {
      numericInput(
        inputId = accrual_input_id("accrual_start", segment),
        label = NULL,
        value = NULL,
        min = 0,
        width = "100%"
      )
    }
    rows[[segment]] <- tags$tr(
      tags$td(segment),
      tags$td(start_time_cell),
      tags$td(numericInput(
        inputId = accrual_input_id("accrual_rate", segment),
        label = NULL,
        value = NULL,
        min = 0,
        width = "100%"
      ))
    )
  }

  tags$table(
    class = "table table-condensed table-bordered",
    style = "table-layout: fixed; width: 100%;",
    tags$thead(tags$tr(
      tags$th("Segment"),
      tags$th("Start time"),
      tags$th("Accrual rate")
    )),
    do.call(tags$tbody, rows)
  )
}

selected_looks <- function(input, type, endpoint, looks) {
  looks[vapply(looks, function(look) {
    isTRUE(input[[look_input_id(type, endpoint, look)]])
  }, logical(1))]
}

collect_look_configuration <- function(input, L) {
  endpoints <- c("PFS", "OS")
  efficacy_looks <- list()
  futility_looks <- list()
  futility_HR <- list()

  for (endpoint in endpoints) {
    efficacy_looks[[endpoint]] <- as.integer(unique(c(
      selected_looks(input, "efficacy", endpoint, seq_len(L)),
      L
    )))
    endpoint_futility_looks <- selected_looks(
      input,
      "futility",
      endpoint,
      seq_len(max(L - 1, 0))
    )
    futility_looks[[endpoint]] <- as.integer(endpoint_futility_looks)
    futility_HR[[endpoint]] <- vapply(endpoint_futility_looks, function(look) {
      futility_hr <- input[[look_input_id("futility_HR", endpoint, look)]]
      if (is.null(futility_hr) || is.na(futility_hr)) {
        stop("Each selected futility look must have a futility HR.")
      }
      futility_hr
    }, numeric(1))
  }

  list(
    efficacy_looks = efficacy_looks,
    futility_looks = futility_looks,
    futility_HR = futility_HR
  )
}

look_checkbox <- function(id, checked = FALSE, disabled = FALSE) {
  attributes <- list(type = "checkbox", id = id, value = "true")
  if (checked) {
    attributes$checked <- "checked"
  }
  if (disabled) {
    attributes$disabled <- "disabled"
  }
  do.call(tags$input, attributes)
}

look_selection_table <- function(L) {
  endpoints <- c("PFS", "OS")
  rows <- list()

  for (endpoint in endpoints) {
    for (look in seq_len(L)) {
      is_final_look <- look == L
      futility_cell <- if (is_final_look) {
        tags$span("-")
      } else {
        look_checkbox(look_input_id("futility", endpoint, look))
      }
      futility_hr_cell <- if (is_final_look) {
        tags$span("-")
      } else {
        numericInput(
          inputId = look_input_id("futility_HR", endpoint, look),
          label = NULL,
          value = 1.2,
          min = 1,
          step = 0.01,
          width = "100%"
        )
      }
      rows[[length(rows) + 1]] <- tags$tr(
        tags$td(endpoint),
        tags$td(look),
        tags$td(look_checkbox(
          look_input_id("efficacy", endpoint, look),
          checked = is_final_look,
          disabled = is_final_look
        )),
        tags$td(futility_cell),
        tags$td(futility_hr_cell)
      )
    }
  }

  tags$table(
    class = "table table-condensed table-bordered",
    tags$thead(tags$tr(
      tags$th("Endpoint"),
      tags$th("Look"),
      tags$th("Efficacy"),
      tags$th("Futility"),
      tags$th("Futility HR")
    )),
    do.call(tags$tbody, rows)
  )
}

collect_pipeline_arguments <- function(input) {
  d_PFS_vec <- parse_numeric_vector(input$d_PFS_vec, "d_PFS_vec")
  look_configuration <- collect_look_configuration(input, length(d_PFS_vec))
  accrual_configuration <- collect_accrual_configuration(input)
  treatment_arguments <- switch(input$effect_mode,
    hr = list(
      HR_PFS = input$HR_PFS,
      HR_OS = input$HR_OS
    ),
    median = list(
      median_PFS_T = input$median_PFS_T,
      median_OS_T = input$median_OS_T
    ),
    stop("Select a treatment effect input mode.")
  )

  hierarchy_arguments <- if (identical(input$hierarchy_primary, "PFS")) {
    list(hierarchy_order = c(primary = "PFS", secondary = "OS"))
  } else if (identical(input$hierarchy_primary, "OS")) {
    list(hierarchy_order = c(primary = "OS", secondary = "PFS"))
  } else {
    list()
  }

  optional_arguments <- list()
  alpha <- optional_number(input$alpha, "alpha")
  tol <- optional_number(input$tol, "tol")
  integration_seed <- optional_number(
    input$integration_seed,
    "integration_seed"
  )
  simulation_seed <- optional_number(
    input$simulation_seed,
    "simulation_seed"
  )
  display_digits <- optional_number(input$display_digits, "display_digits")
  n_sim <- optional_number(input$n_sim, "n_sim")
  if (!is.null(alpha)) {
    optional_arguments$alpha <- alpha
  }
  if (!is.null(tol)) {
    optional_arguments$tol <- tol
  }
  if (!is.null(integration_seed)) {
    optional_arguments$integration_seed <- integration_seed
  }
  if (isTRUE(input$simulation)) {
    optional_arguments$simulation <- TRUE
  }
  if (!is.null(n_sim)) {
    optional_arguments$n_sim <- n_sim
  }
  if (!is.null(simulation_seed)) {
    optional_arguments$simulation_seed <- simulation_seed
  }
  if (!is.null(display_digits)) {
    optional_arguments$display_digits <- display_digits
  }

  arguments <- list(
    n_T = input$n_T,
    n_C = input$n_C,
    d_PFS_vec = d_PFS_vec,
    t_vec = accrual_configuration$t_vec,
    v_vec = accrual_configuration$v_vec,
    median_PFS_C = input$median_PFS_C,
    median_OS_C = input$median_OS_C,
    alpha_spending_PFS = input$alpha_spending_PFS,
    alpha_spending_OS = input$alpha_spending_OS,
    alpha_spending_gamma_PFS = if (
      identical(input$alpha_spending_PFS, "HSD")
    ) {
      parse_single_number(input$alpha_spending_gamma_PFS, "PFS HSD gamma")
    } else {
      NA_real_
    },
    alpha_spending_gamma_OS = if (
      identical(input$alpha_spending_OS, "HSD")
    ) {
      parse_single_number(input$alpha_spending_gamma_OS, "OS HSD gamma")
    } else {
      NA_real_
    },
    efficacy_looks = look_configuration$efficacy_looks,
    futility_looks = look_configuration$futility_looks,
    futility_HR = look_configuration$futility_HR
  )

  c(arguments, optional_arguments, hierarchy_arguments, treatment_arguments)
}

ui <- fluidPage(
  titlePanel("Group Sequential PFS and OS Design"),
  sidebarLayout(
    sidebarPanel(
      width = 4,
      tags$h3("Basic parameters"),
      numericInput(
        "n_T",
        "Treatment group sample size",
        value = NULL,
        min = 1,
        step = 1
      ),
      numericInput(
        "n_C",
        "Control group sample size",
        value = NULL,
        min = 1,
        step = 1
      ),
      numericInput(
        "median_PFS_C",
        "Control group median PFS",
        value = NULL,
        min = 0
      ),
      numericInput(
        "median_OS_C",
        "Control group median OS",
        value = NULL,
        min = 0
      ),
      radioButtons(
        "effect_mode",
        "Input type for treatment effect",
        choices = c(
          "Hazard ratios" = "hr",
          "Median survival times" = "median"
        ),
        selected = character(0),
        inline = TRUE
      ),
      conditionalPanel(
        "input.effect_mode == 'hr'",
        numericInput(
          "HR_PFS",
          "PFS hazard ratio",
          value = NULL,
          min = 0,
          max = 1
        ),
        numericInput(
          "HR_OS",
          "OS hazard ratio",
          value = NULL,
          min = 0,
          max = 1
        )
      ),
      conditionalPanel(
        "input.effect_mode == 'median'",
        numericInput(
          "median_PFS_T",
          "Treatment group median PFS",
          value = NULL,
          min = 0
        ),
        numericInput(
          "median_OS_T",
          "Treatment group median OS",
          value = NULL,
          min = 0
        )
      ),
      tags$hr(),
      tags$h3("Patient accrual"),
      numericInput(
        "accrual_segments",
        "Number of accrual segments",
        value = NULL,
        min = 1,
        step = 1
      ),
      tags$label("Accrual configuration", class = "control-label"),
      uiOutput("accrual_configuration"),
      tags$hr(),
      tags$h3("Sequential design parameters"),
      selectInput(
        "hierarchy_primary",
        "Primary endpoint",
        choices = c("PFS", "OS"),
        selected = "PFS"
      ),
      textInput(
        "d_PFS_vec",
        "Target PFS events by look"
      ),
      numericInput(
        "alpha",
        "One-sided alpha",
        value = 0.025,
        min = 0,
        max = 1
      ),
      selectInput(
        "alpha_spending_PFS",
        "PFS alpha spending function",
        choices = c("Select a spending function" = "", "OF", "Pocock", "HSD"),
        selected = ""
      ),
      selectInput(
        "alpha_spending_OS",
        "OS alpha spending function",
        choices = c("Select a spending function" = "", "OF", "Pocock", "HSD"),
        selected = ""
      ),
      conditionalPanel(
        "input.alpha_spending_PFS == 'HSD'",
        textInput(
          "alpha_spending_gamma_PFS",
          "PFS HSD gamma"
        )
      ),
      conditionalPanel(
        "input.alpha_spending_OS == 'HSD'",
        textInput(
          "alpha_spending_gamma_OS",
          "OS HSD gamma"
        )
      ),
      tags$label("Boundary configuration", class = "control-label"),
      uiOutput("look_selection"),
      tags$hr(),
      tags$h3("Advanced settings"),
      checkboxInput(
        "simulation",
        "Run simulation"
      ),
      conditionalPanel(
        "input.simulation",
        numericInput(
          "n_sim",
          "Simulation replicates",
          value = 500,
          min = 2,
          step = 1
        ),
        numericInput(
          "simulation_seed",
          "Random seed",
          value = 1,
          min = 0,
          step = 1
        )
      ),
      checkboxInput(
        "show_computation_settings",
        "Computational settings"
      ),
      conditionalPanel(
        "input.show_computation_settings",
        numericInput(
          "tol",
          "Numerical tolerance",
          value = 1e-8,
          min = 0
        ),
        numericInput(
          "integration_seed",
          "Random seed",
          value = 1,
          min = 0,
          step = 1
        ),
        numericInput(
          "display_digits",
          "Result decimal places",
          value = 4,
          min = 0,
          max = 10,
          step = 1
        )
      ),
      actionButton(
        "run",
        "Run pipeline",
        class = "btn-primary"
      )
    ),
    mainPanel(
      width = 8,
      textOutput("status"),
      downloadButton("download_tables", "Download tables (CSV)"),
      tabsetPanel(
        tabPanel(
          "Input and design",
          tags$h3("Trial design"),
          tableOutput("trial_design"),
          tags$h3("Accrual"),
          tableOutput("accrual_design"),
          tags$h3("Analysis schedule and testing"),
          tags$h4("Analysis schedule"),
          tableOutput("analysis_schedule"),
          tags$h4("Testing configuration"),
          tableOutput("testing_summary"),
          tableOutput("endpoint_testing")
        ),
        tabPanel(
          "Results",
          tags$h3("Marginal power"),
          tableOutput("marginal_power_results"),
          tags$h3("Joint power"),
          uiOutput("joint_power_results")
        )
      )
    )
  )
)

server <- function(input, output, session) {
  output$accrual_configuration <- renderUI({
    if (is.null(input$accrual_segments) ||
        length(input$accrual_segments) != 1 ||
        is.na(input$accrual_segments)) {
      return(NULL)
    }
    tryCatch(
      {
        segment_count <- parse_accrual_segment_count(input$accrual_segments)
        accrual_table(segment_count)
      },
      error = function(error) {
        tags$p(
          class = "text-danger",
          "Enter the number of accrual segments to configure accrual."
        )
      }
    )
  })

  output$look_selection <- renderUI({
    if (is.null(input$d_PFS_vec) || !nzchar(trimws(input$d_PFS_vec))) {
      return(NULL)
    }
    tryCatch(
      {
        d_PFS_vec <- parse_numeric_vector(input$d_PFS_vec, "d_PFS_vec")
        look_selection_table(length(d_PFS_vec))
      },
      error = function(error) {
        tags$p(
          class = "text-danger",
          "Enter a valid target PFS event vector to configure looks."
        )
      }
    )
  })

  calculation <- eventReactive(
    input$run,
    {
      tryCatch(
        {
          state <- withProgress(
            message = "Running pipeline",
            value = 0,
            {
              arguments <- collect_pipeline_arguments(input)
              incProgress(0.1, detail = "Validating inputs")
              result <- do.call(run_pipeline, arguments)
              incProgress(0.9, detail = "Preparing results")
              result
            }
          )
          list(state = state, error = NULL)
        },
        error = function(error) {
          list(state = NULL, error = conditionMessage(error))
        }
      )
    },
    ignoreInit = TRUE
  )

  current_state <- reactive({
    calculation_result <- calculation()
    req(is.null(calculation_result$error))
    calculation_result$state
  })

  output$status <- renderText({
    calculation_result <- calculation()

    if (!is.null(calculation_result$error)) {
      return(paste("Error:", calculation_result$error))
    }

    "Pipeline completed."
  })

  output$download_tables <- downloadHandler(
    filename = function() {
      paste0("group-sequential-design-", Sys.Date(), ".csv")
    },
    content = function(file) {
      state <- current_state()
      digits <- state$options$display_digits
      connection <- base::file(file, open = "wt", encoding = "UTF-8")
      on.exit(close(connection), add = TRUE)

      write_csv_section(connection, "Trial design", trial_design_table(state))
      write_csv_section(connection, "Accrual", accrual_design_table(state))
      write_csv_section(
        connection,
        "Analysis schedule",
        analysis_schedule_table(state)
      )
      write_csv_section(
        connection,
        "Testing summary",
        testing_summary_table(state)
      )
      write_csv_section(
        connection,
        "Testing configuration",
        endpoint_testing_table(state)
      )
      write_csv_section(
        connection,
        "Marginal power",
        marginal_power_results_table(state)
      )
      write_csv_section(
        connection,
        "Theoretical joint power",
        joint_power_results_table(
          state$theoretical_results$joint_power_matrix,
          digits = digits
        )
      )

      if (isTRUE(state$options$simulation)) {
        write_csv_section(
          connection,
          "Simulation joint power",
          joint_power_results_table(
            state$empirical_results$joint_power_matrix,
            digits = digits
          )
        )
      }
    }
  )

  output$trial_design <- renderTable({
    trial_design_table(current_state())
  }, striped = TRUE, bordered = TRUE, rownames = FALSE)

  output$accrual_design <- renderTable({
    accrual_design_table(current_state())
  }, striped = TRUE, bordered = TRUE, rownames = FALSE)

  output$analysis_schedule <- renderTable({
    analysis_schedule_table(current_state())
  }, striped = TRUE, bordered = TRUE, rownames = FALSE)

  output$testing_summary <- renderTable({
    testing_summary_table(current_state())
  }, striped = TRUE, bordered = TRUE, rownames = FALSE)

  output$endpoint_testing <- renderTable({
    endpoint_testing_table(current_state())
  }, striped = TRUE, bordered = TRUE, rownames = FALSE)

  output$marginal_power_results <- renderTable({
    marginal_power_results_table(current_state())
  }, striped = TRUE, bordered = TRUE, rownames = FALSE)

  output$joint_power_results <- renderUI({
    state <- current_state()

    if (isTRUE(state$options$simulation)) {
      fluidRow(
        column(
          width = 6,
          tags$h4("Theoretical"),
          tableOutput("theoretical_joint_power")
        ),
        column(
          width = 6,
          tags$h4("Simulation"),
          tableOutput("simulation_joint_power")
        )
      )
    } else {
      fluidRow(
        column(
          width = 12,
          tags$h4("Theoretical"),
          tableOutput("theoretical_joint_power")
        )
      )
    }
  })

  output$theoretical_joint_power <- renderTable({
    state <- current_state()

    joint_power_results_table(
      state$theoretical_results$joint_power_matrix,
      digits = state$options$display_digits
    )
  }, striped = TRUE, bordered = TRUE, rownames = FALSE)

  output$simulation_joint_power <- renderTable({
    state <- current_state()
    req(isTRUE(state$options$simulation))

    joint_power_results_table(
      state$empirical_results$joint_power_matrix,
      digits = state$options$display_digits
    )
  }, striped = TRUE, bordered = TRUE, rownames = FALSE)
}

shinyApp(ui = ui, server = server)
