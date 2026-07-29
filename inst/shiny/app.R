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

without_mean_vector <- function(results) {
  results$joint_mean_vector <- NULL
  results
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
          value = NULL,
          min = 1,
          step = 0.01,
          width = "90px"
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
  seed <- optional_number(input$seed, "seed")
  n_sim <- optional_number(input$n_sim, "n_sim")
  if (!is.null(alpha)) {
    optional_arguments$alpha <- alpha
  }
  if (!is.null(tol)) {
    optional_arguments$tol <- tol
  }
  if (!is.null(seed)) {
    optional_arguments$seed <- seed
  }
  if (isTRUE(input$simulation)) {
    optional_arguments$simulation <- TRUE
  }
  if (!is.null(n_sim)) {
    optional_arguments$n_sim <- n_sim
  }

  arguments <- list(
    n_T = input$n_T,
    n_C = input$n_C,
    d_PFS_vec = d_PFS_vec,
    t_vec = parse_numeric_vector(
      input$t_vec,
      "t_vec",
      allow_empty = TRUE
    ),
    v_vec = parse_numeric_vector(
      input$v_vec,
      "v_vec"
    ),
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
      numericInput(
        "n_T",
        "Treatment sample size",
        value = NULL,
        min = 1,
        step = 1
      ),
      numericInput(
        "n_C",
        "Control sample size",
        value = NULL,
        min = 1,
        step = 1
      ),
      textInput(
        "d_PFS_vec",
        "Target PFS events by look"
      ),
      textInput(
        "t_vec",
        "Accrual start times (optional, leave blank for uniform accrual)"
      ),
      textInput(
        "v_vec",
        "Accrual rates by segment"
      ),
      numericInput(
        "median_PFS_C",
        "Control median PFS",
        value = NULL,
        min = 0
      ),
      numericInput(
        "median_OS_C",
        "Control median OS",
        value = NULL,
        min = 0
      ),
      radioButtons(
        "effect_mode",
        "Treatment effect input",
        choices = c(
          "Hazard ratios" = "hr",
          "Treatment medians" = "median"
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
          "Treatment median PFS",
          value = NULL,
          min = 0
        ),
        numericInput(
          "median_OS_T",
          "Treatment median OS",
          value = NULL,
          min = 0
        )
      ),
      selectInput(
        "alpha_spending_PFS",
        "PFS alpha spending",
        choices = c("Select a spending family" = "", "OF", "Pocock", "HSD"),
        selected = ""
      ),
      selectInput(
        "alpha_spending_OS",
        "OS alpha spending",
        choices = c("Select a spending family" = "", "OF", "Pocock", "HSD"),
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
      tags$hr(),
      tags$h4("Look configuration"),
      tags$p("Select efficacy and futility analyses by endpoint and look."),
      uiOutput("look_selection"),
      selectInput(
        "hierarchy_primary",
        "Primary endpoint",
        choices = c("PFS", "OS"),
        selected = "PFS"
      ),
      checkboxInput(
        "show_advanced_settings",
        "Show advanced settings"
      ),
      conditionalPanel(
        "input.show_advanced_settings",
        numericInput(
          "alpha",
          "One-sided alpha",
          value = 0.025,
          min = 0,
          max = 1
        ),
        numericInput(
          "tol",
          "Numerical tolerance",
          value = 1e-8,
          min = 0
        ),
        numericInput(
          "seed",
          "Random seed",
          value = 1,
          min = 0,
          step = 1
        )
      ),
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
      tabsetPanel(
        tabPanel(
          "Input and design",
          verbatimTextOutput("design")
        ),
        tabPanel(
          "Theoretical results",
          verbatimTextOutput("theoretical")
        ),
        tabPanel(
          "Simulation results",
          verbatimTextOutput("simulation_results")
        )
      )
    )
  )
)

server <- function(input, output, session) {
  output$look_selection <- renderUI({
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

  output$design <- renderPrint({
    current_state()$design
  })

  output$theoretical <- renderPrint({
    without_mean_vector(current_state()$theoretical_results)
  })

  output$simulation_results <- renderPrint({
    state <- current_state()

    if (!state$options$simulation) {
      cat("Simulation was not requested.")
      return(invisible(NULL))
    }

    without_mean_vector(state$empirical_results)
  })
}

shinyApp(ui = ui, server = server)
