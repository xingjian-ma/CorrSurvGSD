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

collect_pipeline_arguments <- function(input) {
  treatment_arguments <- if (input$effect_mode == "hr") {
    list(
      HR_PFS = input$HR_PFS,
      HR_OS = input$HR_OS
    )
  } else {
    list(
      median_PFS_T = input$median_PFS_T,
      median_OS_T = input$median_OS_T
    )
  }

  hierarchy_secondary <- if (input$hierarchy_primary == "PFS") {
    "OS"
  } else {
    "PFS"
  }

  arguments <- list(
    n_T = input$n_T,
    n_C = input$n_C,
    d_PFS_vec = parse_numeric_vector(
      input$d_PFS_vec,
      "d_PFS_vec"
    ),
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
    alpha = input$alpha,
    futility_analysis = c(
      PFS = input$futility_PFS,
      OS = input$futility_OS
    ),
    futility_HR = c(
      PFS = input$futility_HR_PFS,
      OS = input$futility_HR_OS
    ),
    efficacy_start = c(
      PFS = if (input$futility_PFS) input$efficacy_start_PFS else 1,
      OS = if (input$futility_OS) input$efficacy_start_OS else 1
    ),
    hierarchy_order = c(
      primary = input$hierarchy_primary,
      secondary = hierarchy_secondary
    ),
    tol = input$tol,
    seed = input$seed,
    simulation = input$simulation,
    n_sim = input$n_sim
  )

  c(arguments, treatment_arguments)
}

ui <- fluidPage(
  titlePanel("Group Sequential PFS and OS Design"),
  sidebarLayout(
    sidebarPanel(
      width = 4,
      numericInput(
        "n_T",
        "Treatment sample size",
        value = 800,
        min = 1,
        step = 1
      ),
      numericInput(
        "n_C",
        "Control sample size",
        value = 1200,
        min = 1,
        step = 1
      ),
      textInput(
        "d_PFS_vec",
        "Target PFS events by look",
        value = "400, 800, 1200"
      ),
      textInput(
        "t_vec",
        "Accrual breakpoints (empty for one segment)",
        value = "6, 10"
      ),
      textInput(
        "v_vec",
        "Accrual rates by segment",
        value = "8.333333, 8.333333, 8.333333"
      ),
      numericInput(
        "median_PFS_C",
        "Control median PFS",
        value = log(2) / 0.25,
        min = 0
      ),
      numericInput(
        "median_OS_C",
        "Control median OS",
        value = log(2) / 0.10,
        min = 0
      ),
      radioButtons(
        "effect_mode",
        "Treatment effect input",
        choices = c(
          "Hazard ratios" = "hr",
          "Treatment medians" = "median"
        ),
        selected = "hr",
        inline = TRUE
      ),
      conditionalPanel(
        "input.effect_mode == 'hr'",
        numericInput(
          "HR_PFS",
          "PFS hazard ratio",
          value = 0.8,
          min = 0,
          max = 1
        ),
        numericInput(
          "HR_OS",
          "OS hazard ratio",
          value = 0.8,
          min = 0,
          max = 1
        )
      ),
      conditionalPanel(
        "input.effect_mode == 'median'",
        numericInput(
          "median_PFS_T",
          "Treatment median PFS",
          value = log(2) / 0.20,
          min = 0
        ),
        numericInput(
          "median_OS_T",
          "Treatment median OS",
          value = log(2) / 0.08,
          min = 0
        )
      ),
      selectInput(
        "alpha_spending_PFS",
        "PFS alpha spending",
        choices = c("WT", "OF", "Pocock"),
        selected = "OF"
      ),
      selectInput(
        "alpha_spending_OS",
        "OS alpha spending",
        choices = c("WT", "OF", "Pocock"),
        selected = "Pocock"
      ),
      checkboxInput(
        "futility_PFS",
        "Enable PFS futility",
        value = TRUE
      ),
      numericInput(
        "futility_HR_PFS",
        "PFS futility HR",
        value = 1.2,
        min = 1
      ),
      checkboxInput(
        "futility_OS",
        "Enable OS futility",
        value = FALSE
      ),
      numericInput(
        "futility_HR_OS",
        "OS futility HR",
        value = 1.2,
        min = 1
      ),
      conditionalPanel(
        "input.futility_PFS",
        numericInput(
          "efficacy_start_PFS",
          "First PFS efficacy look",
          value = 1,
          min = 1,
          step = 1
        )
      ),
      conditionalPanel(
        "input.futility_OS",
        numericInput(
          "efficacy_start_OS",
          "First OS efficacy look",
          value = 1,
          min = 1,
          step = 1
        )
      ),
      selectInput(
        "hierarchy_primary",
        "Primary endpoint",
        choices = c("PFS", "OS"),
        selected = "PFS"
      ),
      checkboxInput(
        "show_advanced_settings",
        "Show advanced settings",
        value = FALSE
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
          value = 777,
          min = 0,
          step = 1
        )
      ),
      checkboxInput(
        "simulation",
        "Run simulation",
        value = FALSE
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
