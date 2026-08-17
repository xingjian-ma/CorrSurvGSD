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
    return("-")
  }

  formatC(value, format = "f", digits = digits)
}

format_boundary_value <- function(value, digits = 4) {
  if (is.na(value)) {
    return("-")
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
    "Incremental power" = vapply(
      seq_len(nrow(theoretical)),
      function(index) {
        format_power_comparison(
          theoretical[index, "incremental"],
          if (is.null(simulation)) NULL else simulation[index, "incremental"],
          digits
        )
      },
      character(1)
    ),
    "Cumulative power" = vapply(
      seq_len(nrow(theoretical)),
      function(index) {
        format_power_comparison(
          theoretical[index, "cumulative"],
          if (is.null(simulation)) NULL else simulation[index, "cumulative"],
          digits
        )
      },
      character(1)
    ),
    check.names = FALSE
  )

  result
}

format_power_comparison <- function(closed_form, simulation = NULL,
                                    digits = 4) {
  closed_form_value <- format_result_value(closed_form, digits)

  if (is.null(simulation) || is.na(closed_form)) {
    return(closed_form_value)
  }

  paste0(
    closed_form_value,
    " (",
    format_result_value(simulation, digits),
    ")"
  )
}

joint_power_results_table <- function(power, simulation_power = NULL,
                                      digits = 4) {
  formatted_power <- matrix(
    NA_character_,
    nrow = nrow(power),
    ncol = ncol(power),
    dimnames = dimnames(power)
  )
  for (row in seq_len(nrow(power))) {
    for (column in seq_len(ncol(power))) {
      simulation_value <- if (is.null(simulation_power)) {
        NULL
      } else {
        simulation_power[row, column]
      }
      formatted_power[row, column] <- format_power_comparison(
        power[row, column],
        simulation_value,
        digits
      )
    }
  }
  row_labels <- gsub("_", ", look ", rownames(power), fixed = TRUE)
  column_labels <- gsub("_", ", look ", colnames(power), fixed = TRUE)

  result <- data.frame(
    row_label = row_labels,
    formatted_power,
    check.names = FALSE
  )
  names(result)[1] <- "Endpoint, look"
  names(result)[-1] <- column_labels
  result
}

trial_design_table <- function(state) {
  design <- state$design
  digits <- state$options$display_digits

  data.frame(
    Characteristic = c(
      "Sample size",
      "Median PFS",
      "Median OS"
    ),
    Treatment = c(
      as.character(design$n_T),
      format_result_value(log(2) / design$Lambda_T, digits),
      format_result_value(log(2) / design$lambda_2_T, digits)
    ),
    Control = c(
      as.character(design$n_C),
      format_result_value(log(2) / design$Lambda_C, digits),
      format_result_value(log(2) / design$lambda_2_C, digits)
    ),
    "Comparative measure" = c(
      paste0("Allocation ratio (T/C): ", format_result_value(design$r, digits)),
      paste0("PFS HR (T/C): ", format_result_value(design$HR_PFS, digits)),
      paste0("OS HR (T/C): ", format_result_value(design$HR_OS, digits))
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
    boundary_HR <- design$boundary_HR[
      paste0(endpoint, "_", looks),
      ,
      drop = FALSE
    ]
    boundary_p <- design$boundary_p[
      paste0(endpoint, "_", looks),
      ,
      drop = FALSE
    ]

    data.frame(
      Endpoint = endpoint,
      Look = looks,
      "Futility HR threshold" = vapply(futility_positions, function(position) {
        if (is.na(position)) {
          return("-")
        }
        format_result_value(design$futility_HR[[endpoint]][position], digits)
      }, character(1)),
      "Futility boundary (Z)" = vapply(
        boundary[, "futility"], format_boundary_value,
        character(1), digits = digits
      ),
      "Efficacy boundary (Z)" = vapply(
        boundary[, "efficacy"], format_boundary_value,
        character(1), digits = digits
      ),
      "Futility boundary (HR)" = vapply(
        boundary_HR[, "futility"], format_boundary_value,
        character(1), digits = digits
      ),
      "Efficacy boundary (HR)" = vapply(
        boundary_HR[, "efficacy"], format_boundary_value,
        character(1), digits = digits
      ),
      "Futility boundary (p)" = vapply(
        boundary_p[, "futility"], format_boundary_value,
        character(1), digits = digits
      ),
      "Efficacy boundary (p)" = vapply(
        boundary_p[, "efficacy"], format_boundary_value,
        character(1), digits = digits
      ),
      check.names = FALSE
    )
  }))
}

boundary_results_table <- function(state, scale = c("z", "hr", "p")) {
  scale <- match.arg(scale)
  design <- state$design
  digits <- state$options$display_digits
  boundary <- switch(
    scale,
    z = design$boundary,
    hr = design$boundary_HR,
    p = design$boundary_p
  )

  do.call(rbind, lapply(c("PFS", "OS"), function(endpoint) {
    looks <- seq_len(design$L)
    futility_looks <- design$futility_looks[[endpoint]]
    futility_positions <- match(looks, futility_looks)
    endpoint_boundary <- boundary[
      paste0(endpoint, "_", looks),
      ,
      drop = FALSE
    ]

    data.frame(
      Endpoint = endpoint,
      Look = looks,
      "Futility HR" = vapply(
        futility_positions,
        function(position) {
          if (is.na(position)) {
            return("-")
          }
          format_result_value(design$futility_HR[[endpoint]][position], digits)
        },
        character(1)
      ),
      "Futility boundary" = vapply(
        endpoint_boundary[, "futility"],
        format_boundary_value,
        character(1),
        digits = digits
      ),
      "Efficacy boundary" = vapply(
        endpoint_boundary[, "efficacy"],
        format_boundary_value,
        character(1),
        digits = digits
      ),
      check.names = FALSE
    )
  }))
}

alpha_spending_label <- function(design, endpoint, digits) {
  spending <- design[[paste0("alpha_spending_", endpoint)]]
  gamma <- design[[paste0("alpha_spending_gamma_", endpoint)]]

  if (identical(spending, "HSD")) {
    return(paste0("HSD (gamma = ", format_result_value(gamma, digits), ")"))
  }

  spending
}

testing_framework_details <- function(state) {
  design <- state$design
  digits <- state$options$display_digits

  data.frame(
    Label = c(
      "One-sided alpha",
      "Gatekeeping order",
      "PFS alpha spending",
      "OS alpha spending"
    ),
    Value = c(
      format_result_value(design$alpha, digits),
      paste(design$hierarchy_order, collapse = " -> "),
      alpha_spending_label(design, "PFS", digits),
      alpha_spending_label(design, "OS", digits)
    ),
    check.names = FALSE
  )
}

result_context_text <- function(state) {
  design <- state$design

  paste(
    paste0(design$n_T + design$n_C, " participants"),
    paste0(design$L, " analysis looks"),
    paste(design$hierarchy_order, collapse = " -> "),
    sep = " · "
  )
}

testing_summary_table <- function(state) {
  design <- state$design
  digits <- state$options$display_digits

  data.frame(
    Item = c("One-sided Type I error", "Gatekeeping order"),
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
      tags$th("Start time (calendar time)"),
      tags$th("Accrual rate (participants per time unit)")
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
      tags$th("Efficacy analysis"),
      tags$th("Futility analysis"),
      tags$th("Futility HR threshold")
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

is_missing_input <- function(value) {
  is.null(value) || length(value) != 1 || is.na(value) ||
    (is.character(value) && !nzchar(trimws(value)))
}

required_input_labels <- function(input) {
  missing <- list()
  add_missing <- function(id, label) {
    missing[[id]] <<- label
  }

  if (is_missing_input(input$n_T)) {
    add_missing("n_T", "Treatment group sample size")
  }
  if (is_missing_input(input$n_C)) {
    add_missing("n_C", "Control group sample size")
  }
  if (is_missing_input(input$median_PFS_C)) {
    add_missing("median_PFS_C", "Control group median PFS")
  }
  if (is_missing_input(input$median_OS_C)) {
    add_missing("median_OS_C", "Control group median OS")
  }
  if (is_missing_input(input$effect_mode)) {
    add_missing("effect_mode", "Treatment effect input type")
  }
  if (identical(input$effect_mode, "hr")) {
    if (is_missing_input(input$HR_PFS)) {
      add_missing("HR_PFS", "PFS hazard ratio")
    }
    if (is_missing_input(input$HR_OS)) {
      add_missing("HR_OS", "OS hazard ratio")
    }
  }
  if (identical(input$effect_mode, "median")) {
    if (is_missing_input(input$median_PFS_T)) {
      add_missing("median_PFS_T", "Treatment group median PFS")
    }
    if (is_missing_input(input$median_OS_T)) {
      add_missing("median_OS_T", "Treatment group median OS")
    }
  }
  if (is_missing_input(input$accrual_segments)) {
    add_missing("accrual_segments", "Number of accrual segments")
  } else {
    segment_count <- suppressWarnings(as.integer(input$accrual_segments))
    has_segment_table <- !is.na(segment_count) && segment_count > 0 &&
      input$accrual_segments == segment_count
    if (has_segment_table) {
      for (segment in seq_len(segment_count)) {
        rate_id <- accrual_input_id("accrual_rate", segment)
        if (is_missing_input(input[[rate_id]])) {
          add_missing(rate_id, paste("Accrual rate for segment", segment))
        }
        if (segment > 1) {
          start_id <- accrual_input_id("accrual_start", segment)
          if (is_missing_input(input[[start_id]])) {
            add_missing(start_id, paste("Start time for segment", segment))
          }
        }
      }
    }
  }
  if (is_missing_input(input$d_PFS_vec)) {
    add_missing("d_PFS_vec", "Target PFS events by look")
  }
  if (is_missing_input(input$alpha_spending_PFS)) {
    add_missing("alpha_spending_PFS", "PFS alpha-spending function")
  }
  if (is_missing_input(input$alpha_spending_OS)) {
    add_missing("alpha_spending_OS", "OS alpha-spending function")
  }
  if (identical(input$alpha_spending_PFS, "HSD") &&
      is_missing_input(input$alpha_spending_gamma_PFS)) {
    add_missing("alpha_spending_gamma_PFS", "PFS HSD gamma")
  }
  if (identical(input$alpha_spending_OS, "HSD") &&
      is_missing_input(input$alpha_spending_gamma_OS)) {
    add_missing("alpha_spending_gamma_OS", "OS HSD gamma")
  }

  missing
}

translate_pipeline_error <- function(message) {
  mappings <- list(
    list("d_PFS_vec must be a comma-separated", paste(
      "Enter target PFS events as comma-separated numbers, for example: ",
      "50, 100."
    )),
    list("Select a treatment effect input mode", paste(
      "Choose either Hazard ratios or Median survival times for the ",
      "treatment effect."
    )),
    list("Number of accrual segments", paste(
      "Enter the number of accrual segments as a positive whole number."
    )),
    list("PFS HSD gamma must be a finite number", paste(
      "Enter PFS HSD gamma as a finite number."
    )),
    list("OS HSD gamma must be a finite number", paste(
      "Enter OS HSD gamma as a finite number."
    )),
    list("Each accrual segment must have an accrual rate", paste(
      "Enter an accrual rate for every configured accrual segment."
    )),
    list("Each accrual segment after the first must have a start time", paste(
      "Enter a start time for every accrual segment after the first."
    )),
    list("Each selected futility look must have a futility HR", paste(
      "Enter a futility HR for every selected futility look."
    )),
    list("median_PFS_C must be less than median_OS_C", paste(
      "Control median PFS must be shorter than control median OS."
    )),
    list("median_PFS_T must be less than median_OS_T", paste(
      "Treatment median PFS must be shorter than treatment median OS."
    )),
    list("lambda_1_T must be positive", paste(
      "The selected PFS and OS treatment effects are incompatible with ",
      "the illness-to-death model. Adjust the treatment assumptions."
    )),
    list("R_s must be positive", paste(
      "The accrual schedule cannot enroll the planned sample size. Adjust ",
      "the segment start times or accrual rates."
    )),
    list("t_vec must be strictly increasing", paste(
      "Accrual segment start times must be strictly increasing."
    )),
    list("t_vec must have length", paste(
      "Configure one start time for each accrual segment after the first."
    )),
    list("d_PFS_ell must not exceed", paste(
      "Target PFS events by look must not exceed the total sample ",
      "size."
    )),
    list("d_PFS_vec", paste(
      "Target PFS events by look must be strictly increasing positive whole ",
      "numbers."
    )),
    list("efficacy_looks", paste(
      "Each endpoint must include the final look in its efficacy schedule."
    )),
    list("futility_looks", paste(
      "Futility looks must be selected before the final analysis only."
    )),
    list("futility_HR", paste(
      "Each futility hazard ratio must be a finite number greater than 1."
    )),
    list("alpha_spending_gamma", paste(
      "Each HSD gamma must be a finite number from -40 inclusive to 40 ",
      "exclusive."
    )),
    list("alpha_spending", paste(
      "Select OF, Pocock, or HSD for each alpha-spending function."
    )),
    list("n_sim must be at least 2", paste(
      "Enter at least 2 simulation replicates when simulation is enabled."
    )),
    list("display_digits", paste(
      "Displayed decimal places must be a whole number from 0 to 10."
    )),
    list("Assertion on 'n_T'", paste(
      "Enter treatment group sample size as a positive whole number."
    )),
    list("Assertion on 'n_C'", paste(
      "Enter control group sample size as a positive whole number."
    )),
    list("Assertion on 'median_PFS_C'", paste(
      "Enter control median PFS as a positive finite number."
    )),
    list("Assertion on 'median_OS_C'", paste(
      "Enter control median OS as a positive finite number."
    )),
    list("Assertion on 'median_PFS_T'", paste(
      "Enter treatment median PFS as a positive finite number."
    )),
    list("Assertion on 'median_OS_T'", paste(
      "Enter treatment median OS as a positive finite number."
    )),
    list("Assertion on 'HR_PFS'", paste(
      "Enter PFS hazard ratio as a number greater than 0 and no greater ",
      "than 1."
    )),
    list("Assertion on 'HR_OS'", paste(
      "Enter OS hazard ratio as a number greater than 0 and no greater ",
      "than 1."
    )),
    list("Assertion on 'v_vec'", paste(
      "Enter a positive accrual rate for every accrual segment."
    )),
    list("Assertion on 'alpha'", paste(
      "Enter one-sided alpha as a finite number strictly between 0 and 1."
    )),
    list("alpha must be strictly between", paste(
      "Enter one-sided alpha as a finite number strictly between 0 and 1."
    )),
    list("Assertion on 'tol'", paste(
      "Enter numerical tolerance as a positive finite number."
    )),
    list("Assertion on 'integration_seed'", paste(
      "Enter the integration random seed as a non-negative whole number."
    )),
    list("Assertion on 'simulation_seed'", paste(
      "Enter the simulation random seed as a non-negative whole number."
    )),
    list("Assertion on 'n_sim'", paste(
      "Enter simulation replicates as a whole number of at least 2."
    ))
  )

  for (mapping in mappings) {
    if (grepl(mapping[[1]], message, fixed = TRUE)) {
      return(gsub(" {2,}", " ", mapping[[2]]))
    }
  }

  paste0(
    "The trial design could not be evaluated. Check the entered inputs ",
    "and try again."
  )
}

pipeline_error_input_ids <- function(message) {
  patterns <- list(
    "d_PFS_vec" = "d_PFS_vec",
    "median_PFS_C" = c("median_PFS_C", "median_OS_C"),
    "median_OS_C" = c("median_PFS_C", "median_OS_C"),
    "median_PFS_T" = c("median_PFS_T", "median_OS_T"),
    "median_OS_T" = c("median_PFS_T", "median_OS_T"),
    "lambda_1_T" = c("HR_PFS", "HR_OS", "median_PFS_T", "median_OS_T"),
    "HR_PFS" = "HR_PFS",
    "HR_OS" = "HR_OS",
    "n_T" = "n_T",
    "n_C" = "n_C",
    "v_vec" = "accrual_segments",
    "t_vec" = "accrual_segments",
    "R_s" = "accrual_segments",
    "alpha_spending_gamma_PFS" = "alpha_spending_gamma_PFS",
    "alpha_spending_gamma_OS" = "alpha_spending_gamma_OS",
    "alpha_spending" = c("alpha_spending_PFS", "alpha_spending_OS"),
    "efficacy_looks" = "look_selection",
    "futility_looks" = "look_selection",
    "futility_HR" = "look_selection",
    "alpha" = "alpha",
    "tol" = "tol",
    "integration_seed" = "integration_seed",
    "simulation_seed" = "simulation_seed",
    "n_sim" = "n_sim",
    "display_digits" = "display_digits"
  )

  for (pattern in names(patterns)) {
    if (grepl(pattern, message, fixed = TRUE)) {
      return(patterns[[pattern]])
    }
  }

  character(0)
}

ui <- fluidPage(
  tags$head(
    tags$style(HTML(
      "
      body {
        background-color: #f7f5f5;
        color: #31292b;
        font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
      }
      .app-title {
        margin: 0 -15px 20px;
        padding: 24px 30px 20px;
        background: linear-gradient(135deg, #7d1020, #b52239);
        color: #ffffff;
        box-shadow: 0 3px 12px rgba(92, 15, 28, 0.2);
      }
      .app-title h2 {
        margin: 0;
        font-weight: 600;
      }
      .well {
        border: 1px solid #e4dadd;
        border-radius: 10px;
        background-color: #ffffff;
        box-shadow: 0 2px 8px rgba(49, 41, 43, 0.06);
      }
      .sidebar-section-title {
        margin-top: 20px;
        padding-bottom: 8px;
        border-bottom: 2px solid #f0d8dc;
        color: #7d1020;
        font-size: 17px;
        font-weight: 600;
      }
      .required-note {
        float: right;
        color: #5b7180;
        font-size: 12px;
        font-weight: 400;
      }
      .app-run-button {
        width: 100%;
        margin-top: 12px;
        border: 0;
        border-radius: 6px;
        background-color: #9f1d32;
        font-weight: 600;
      }
      .app-run-button:hover,
      .app-run-button:focus {
        background-color: #7d1020;
      }
      .app-download-button {
        margin-bottom: 14px;
      }
      .app-status {
        min-height: 36px;
        margin-bottom: 14px;
        padding: 9px 13px;
        border-left: 4px solid #9f1d32;
        border-radius: 4px;
        background-color: #fdf1f3;
        color: #7d1020;
        font-weight: 600;
      }
      .app-status-error {
        border-left-color: #b44747;
        background-color: #fff4f4;
        color: #923535;
      }
      .app-status-success {
        border-left-color: #2c7a56;
        background-color: #eff9f3;
        color: #236345;
      }
      .required-marker {
        margin-left: 4px;
        color: #b44747;
      }
      .input-help {
        margin-top: -10px;
        margin-bottom: 10px;
        color: #5b7180;
        font-size: 12px;
      }
      .configuration-label {
        display: block;
        margin-top: 8px;
        margin-bottom: 8px;
      }
      .configuration-help {
        margin-top: 0;
      }
      .app-help {
        margin: 4px 0 12px;
        color: #406070;
        font-size: 13px;
      }
      .required-inputs {
        margin: 14px 0 4px;
        padding: 10px 12px;
        border: 1px solid #ead0a7;
        border-radius: 6px;
        background-color: #fffbf2;
        color: #76511d;
        font-size: 13px;
      }
      .required-inputs p {
        margin: 0 0 5px;
        font-weight: 600;
      }
      .required-inputs ul {
        margin-bottom: 0;
        padding-left: 18px;
      }
      .required-input-link {
        color: #76511d;
        text-decoration: underline;
      }
      .input-complete {
        margin: 14px 0 4px;
        padding: 10px 12px;
        border: 1px solid #b7dbc6;
        border-radius: 6px;
        background-color: #f2fbf5;
        color: #236345;
        font-size: 13px;
      }
      .app-input-error input,
      .app-input-error select {
        border-color: #b44747;
        box-shadow: 0 0 0 2px rgba(180, 71, 71, 0.14);
      }
      td.app-input-error {
        background-color: #fff4f4;
      }
      .empty-state {
        margin: 24px 0;
        padding: 28px;
        border: 1px dashed #b9cbd6;
        border-radius: 8px;
        background-color: #f8fbfc;
        color: #5b7180;
        text-align: center;
      }
      .empty-state-error {
        border-color: #d9a7a7;
        background-color: #fff8f8;
        color: #9b3d3d;
      }
      .intro-hero {
        margin-bottom: 20px;
        padding: 30px 32px;
        border-radius: 10px;
        background: linear-gradient(135deg, #fff4f5, #ffffff);
        border: 1px solid #edcfd4;
      }
      .intro-hero h2 {
        margin-top: 0;
        color: #7d1020;
        font-weight: 600;
      }
      .intro-grid {
        display: grid;
        grid-template-columns: repeat(auto-fit, minmax(210px, 1fr));
        gap: 16px;
        margin-bottom: 20px;
      }
      .intro-card {
        min-height: 162px;
        padding: 18px;
        border: 1px solid #e4dadd;
        border-radius: 8px;
        background-color: #ffffff;
        box-shadow: 0 2px 6px rgba(49, 41, 43, 0.05);
      }
      .intro-card h4 {
        margin-top: 0;
        color: #9f1d32;
        font-weight: 600;
      }
      .intro-step {
        display: inline-block;
        width: 26px;
        height: 26px;
        margin-right: 6px;
        border-radius: 50%;
        background-color: #9f1d32;
        color: #ffffff;
        line-height: 26px;
        text-align: center;
        font-weight: 600;
      }
      .tab-content {
        padding-top: 18px;
      }
      .intro-eyebrow {
        margin-bottom: 8px;
        color: #9f1d32;
        font-size: 12px;
        font-weight: 700;
        letter-spacing: 0.08em;
        text-transform: uppercase;
      }
      .intro-lead {
        max-width: 760px;
        margin-bottom: 0;
        color: #53464a;
        font-size: 16px;
        line-height: 1.6;
      }
      .intro-output-list {
        display: grid;
        grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
        gap: 12px;
        margin: 0;
        padding: 0;
        list-style: none;
      }
      .intro-output-list li {
        padding: 14px 16px;
        border-left: 3px solid #b52239;
        border-radius: 4px;
        background-color: #ffffff;
        box-shadow: 0 1px 4px rgba(49, 41, 43, 0.05);
      }
      .result-context {
        margin-bottom: 18px;
        padding: 13px 16px;
        border: 1px solid #eadfe1;
        border-radius: 6px;
        background-color: #fcfaf9;
        color: #68565a;
      }
      .result-context-title {
        margin-bottom: 5px;
        color: #7d1020;
        font-size: 18px;
        font-weight: 700;
      }
      .result-context-details {
        font-size: 13px;
      }
      .result-section-label {
        display: block;
        margin: 18px 0 8px;
        color: #7d1020;
        font-size: 14px;
        font-weight: 700;
        letter-spacing: 0.04em;
        text-transform: uppercase;
      }
      .result-sections .tab-pane > .result-section-label:first-child {
        margin-top: 0;
      }
      .result-sections > .nav-tabs > li.active > a,
      .result-sections > .nav-tabs > li.active > a:focus,
      .result-sections > .nav-tabs > li.active > a:hover {
        border-top: 2px solid #9f1d32;
        color: #7d1020;
      }
      @media (max-width: 767px) {
        .app-title {
          margin-bottom: 14px;
          padding: 20px;
        }
        .intro-hero {
          padding: 22px;
        }
      }
      table {
        background-color: #ffffff;
      }
      "
    )),
    tags$script(HTML(
      "Shiny.addCustomMessageHandler('corrsurvInputError', function(ids) {
        $('.app-input-error').removeClass('app-input-error');
        (Array.isArray(ids) ? ids : []).forEach(function(id) {
          var element = $('#' + id);
          element.closest('.form-group').addClass('app-input-error');
          element.closest('td').addClass('app-input-error');
          element.addClass('app-input-error');
        });
      });"
    ))
  ),
  div(class = "app-title", titlePanel("Group-Sequential PFS and OS Design")),
  sidebarLayout(
    sidebarPanel(
      width = 4,
      tags$h3(
        tags$span("Basic parameters"),
        tags$span("* Required", class = "required-note"),
        class = "sidebar-section-title"
      ),
      numericInput(
        "n_T",
        tagList("Treatment group sample size", tags$span(
          "*", class = "required-marker"
        )),
        value = NULL,
        min = 1,
        step = 1
      ),
      tags$p("Number of participants assigned to treatment.",
             class = "input-help"),
      numericInput(
        "n_C",
        tagList("Control group sample size", tags$span(
          "*", class = "required-marker"
        )),
        value = NULL,
        min = 1,
        step = 1
      ),
      tags$p("Number of participants assigned to control.",
             class = "input-help"),
      numericInput(
        "median_PFS_C",
        tagList("Control group median PFS", tags$span(
          "*", class = "required-marker"
        )),
        value = NULL,
        min = 0
      ),
      tags$p(
        "Median PFS expected in the control group. Use the same time unit ",
        "for all survival medians.",
        class = "input-help"
      ),
      numericInput(
        "median_OS_C",
        tagList("Control group median OS", tags$span(
          "*", class = "required-marker"
        )),
        value = NULL,
        min = 0
      ),
      tags$p(
        "Median OS expected in the control group. Use the same time unit ",
        "for all survival medians.",
        class = "input-help"
      ),
      radioButtons(
        "effect_mode",
        tagList("Input type for treatment effect", tags$span(
          "*", class = "required-marker"
        )),
        choices = c(
          "Hazard ratios" = "hr",
          "Median survival times" = "median"
        ),
        selected = character(0),
        inline = TRUE
      ),
      tags$p(
        "Choose Hazard ratios to specify treatment-to-control effects, or ",
        "Median survival times to specify treatment-group medians directly.",
        class = "input-help"
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
        tags$p(
          "Treatment-to-control PFS hazard ratio. Values below 1 indicate ",
          "a lower progression or death hazard with treatment.",
          class = "input-help"
        ),
        numericInput(
          "HR_OS",
          "OS hazard ratio",
          value = NULL,
          min = 0,
          max = 1
        ),
        tags$p(
          "Treatment-to-control OS hazard ratio. Values below 1 indicate ",
          "a lower death hazard with treatment.",
          class = "input-help"
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
        tags$p(
          "Median PFS expected in the treatment group. Use the same time ",
          "unit for all survival medians.",
          class = "input-help"
        ),
        numericInput(
          "median_OS_T",
          "Treatment group median OS",
          value = NULL,
          min = 0
        ),
        tags$p(
          "Median OS expected in the treatment group. Use the same time ",
          "unit for all survival medians.",
          class = "input-help"
        )
      ),
      tags$hr(),
      tags$h3("Patient accrual settings", class = "sidebar-section-title"),
      numericInput(
        "accrual_segments",
        tagList("Number of accrual segments", tags$span(
          "*", class = "required-marker"
        )),
        value = NULL,
        min = 1,
        step = 1
      ),
      tags$p(
        "Use one segment for constant accrual, or multiple segments when ",
        "the accrual rate changes over calendar time.",
        class = "input-help"
      ),
      tags$label(
        "Accrual configuration",
        class = "control-label configuration-label"
      ),
      tags$p(
        "Start time is 0 for the first segment. Later start times must use ",
        "the same calendar-time unit as the survival medians.",
        class = "input-help configuration-help"
      ),
      uiOutput("accrual_configuration"),
      tags$hr(),
      tags$h3(
        "Sequential design parameters",
        class = "sidebar-section-title"
      ),
      selectInput(
        "hierarchy_primary",
        "Primary endpoint",
        choices = c("PFS", "OS"),
        selected = "PFS"
      ),
      tags$p(
        "The primary endpoint is tested first. The secondary endpoint is ",
        "tested only after the primary endpoint is rejected.",
        class = "input-help"
      ),
      textInput(
        "d_PFS_vec",
        tagList("Target PFS events by look", tags$span(
          "*", class = "required-marker"
        ))
      ),
      tags$p("Enter increasing event targets separated by commas, for ",
             "example: 50, 100.", class = "input-help"),
      numericInput(
        "alpha",
        "One-sided alpha",
        value = 0.025,
        min = 0,
        max = 1
      ),
      tags$p(
        "One-sided Type I error level used to construct endpoint-specific ",
        "efficacy boundaries.",
        class = "input-help"
      ),
      selectInput(
        "alpha_spending_PFS",
        tagList("PFS alpha-spending function", tags$span(
          "*", class = "required-marker"
        )),
        choices = c("Select a spending function" = "", "OF", "Pocock", "HSD"),
        selected = ""
      ),
      tags$p(
        "Select the alpha-spending function used for PFS efficacy ",
        "boundaries.",
        class = "input-help"
      ),
      selectInput(
        "alpha_spending_OS",
        tagList("OS alpha-spending function", tags$span(
          "*", class = "required-marker"
        )),
        choices = c("Select a spending function" = "", "OF", "Pocock", "HSD"),
        selected = ""
      ),
      tags$p(
        "Select the alpha-spending function used for OS efficacy ",
        "boundaries.",
        class = "input-help"
      ),
      conditionalPanel(
        "input.alpha_spending_PFS == 'HSD'",
        textInput(
          "alpha_spending_gamma_PFS",
          tagList("PFS HSD gamma", tags$span(
            "*", class = "required-marker"
          ))
        ),
        tags$p(
          "HSD gamma for PFS. Use a finite value in [-40, 40). Negative ",
          "values spend more alpha later; positive values spend more alpha ",
          "earlier.",
          class = "input-help"
        )
      ),
      conditionalPanel(
        "input.alpha_spending_OS == 'HSD'",
        textInput(
          "alpha_spending_gamma_OS",
          tagList("OS HSD gamma", tags$span(
            "*", class = "required-marker"
          ))
        ),
        tags$p(
          "HSD gamma for OS. Use a finite value in [-40, 40). Negative ",
          "values spend more alpha later; positive values spend more alpha ",
          "earlier.",
          class = "input-help"
        )
      ),
      tags$label(
        "Boundary configuration",
        class = "control-label configuration-label"
      ),
      uiOutput("look_selection"),
      tags$p(
        "The final efficacy look is always included and preselected; select ",
        "any additional efficacy looks. Futility analyses are off by default ",
        "and available before the final look only. For a selected futility ",
        "look, set an HR threshold above 1 (default 1.2). HSD requires gamma.",
        class = "input-help configuration-help"
      ),
      tags$hr(),
      tags$h3("Advanced settings", class = "sidebar-section-title"),
      checkboxInput(
        "simulation",
        "Run simulation"
      ),
      tags$p(
        "Run a Monte Carlo simulation to compare simulation estimates with ",
        "closed-form estimates. This increases runtime.",
        class = "input-help"
      ),
      conditionalPanel(
        "input.simulation",
        numericInput(
          "n_sim",
          "Monte Carlo replicates",
          value = 500,
          min = 2,
          step = 1
        ),
        tags$p(
          "Number of valid Monte Carlo replicates used for the empirical ",
          "power estimates.",
          class = "input-help"
        ),
        numericInput(
          "simulation_seed",
          "Simulation random seed",
          value = 1,
          min = 0,
          step = 1
        ),
        tags$p(
          "Seed for reproducible simulation results. Changing it changes ",
          "the simulated sample.",
          class = "input-help"
        )
      ),
      checkboxInput(
        "show_computation_settings",
        "Computational settings"
      ),
      tags$p(
        "Show numerical settings. The defaults are suitable for routine ",
        "design evaluation.",
        class = "input-help"
      ),
      conditionalPanel(
        "input.show_computation_settings",
        numericInput(
          "tol",
          "Numerical tolerance",
          value = 1e-8,
          min = 0
        ),
        tags$p(
          "Tolerance used when solving calendar cutoff times. Smaller values ",
          "request greater numerical precision.",
          class = "input-help"
        ),
        numericInput(
          "integration_seed",
          "Integration random seed",
          value = 1,
          min = 0,
          step = 1
        ),
        tags$p(
          "Seed for reproducible multivariate-normal integration in the ",
          "closed-form calculation.",
          class = "input-help"
        ),
        numericInput(
          "display_digits",
          "Displayed decimal places",
          value = 4,
          min = 0,
          max = 10,
          step = 1
        ),
        tags$p(
          "Number of decimal places shown in result tables and CSV exports. ",
          "This does not change calculation precision.",
          class = "input-help"
        )
      ),
      uiOutput("run_control")
    ),
    mainPanel(
      width = 8,
      uiOutput("status"),
      uiOutput("download_control"),
      tabsetPanel(
        id = "main_tabs",
        tabPanel(
          "Introduction",
          div(
            class = "intro-hero",
            tags$div(
              "Clinical trial design evaluation",
              class = "intro-eyebrow"
            ),
            tags$h2("Closed-form evaluation for PFS and OS designs"),
            tags$p(
              "CorrSurvGSD evaluates correlated group sequential designs for "
                , "progression-free survival (PFS) and overall survival (OS) "
                , "under the Fleischer model.",
              class = "intro-lead"
            )
          ),
          div(
            class = "intro-grid",
            div(
              class = "intro-card",
              tags$h4("1. Configure"),
              tags$p(
                "Define sample sizes, control-group medians, accrual, planned "
                  , "PFS events, and alpha-spending choices."
              )
            ),
            div(
              class = "intro-card",
              tags$h4("2. Evaluate"),
              tags$p(
                "Specify treatment effects using hazard ratios or treatment "
                  , "medians, then run the closed-form pipeline."
              )
            ),
            div(
              class = "intro-card",
              tags$h4("3. Review"),
              tags$p(
                "Review design assumptions, decision boundaries, and power "
                  , "summaries. Compare simulation when it is enabled."
              )
            )
          ),
        ),
        tabPanel(
          "Results",
          uiOutput("results_panel")
        )
      )
    )
  )
)

server <- function(input, output, session) {
  required_inputs <- reactive({
    required_input_labels(input)
  })

  output$run_control <- renderUI({
    missing <- required_inputs()

    if (length(missing) > 0) {
      completion_message <- tags$div(
        class = "required-inputs",
        tags$p(
          "Complete the required trial-design inputs before running the ",
          "pipeline."
        ),
        tags$ul(lapply(names(missing), function(id) {
          tags$li(tags$a(
            missing[[id]],
            href = paste0("#", id),
            class = "required-input-link"
          ))
        }))
      )
    } else {
      completion_message <- tags$div(
        class = "input-complete",
        "All required trial-design inputs are complete. Run the pipeline."
      )
    }

    tagList(
      completion_message,
      actionButton(
        "run",
        "Run pipeline",
        class = "btn-primary app-run-button",
        disabled = if (length(missing) > 0) "disabled" else NULL
      )
    )
  })

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
              result <- do.call(CorrSurvGSD::run_pipeline, arguments)
              incProgress(0.9, detail = "Preparing results")
              result
            }
          )
          list(state = state, error = NULL)
        },
        error = function(error) {
          error_message <- conditionMessage(error)
          list(
            state = NULL,
            error = translate_pipeline_error(error_message),
            error_input_ids = pipeline_error_input_ids(error_message)
          )
        }
      )
    },
    ignoreInit = TRUE
  )

  last_successful_state <- reactiveVal(NULL)

  current_state <- reactive({
    state <- last_successful_state()
    req(!is.null(state))
    state
  })

  observeEvent(calculation(), {
    calculation_result <- calculation()
    error_input_ids <- if (is.null(calculation_result$error)) {
      character(0)
    } else {
      calculation_result$error_input_ids
    }

    session$sendCustomMessage("corrsurvInputError", error_input_ids)
    if (is.null(calculation_result$error)) {
      last_successful_state(calculation_result$state)
      updateTabsetPanel(session, "main_tabs", selected = "Results")
    }
  }, ignoreNULL = TRUE)

  empty_state <- function(message, error = FALSE) {
    class_name <- if (error) {
      "empty-state empty-state-error"
    } else {
      "empty-state"
    }
    tags$div(class = class_name, tags$p(message))
  }

  output$status <- renderUI({
    calculation_result <- calculation()

    if (is.null(calculation_result)) {
      return(tags$div(
        class = "app-status",
        "Complete the required trial-design inputs, then run the pipeline."
      ))
    }

    if (!is.null(calculation_result$error)) {
      return(tags$div(
        class = "app-status app-status-error",
        calculation_result$error
      ))
    }

    tags$div(
      class = "app-status app-status-success",
      "Pipeline completed. Review the results or download the tables."
    )
  })

  output$download_control <- renderUI({
    is_complete <- !is.null(last_successful_state())

    if (!is_complete) {
      return(tags$button(
        "Download tables (CSV)",
        class = "btn btn-default app-download-button",
        disabled = "disabled"
      ))
    }

    downloadButton(
      "download_tables",
      "Download tables (CSV)",
      class = "app-download-button"
    )
  })

  output$results_panel <- renderUI({
    state <- last_successful_state()

    if (is.null(state)) {
      return(empty_state(
        "Run the pipeline to view marginal and joint power results."
      ))
    }

    tagList(
      uiOutput("result_context"),
      div(
        class = "result-sections",
        tabsetPanel(
          id = "result_sections",
          tabPanel(
            "Design summary",
            tags$div("Trial characteristics", class = "result-section-label"),
            div(class = "table-responsive", tableOutput("trial_design")),
            tags$div("Accrual schedule", class = "result-section-label"),
            div(class = "table-responsive", tableOutput("accrual_design")),
            tags$div("Analysis schedule", class = "result-section-label"),
            div(class = "table-responsive", tableOutput("analysis_schedule"))
          ),
          tabPanel(
            "Boundaries",
            tags$div("Testing framework", class = "result-section-label"),
            div(
              class = "table-responsive",
              tableOutput("testing_framework")
            ),
            tags$div("Boundary scale", class = "result-section-label"),
            radioButtons(
              "boundary_scale",
              NULL,
              choices = c(
                "Z-score" = "z",
                "Hazard ratio" = "hr",
                "One-sided p-value" = "p"
              ),
              selected = "z",
              inline = TRUE
            ),
            tags$p(
              "Futility and efficacy boundaries are shown on the selected "
                , "scale.",
              class = "app-help"
            ),
            div(class = "table-responsive", tableOutput("boundary_results"))
          ),
          tabPanel(
            "Power",
            tags$div("Power summary", class = "result-section-label"),
            radioButtons(
              "power_view",
              NULL,
              choices = c(
                "Marginal power" = "marginal",
                "Joint power" = "joint"
              ),
              selected = "marginal",
              inline = TRUE
            ),
            uiOutput("power_details")
          )
        )
      )
    )
  })

  output$download_tables <- downloadHandler(
    filename = function() {
      paste0("CorrSurvGSD_results_", format(Sys.time(), "%Y-%m-%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      state <- current_state()
      digits <- state$options$display_digits
      connection <- base::file(file, open = "wt", encoding = "UTF-8")
      on.exit(close(connection), add = TRUE)

      sections <- list(
        "Trial design" = trial_design_table(state),
        "Accrual" = accrual_design_table(state),
        "Analysis schedule" = analysis_schedule_table(state),
        "Testing summary" = testing_summary_table(state),
        "Testing configuration" = endpoint_testing_table(state),
        "Marginal power" = marginal_power_results_table(state),
        "Closed-form joint power" = joint_power_results_table(
          state$theoretical_results$joint_power_matrix,
          digits = digits
        )
      )

      if (isTRUE(state$options$simulation)) {
        sections[["Simulation joint power"]] <- joint_power_results_table(
          state$empirical_results$joint_power_matrix,
          digits = digits
        )
      }

      for (title in names(sections)) {
        write_csv_section(connection, title, sections[[title]])
      }
    }
  )

  render_result_table <- function(table_function) {
    renderTable(
      table_function(current_state()),
      striped = TRUE,
      bordered = TRUE,
      rownames = FALSE
    )
  }

  output$trial_design <- render_result_table(trial_design_table)
  output$accrual_design <- render_result_table(accrual_design_table)
  output$analysis_schedule <- render_result_table(analysis_schedule_table)
  output$result_context <- renderUI({
    tags$div(
      tags$div("Overview", class = "result-context-title"),
      tags$div(
        result_context_text(current_state()),
        class = "result-context-details"
      ),
      class = "result-context"
    )
  })
  output$testing_framework <- render_result_table(testing_framework_details)
  output$boundary_results <- render_result_table(function(state) {
    scale <- if (is.null(input$boundary_scale)) "z" else input$boundary_scale
    boundary_results_table(state, scale = scale)
  })
  output$marginal_power_results <- render_result_table(
    marginal_power_results_table
  )

  output$power_details <- renderUI({
    state <- current_state()
    power_view <- if (is.null(input$power_view)) "marginal" else {
      input$power_view
    }
    comparison_note <- if (isTRUE(state$options$simulation)) {
      "The value in parentheses is the simulation estimate."
    } else {
      NULL
    }

    if (identical(power_view, "marginal")) {
      return(tagList(
        tags$p(
          "Endpoint-level efficacy probability by analysis look. ",
          comparison_note,
          class = "app-help"
        ),
        div(
          class = "table-responsive",
          tableOutput("marginal_power_results")
        )
      ))
    }

    tagList(
      tags$p(
        "Rows identify the primary endpoint's first efficacy crossing, and "
          , "columns identify the secondary endpoint's first efficacy "
          , "crossing. Each cell is the joint probability for that pair of "
          , "looks; - marks an impossible ordering. ",
        comparison_note,
        class = "app-help"
      ),
      div(class = "table-responsive", tableOutput("joint_power_results"))
    )
  })

  output$joint_power_results <- renderTable({
    state <- current_state()
    simulation_power <- if (isTRUE(state$options$simulation)) {
      state$empirical_results$joint_power_matrix
    } else {
      NULL
    }

    joint_power_results_table(
      state$theoretical_results$joint_power_matrix,
      simulation_power = simulation_power,
      digits = state$options$display_digits
    )
  }, striped = TRUE, bordered = TRUE, rownames = FALSE)
}

shinyApp(ui = ui, server = server)
