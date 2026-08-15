# validate_trial_input.R — input validation and derived parameters.
#
# Returns a `trial_state` list with `design` and `options`.
validate_trial_input <- function(n_T, n_C, d_PFS_vec,
                                 t_vec, v_vec,
                                 median_PFS_C, median_OS_C,
                                 median_PFS_T, median_OS_T,
                                 HR_PFS, HR_OS,
                                 alpha_spending_PFS,
                                 alpha_spending_OS,
                                 alpha_spending_gamma_PFS,
                                 alpha_spending_gamma_OS,
                                 alpha,
                                 efficacy_looks,
                                 futility_looks,
                                 futility_HR,
                                 hierarchy_order,
                                 tol,
                                 integration_seed,
                                 simulation_seed,
                                 simulation,
                                 n_sim,
                                 display_digits) {
  # Validate sample sizes.
  checkmate::assert_count(n_T, positive = TRUE, .var.name = "n_T")
  checkmate::assert_count(n_C, positive = TRUE, .var.name = "n_C")

  # Validate control-arm medians and derive control hazards.
  checkmate::assert_number(median_PFS_C, lower = .Machine$double.eps, finite = TRUE,
                .var.name = "median_PFS_C")
  checkmate::assert_number(median_OS_C, lower = .Machine$double.eps, finite = TRUE,
                .var.name = "median_OS_C")
  if (median_PFS_C >= median_OS_C) {
    stop("median_PFS_C must be less than median_OS_C.")
  }

  Lambda_C   <- log(2) / median_PFS_C
  lambda_2_C <- log(2) / median_OS_C
  lambda_1_C <- Lambda_C - lambda_2_C

  has_treatment_medians <- !is.null(median_PFS_T) || !is.null(median_OS_T)
  has_treatment_hr      <- !is.null(HR_PFS) || !is.null(HR_OS)

  if (has_treatment_medians && has_treatment_hr) {
    stop("Provide either treatment medians or HR_PFS/HR_OS, not both.")
  }
  if (!has_treatment_medians && !has_treatment_hr) {
    stop("Provide either treatment medians or HR_PFS/HR_OS.")
  }

  if (has_treatment_medians) {
    checkmate::assert_number(median_PFS_T, lower = .Machine$double.eps, finite = TRUE,
                  .var.name = "median_PFS_T")
    checkmate::assert_number(median_OS_T, lower = .Machine$double.eps, finite = TRUE,
                  .var.name = "median_OS_T")
    if (median_PFS_T >= median_OS_T) {
      stop("median_PFS_T must be less than median_OS_T.")
    }

    Lambda_T   <- log(2) / median_PFS_T
    lambda_2_T <- log(2) / median_OS_T
    lambda_1_T <- Lambda_T - lambda_2_T
    HR_PFS     <- Lambda_T / Lambda_C
    HR_OS      <- lambda_2_T / lambda_2_C
  } else {
    checkmate::assert_number(HR_PFS, lower = .Machine$double.eps, upper = 1,
                  finite = TRUE,
                  .var.name = "HR_PFS")
    checkmate::assert_number(HR_OS, lower = .Machine$double.eps, upper = 1,
                  finite = TRUE,
                  .var.name = "HR_OS")

    Lambda_T   <- HR_PFS * Lambda_C
    lambda_2_T <- HR_OS * lambda_2_C
    lambda_1_T <- Lambda_T - lambda_2_T
    if (lambda_1_T <= 0) {
      stop("lambda_1_T must be positive after applying HR_PFS and HR_OS.")
    }
  }

  # Validate accrual rates.
  checkmate::assert_numeric(v_vec, lower = .Machine$double.eps,
                 any.missing = FALSE,
                 .var.name = "v_vec")

  # Validate accrual segment boundaries.
  n     <- n_T + n_C
  s     <- length(v_vec)

  if (s == 1) {
    if (!is.null(t_vec))
      stop("t_vec must be NULL for single-segment (s=1) accrual.")
    R_vec <- n / v_vec[1]
  } else {
    if (is.null(t_vec) || length(t_vec) != s - 1)
      stop("t_vec must have length s-1 = ", s - 1,
           " for ", s, " segments.")
    checkmate::assert_numeric(t_vec, lower = .Machine$double.eps,
                   any.missing = FALSE,
                   .var.name = "t_vec")
    if (is.unsorted(t_vec, strictly = TRUE))
      stop("t_vec must be strictly increasing.")

    R_vec <- numeric(s)
    R_vec[1] <- t_vec[1]
    if (s > 2)
      for (i in 2:(s - 1))
        R_vec[i] <- t_vec[i] - t_vec[i - 1]
    R_vec[s] <- (n - sum(v_vec[1:(s - 1)] * R_vec[1:(s - 1)])) / v_vec[s]
    if (R_vec[s] <= 0)
      stop("R_s must be positive; check t_vec and v_vec.")
  }

  t_vec <- cumsum(R_vec)

  # Derive accrual quantities.
  R     <- sum(R_vec)
  p_vec <- v_vec / n

  # Check the total accrued sample size.
  if (abs(sum(v_vec * R_vec) - n) > 1e-8 * n)
    stop("sum(v_i * R_i) must equal n (total sample size).")

  # Validate target PFS event counts.
  checkmate::assert_numeric(d_PFS_vec, lower = 1, any.missing = FALSE,
                 min.len = 1, sorted = TRUE,
                 .var.name = "d_PFS_vec")
  checkmate::assert_integerish(d_PFS_vec, .var.name = "d_PFS_vec")

  if (any(duplicated(d_PFS_vec))) {
    stop("Assertion on 'd_PFS_vec' failed: ",
         "Must be strictly increasing (no duplicate values)")
  }

  if (any(d_PFS_vec > n))
    stop("d_PFS_ell must not exceed total sample size n.")

  L <- length(d_PFS_vec)

  # Validate endpoint look schedules.
  endpoints <- c("PFS", "OS")
  if (is.null(efficacy_looks)) {
    efficacy_looks <- list(PFS = L, OS = L)
  }
  if (is.null(futility_looks)) {
    futility_looks <- list(PFS = integer(0), OS = integer(0))
  }
  if (!is.list(efficacy_looks) || length(efficacy_looks) != 2 ||
      !identical(names(efficacy_looks), endpoints)) {
    stop("efficacy_looks must be a named list with elements PFS and OS.")
  }
  if (!is.list(futility_looks) || length(futility_looks) != 2 ||
      !identical(names(futility_looks), endpoints)) {
    stop("futility_looks must be a named list with elements PFS and OS.")
  }

  for (endpoint in endpoints) {
    endpoint_efficacy_looks <- efficacy_looks[[endpoint]]
    if (!is.numeric(endpoint_efficacy_looks) ||
        anyNA(endpoint_efficacy_looks) ||
        any(!is.finite(endpoint_efficacy_looks)) ||
        any(endpoint_efficacy_looks != as.integer(endpoint_efficacy_looks)) ||
        any(endpoint_efficacy_looks < 1 | endpoint_efficacy_looks > L) ||
        is.unsorted(endpoint_efficacy_looks, strictly = TRUE) ||
        length(endpoint_efficacy_looks) == 0 ||
        utils::tail(endpoint_efficacy_looks, 1) != L) {
      stop(
        "efficacy_looks must contain strictly increasing looks and include L."
      )
    }
    efficacy_looks[[endpoint]] <- as.integer(endpoint_efficacy_looks)

    endpoint_futility_looks <- futility_looks[[endpoint]]
    if (!is.numeric(endpoint_futility_looks) ||
        anyNA(endpoint_futility_looks) ||
        any(!is.finite(endpoint_futility_looks)) ||
        any(endpoint_futility_looks != as.integer(endpoint_futility_looks)) ||
        any(endpoint_futility_looks < 1 |
            endpoint_futility_looks > L - 1) ||
        is.unsorted(endpoint_futility_looks, strictly = TRUE)) {
      stop(
        "futility_looks must contain strictly increasing looks in 1:(L - 1)."
      )
    }
    futility_looks[[endpoint]] <- as.integer(endpoint_futility_looks)
  }

  # Validate and expand futility hazard ratios.
  if (is.numeric(futility_HR) && length(futility_HR) == 1) {
    if (is.na(futility_HR) || !is.finite(futility_HR) ||
        futility_HR <= 1) {
      stop("futility_HR must be finite and greater than 1.")
    }
    futility_HR <- lapply(futility_looks, function(looks) {
      rep(futility_HR, length(looks))
    })
  } else if (is.numeric(futility_HR) && length(futility_HR) == 2) {
    if (!identical(names(futility_HR), endpoints) || anyNA(futility_HR) ||
        any(!is.finite(futility_HR)) || any(futility_HR <= 1)) {
      stop(
        "Length-two futility_HR must be named PFS and OS with values above 1."
      )
    }
    endpoint_hr <- futility_HR
    futility_HR <- lapply(endpoints, function(endpoint) {
      rep(endpoint_hr[[endpoint]], length(futility_looks[[endpoint]]))
    })
    names(futility_HR) <- endpoints
  } else if (is.list(futility_HR) && length(futility_HR) == 2) {
    if (!identical(names(futility_HR), endpoints)) {
      stop("List futility_HR must be named PFS and OS.")
    }
    for (endpoint in endpoints) {
      endpoint_hr <- futility_HR[[endpoint]]
      if (!is.numeric(endpoint_hr) ||
          length(endpoint_hr) != length(futility_looks[[endpoint]]) ||
          anyNA(endpoint_hr) || any(!is.finite(endpoint_hr)) ||
          any(endpoint_hr <= 1)) {
        stop(
          "Each futility_HR list element must match its futility_looks length."
        )
      }
    }
  } else {
    stop(
      "futility_HR must be a scalar, a named endpoint vector, or a named list."
    )
  }

  if (!is.character(hierarchy_order) || length(hierarchy_order) != 2 ||
      !identical(names(hierarchy_order), c("primary", "secondary")) ||
      !setequal(hierarchy_order, endpoints)) {
    stop("hierarchy_order must name primary and secondary, with values PFS and OS.")
  }

  # Validate alpha-spending configurations.
  checkmate::assert_choice(alpha_spending_PFS, c("OF", "Pocock", "HSD"),
                .var.name = "alpha_spending_PFS")
  checkmate::assert_choice(alpha_spending_OS, c("OF", "Pocock", "HSD"),
                .var.name = "alpha_spending_OS")
  alpha_spending <- c(PFS = alpha_spending_PFS, OS = alpha_spending_OS)
  alpha_spending_gamma <- list(
    PFS = alpha_spending_gamma_PFS,
    OS = alpha_spending_gamma_OS
  )
  for (endpoint in endpoints) {
    if (identical(alpha_spending[[endpoint]], "HSD")) {
      endpoint_gamma <- alpha_spending_gamma[[endpoint]]
      if (!is.numeric(endpoint_gamma) || length(endpoint_gamma) != 1 ||
          is.na(endpoint_gamma) || !is.finite(endpoint_gamma) ||
          endpoint_gamma < -40 || endpoint_gamma >= 40) {
        stop(
          "alpha_spending_gamma must be finite and in the interval [-40, 40)."
        )
      }
    } else {
      alpha_spending_gamma[[endpoint]] <- NA_real_
    }
  }

  # Validate the type-I error level.
  checkmate::assert_number(alpha, lower = 0, upper = 1, finite = TRUE,
                .var.name = "alpha")
  if (alpha <= 0 || alpha >= 1) {
    stop("alpha must be strictly between 0 and 1.")
  }

  # Validate numerical and simulation options.
  checkmate::assert_number(tol, lower = 0, finite = TRUE,
                .var.name = "tol")

  checkmate::assert_count(integration_seed, positive = FALSE,
               .var.name = "integration_seed")
  checkmate::assert_count(simulation_seed, positive = FALSE,
               .var.name = "simulation_seed")

  checkmate::assert_flag(simulation, .var.name = "simulation")
  if (simulation) {
    checkmate::assert_count(n_sim, positive = TRUE, .var.name = "n_sim")
    if (n_sim < 2) {
      stop("n_sim must be at least 2 when simulation is enabled.")
    }
  }

  checkmate::assert_count(display_digits, positive = FALSE,
               .var.name = "display_digits")
  if (display_digits > 10) {
    stop("display_digits must be at most 10.")
  }

  # Assemble the validated state.
  r          <- n_T / n_C

  state <- list(
    design = list(
      n_T        = n_T,
      n_C        = n_C,
      n          = n,
      R          = R,
      R_vec      = R_vec,
      t_vec      = t_vec,
      v_vec      = v_vec,
      p_vec      = p_vec,
      lambda_1_T = lambda_1_T,
      lambda_2_T = lambda_2_T,
      lambda_1_C = lambda_1_C,
      lambda_2_C = lambda_2_C,
      Lambda_T   = Lambda_T,
      Lambda_C   = Lambda_C,
      d_PFS_vec  = d_PFS_vec,
      L          = L,
      r          = r,
      HR_PFS     = HR_PFS,
      HR_OS      = HR_OS,
      efficacy_looks    = efficacy_looks,
      futility_looks    = futility_looks,
      futility_HR       = futility_HR,
      hierarchy_order   = hierarchy_order,
      alpha_spending_PFS = alpha_spending_PFS,
      alpha_spending_OS  = alpha_spending_OS,
      alpha_spending_gamma_PFS = alpha_spending_gamma[["PFS"]],
      alpha_spending_gamma_OS  = alpha_spending_gamma[["OS"]],
      alpha      = alpha
    ),
    options = list(
      tol              = tol,
      integration_seed = integration_seed,
      simulation_seed  = simulation_seed,
      simulation       = simulation,
      n_sim            = n_sim,
      display_digits   = display_digits
    )
  )
  class(state) <- "trial_state"
  state
}
