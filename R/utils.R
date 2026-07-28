# utils.R — Module 0: input validation and derived parameters
#
# validate_trial_input() validates all numeric inputs and returns a
# `state` list with element `design` containing derived quantities.
#
# Dependencies: checkmate

library(checkmate)

# -----------------------------------------------------------------
#
# Args:
#   n_T         — treatment arm sample size
#   n_C         — control arm sample size
#   d_PFS_vec   — target PFS event counts per look
#   t_vec       — segment start nodes (NULL for single-segment s=1)
#   v_vec       — accrual rates per segment
#   median_PFS_C, median_OS_C  — control arm median times
#   median_PFS_T, median_OS_T  — treatment arm median times
#   HR_PFS, HR_OS              — hazard ratio under H1 (PFS, OS)
#   alpha_spending_PFS     — spending function, one of "WT" "OF" "Pocock"
#   alpha_spending_OS      — spending function, one of "WT" "OF" "Pocock"
#   alpha                  — overall alpha level (default 0.025)
#   futility_analysis      — named endpoint-level futility indicators
#   futility_HR            — named endpoint-level futility HR thresholds
#   efficacy_start         — named endpoint-level first efficacy looks
#   hierarchy_order        — named primary-to-secondary endpoint order
#   tol                    — shared numerical tolerance setting
#   seed                   — shared random seed for numerical routines
#   simulation             — whether to run Module 5 simulation
#   n_sim                  — number of valid simulation replicates
#
# Returns a list with class "trial_state", containing elements
# `design` and `options`.
validate_trial_input <- function(n_T, n_C, d_PFS_vec,
                                 t_vec = NULL, v_vec,
                                 median_PFS_C, median_OS_C,
                                 median_PFS_T = NULL, median_OS_T = NULL,
                                 HR_PFS = NULL, HR_OS = NULL,
                                 alpha_spending_PFS,
                                 alpha_spending_OS,
                                 alpha = 0.025,
                                 futility_analysis = c(PFS = TRUE, OS = FALSE),
                                 futility_HR = c(PFS = 1.2, OS = 1.2),
                                 efficacy_start = c(PFS = 1, OS = 1),
                                 hierarchy_order = c(primary = "PFS", secondary = "OS"),
                                 tol = 1e-8,
                                 seed = 777,
                                 simulation = FALSE,
                                 n_sim = 500) {
  # --- n_T, n_C: sample sizes (integer, strictly positive) ---
  assert_count(n_T, positive = TRUE, .var.name = "n_T")
  assert_count(n_C, positive = TRUE, .var.name = "n_C")

  # --- control medians: required, positive, PFS before OS ---
  assert_number(median_PFS_C, lower = .Machine$double.eps, finite = TRUE,
                .var.name = "median_PFS_C")
  assert_number(median_OS_C, lower = .Machine$double.eps, finite = TRUE,
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
    assert_number(median_PFS_T, lower = .Machine$double.eps, finite = TRUE,
                  .var.name = "median_PFS_T")
    assert_number(median_OS_T, lower = .Machine$double.eps, finite = TRUE,
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
    assert_number(HR_PFS, lower = .Machine$double.eps, finite = TRUE,
                  .var.name = "HR_PFS")
    assert_true(HR_PFS < 1)
    assert_number(HR_OS, lower = .Machine$double.eps, finite = TRUE,
                  .var.name = "HR_OS")
    assert_true(HR_OS < 1)

    Lambda_T   <- HR_PFS * Lambda_C
    lambda_2_T <- HR_OS * lambda_2_C
    lambda_1_T <- Lambda_T - lambda_2_T
    if (lambda_1_T <= 0) {
      stop("lambda_1_T must be positive after applying HR_PFS and HR_OS.")
    }
  }

  # --- v_vec: per-segment accrual rates ---
  assert_numeric(v_vec, lower = .Machine$double.eps,
                 any.missing = FALSE,
                 .var.name = "v_vec")

  # --- t_vec: segment start nodes t_1,...,t_{s-1} ---
  #                             (NULL for single-segment s=1)
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
    assert_numeric(t_vec, lower = .Machine$double.eps,
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

  # --- derived accrual quantities ---
  R     <- sum(R_vec)
  p_vec <- v_vec / n

  # --- constraint check ---
  if (abs(sum(v_vec * R_vec) - n) > 1e-8 * n)
    stop("sum(v_i * R_i) must equal n (total sample size).")

  # --- d_PFS_vec: target PFS event counts ---
  #     integer, at least 1 look, strictly increasing
  assert_numeric(d_PFS_vec, lower = 1, any.missing = FALSE,
                 min.len = 1, sorted = TRUE,
                 .var.name = "d_PFS_vec")
  assert_integerish(d_PFS_vec, .var.name = "d_PFS_vec")

  if (any(duplicated(d_PFS_vec))) {
    stop("Assertion on 'd_PFS_vec' failed: ",
         "Must be strictly increasing (no duplicate values)")
  }

  if (any(d_PFS_vec > n))
    stop("d_PFS_ell must not exceed total sample size n.")

  L <- length(d_PFS_vec)

  # --- endpoint-level futility and efficacy configuration ---
  endpoints <- c("PFS", "OS")
  if (!is.logical(futility_analysis) || length(futility_analysis) != 2 ||
      !identical(names(futility_analysis), endpoints) || anyNA(futility_analysis)) {
    stop("futility_analysis must be a named logical vector with names PFS and OS.")
  }
  if (!is.numeric(futility_HR) || length(futility_HR) != 2 ||
      !identical(names(futility_HR), endpoints)) {
    stop("futility_HR must be a named numeric vector with names PFS and OS.")
  }
  if (!is.numeric(efficacy_start) || length(efficacy_start) != 2 ||
      !identical(names(efficacy_start), endpoints) || anyNA(efficacy_start) ||
      any(efficacy_start != as.integer(efficacy_start)) ||
      any(efficacy_start < 1 | efficacy_start > L)) {
    stop("efficacy_start must be a named integer vector with values in 1:L.")
  }
  efficacy_start <- as.integer(efficacy_start)
  names(efficacy_start) <- endpoints

  if (!is.character(hierarchy_order) || length(hierarchy_order) != 2 ||
      !identical(names(hierarchy_order), c("primary", "secondary")) ||
      !setequal(hierarchy_order, endpoints)) {
    stop("hierarchy_order must name primary and secondary, with values PFS and OS.")
  }

  for (endpoint in endpoints) {
    if (futility_analysis[[endpoint]]) {
      if (is.na(futility_HR[[endpoint]]) ||
          !is.finite(futility_HR[[endpoint]]) ||
          futility_HR[[endpoint]] <= 1) {
        stop("futility_HR must be finite and greater than 1 when futility is enabled.")
      }
    } else {
      futility_HR[[endpoint]] <- NA_real_
      if (efficacy_start[[endpoint]] != 1) {
        stop("efficacy_start must be 1 when futility is disabled.")
      }
    }
  }

  # --- alpha_spending_PFS: spending function string for PFS ---
  assert_choice(alpha_spending_PFS, c("WT", "OF", "Pocock"), .var.name = "alpha_spending_PFS")

  # --- alpha_spending_OS: spending function string for OS ---
  assert_choice(alpha_spending_OS, c("WT", "OF", "Pocock"), .var.name = "alpha_spending_OS")

  # --- alpha: overall alpha level ---
  assert_number(alpha, lower = 0, upper = 1,
                .var.name = "alpha")

  # --- shared numerical options ---
  assert_number(tol, lower = 0, finite = TRUE,
                .var.name = "tol")

  assert_count(seed, positive = FALSE, .var.name = "seed")

  assert_flag(simulation, .var.name = "simulation")
  if (simulation) {
    assert_count(n_sim, positive = TRUE, .var.name = "n_sim")
    if (n_sim < 2) {
      stop("n_sim must be at least 2 when simulation is enabled.")
    }
  }

  # --- assemble and return ---
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
      futility_analysis = futility_analysis,
      futility_HR       = futility_HR,
      efficacy_start    = efficacy_start,
      hierarchy_order   = hierarchy_order,
      alpha_spending_PFS = alpha_spending_PFS,
      alpha_spending_OS  = alpha_spending_OS,
      alpha      = alpha
    ),
    options = list(
      tol        = tol,
      seed       = seed,
      simulation = simulation,
      n_sim      = n_sim
    )
  )
  class(state) <- "trial_state"
  state
}
