# Shared fixtures for CorrSurvGSD tests.

validate_test_input <- function(...) {
  supplied <- list(...)
  d_PFS_vec <- if (is.null(supplied$d_PFS_vec)) {
    c(50)
  } else {
    supplied$d_PFS_vec
  }

  defaults <- list(
    n_T = 100,
    n_C = 100,
    d_PFS_vec = d_PFS_vec,
    t_vec = NULL,
    v_vec = c(200 / 12),
    median_PFS_C = log(2) / 0.15,
    median_OS_C = log(2) / 0.05,
    median_PFS_T = NULL,
    median_OS_T = NULL,
    HR_PFS = NULL,
    HR_OS = NULL,
    alpha_spending_PFS = "OF",
    alpha_spending_OS = "Pocock",
    alpha_spending_gamma_PFS = NA_real_,
    alpha_spending_gamma_OS = NA_real_,
    alpha = 0.025,
    efficacy_looks = list(
      PFS = seq_along(d_PFS_vec),
      OS = seq_along(d_PFS_vec)
    ),
    futility_looks = list(PFS = integer(0), OS = integer(0)),
    futility_HR = 1.2,
    hierarchy_order = c(primary = "PFS", secondary = "OS"),
    tol = 1e-8,
    integration_seed = 1,
    simulation_seed = 1,
    simulation = FALSE,
    n_sim = 500,
    display_digits = 4
  )

  do.call(validate_trial_input, modifyList(defaults, supplied, keep.null = TRUE))
}

make_test_state <- function(stage = "correlation", ...) {
  stage <- match.arg(stage, c("validated", "moments", "correlation", "closed"))
  state <- validate_test_input(...)

  if (stage %in% c("moments", "correlation", "closed")) {
    state <- calculate_calendar_cutoffs(state)
  }
  if (stage %in% c("moments", "correlation", "closed")) {
    state <- compute_per_subject_moments(state)
  }
  if (stage %in% c("correlation", "closed")) {
    state <- construct_joint_correlation_matrix(state)
  }
  if (stage == "closed") {
    state <- calculate_closed_form_results(state)
  }

  state
}

local_moments_fixture <- function() {
  state <- make_test_state(
    stage = "moments",
    n_T = 80,
    n_C = 120,
    d_PFS_vec = c(30, 60),
    t_vec = c(6),
    v_vec = c(200 / 12, 200 / 12),
    median_PFS_C = log(2) / 0.25,
    median_OS_C = log(2) / 0.10,
    median_PFS_T = log(2) / 0.08,
    median_OS_T = log(2) / 0.03,
    alpha_spending_PFS = "OF",
    alpha_spending_OS = "Pocock"
  )
  design <- state$design
  fixture <- list(
    state = state,
    design = design,
    moments = state$moments,
    A_vec = design$A_vec,
    R = design$R,
    t_vec = design$t_vec,
    p_vec = design$p_vec,
    L = design$L,
    groups = c("T", "C"),
    types = c("delta", "time"),
    eps = c("PFS", "OS")
  )
  fixture
}

local_correlation_fixture <- function() {
  state <- make_test_state(
    n_T = 80,
    n_C = 120,
    d_PFS_vec = c(30, 60),
    t_vec = c(6),
    v_vec = c(200 / 12, 200 / 12),
    median_PFS_C = log(2) / 0.25,
    median_OS_C = log(2) / 0.10,
    median_PFS_T = log(2) / 0.08,
    median_OS_T = log(2) / 0.03,
    alpha_spending_PFS = "OF",
    alpha_spending_OS = "Pocock"
  )
  design <- state$design
  moments <- state$moments
  L <- design$L
  sm <- moments$single_look_means
  d_OS <- vapply(
    seq_len(L),
    function(ell) {
      design$n_T * sm["T", "delta", "OS", ell] +
        design$n_C * sm["C", "delta", "OS", ell]
    },
    numeric(1)
  )
  fixture <- list(
    state = state,
    design = design,
    moments = moments,
    R_mat = state$theoretical_results$joint_correlation_matrix,
    L = L,
    d_PFS = design$d_PFS_vec,
    n_T = design$n_T,
    n_C = design$n_C,
    sm = sm,
    sc = moments$single_look_covariances,
    cc = moments$cross_look_covariances,
    groups = c("T", "C"),
    d_OS = d_OS
  )
  fixture
}
