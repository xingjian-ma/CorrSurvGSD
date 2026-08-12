# run_pipeline.R — end-to-end theoretical and simulation pipeline.

#' Run the complete correlated PFS and OS group sequential analysis.
#'
#' The pipeline validates the trial design, computes calendar cutoffs and
#' per-subject moments, constructs the joint correlation matrix, calculates
#' theoretical power, and optionally runs Monte Carlo simulation.
#'
#' @param n_T Positive treatment-arm sample size.
#' @param n_C Positive control-arm sample size.
#' @param d_PFS_vec Increasing target PFS event counts by look.
#' @param t_vec Accrual segment end times. Use `NULL` for one segment.
#' @param v_vec Accrual rates for the segments.
#' @param median_PFS_C Control-arm PFS median time.
#' @param median_OS_C Control-arm OS median time.
#' @param median_PFS_T Treatment-arm PFS median time. Use this and
#'   `median_OS_T` instead of the treatment-arm hazard ratios.
#' @param median_OS_T Treatment-arm OS median time.
#' @param HR_PFS Treatment-to-control PFS hazard ratio. Use this and
#'   `HR_OS` instead of the treatment-arm median times.
#' @param HR_OS Treatment-to-control OS hazard ratio.
#' @param alpha_spending_PFS PFS efficacy alpha-spending family: `"OF"`,
#'   `"Pocock"`, or `"HSD"`.
#' @param alpha_spending_OS OS efficacy alpha-spending family: `"OF"`,
#'   `"Pocock"`, or `"HSD"`.
#' @param alpha_spending_gamma_PFS HSD gamma parameter for PFS. Required
#'   only when `alpha_spending_PFS = "HSD"`.
#' @param alpha_spending_gamma_OS HSD gamma parameter for OS. Required only
#'   when `alpha_spending_OS = "HSD"`.
#' @param alpha One-sided type-I error level. Defaults to `0.025`.
#' @param efficacy_looks Named list of PFS and OS efficacy-look schedules.
#'   `NULL` uses every look for both endpoints.
#' @param futility_looks Named list of PFS and OS futility-look schedules.
#'   `NULL` disables futility analysis for both endpoints.
#' @param futility_HR Futility HR specification: one scalar, a named PFS/OS
#'   vector, or a named list of endpoint-specific look vectors.
#' @param hierarchy_order Named primary and secondary endpoint order.
#' @param tol Numerical root-finding tolerance. Defaults to `1e-8`.
#' @param integration_seed Random seed for theoretical integration.
#' @param simulation Whether to run Monte Carlo simulation.
#' @param simulation_seed Random seed for Monte Carlo simulation.
#' @param n_sim Number of valid Monte Carlo replicates when simulation is on.
#' @param display_digits Number of digits used in displayed result tables.
#'
#' @return A `trial_state` object containing the validated design, theoretical
#' results, and, when requested, empirical simulation results.
#'
#' @examples
#' \dontrun{
#' result <- run_pipeline(
#'   n_T = 100,
#'   n_C = 100,
#'   d_PFS_vec = c(50, 100),
#'   v_vec = c(200 / 12),
#'   median_PFS_C = log(2) / 0.15,
#'   median_OS_C = log(2) / 0.05,
#'   median_PFS_T = log(2) / 0.15,
#'   median_OS_T = log(2) / 0.05,
#'   alpha_spending_PFS = "OF",
#'   alpha_spending_OS = "Pocock"
#' )
#' result$theoretical_results$joint_power
#' }
#'
#' @export

run_pipeline <- function(n_T, n_C, d_PFS_vec,
                         t_vec = NULL, v_vec,
                         median_PFS_C, median_OS_C,
                         median_PFS_T = NULL, median_OS_T = NULL,
                         HR_PFS = NULL, HR_OS = NULL,
                         alpha_spending_PFS,
                         alpha_spending_OS,
                         alpha_spending_gamma_PFS = NA_real_,
                         alpha_spending_gamma_OS = NA_real_,
                         alpha = 0.025,
                         efficacy_looks = NULL,
                         futility_looks = NULL,
                         futility_HR = 1.2,
                         hierarchy_order = c(
                           primary = "PFS",
                           secondary = "OS"
                         ),
                         tol = 1e-8,
                         integration_seed = 1,
                         simulation = FALSE,
                         simulation_seed = 1,
                         n_sim = 500,
                         display_digits = 4) {
  state <- validate_trial_input(
    n_T = n_T,
    n_C = n_C,
    d_PFS_vec = d_PFS_vec,
    t_vec = t_vec,
    v_vec = v_vec,
    median_PFS_C = median_PFS_C,
    median_OS_C = median_OS_C,
    median_PFS_T = median_PFS_T,
    median_OS_T = median_OS_T,
    HR_PFS = HR_PFS,
    HR_OS = HR_OS,
    alpha_spending_PFS = alpha_spending_PFS,
    alpha_spending_OS = alpha_spending_OS,
    alpha_spending_gamma_PFS = alpha_spending_gamma_PFS,
    alpha_spending_gamma_OS = alpha_spending_gamma_OS,
    alpha = alpha,
    efficacy_looks = efficacy_looks,
    futility_looks = futility_looks,
    futility_HR = futility_HR,
    hierarchy_order = hierarchy_order,
    tol = tol,
    integration_seed = integration_seed,
    simulation = simulation,
    simulation_seed = simulation_seed,
    n_sim = n_sim,
    display_digits = display_digits
  )
  state <- calculate_calendar_cutoffs(state)
  state <- compute_per_subject_moments(state)
  state <- construct_joint_correlation_matrix(state)
  state <- calculate_closed_form_results(state)

  if (state$options$simulation) {
    state <- run_simulation(state)
  }

  state
}
