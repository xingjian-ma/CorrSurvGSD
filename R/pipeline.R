# pipeline.R — end-to-end theoretical and simulation pipeline.

# -----------------------------------------------------------------
#' Run the complete correlated PFS and OS group sequential analysis.
#'
#' The pipeline validates the trial design, computes calendar cutoffs and
#' per-subject moments, constructs the joint correlation matrix, calculates
#' theoretical power, and optionally runs Monte Carlo simulation.
#'
#' @param ... Arguments forwarded to the trial-design validation stage.
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

run_pipeline <- function(...) {
  state <- validate_trial_input(...)
  state <- calendar_cutoff(state)
  state <- per_subject_moments(state)
  state <- joint_cor_matrix(state)
  state <- closed_gsd_os_and_pfs(state)

  if (state$options$simulation) {
    state <- run_simulation(state)
  }

  state
}
