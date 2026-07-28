# pipeline.R — end-to-end theoretical and simulation pipeline.

# -----------------------------------------------------------------
# run_pipeline — run Modules 0 to 5 and print simulation results
#
# All arguments are forwarded to validate_trial_input().

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
