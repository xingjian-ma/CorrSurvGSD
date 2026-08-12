# run_example.R — execute example group sequential design scenarios.

# Run this script with CorrSurvGSD/ as the working directory.
devtools::load_all(".")

print_pipeline_summary <- function(name, state) {
  digits <- state$options$display_digits

  cat("\n=== ", name, " ===\n", sep = "")
  cat("Calendar cutoffs:\n")
  print(state$design$A_vec)

  cat("Theoretical joint power:\n")
  print(state$theoretical_results$joint_power)

  cat("HR-scale boundaries:\n")
  print(round(state$design$boundary_HR, digits = digits))

  cat("p-scale boundaries:\n")
  print(round(state$design$boundary_p, digits = digits))

  if (isTRUE(state$options$simulation)) {
    cat("Empirical joint power:\n")
    print(state$empirical_results$joint_power)
  }

  invisible(state)
}

base_args <- list(
  n_T = 800,
  n_C = 1200,
  d_PFS_vec = c(400, 800, 1200),
  t_vec = c(6, 10),
  v_vec = c(100 / 12, 100 / 12, 100 / 12),
  median_PFS_C = log(2) / 0.25,
  median_OS_C = log(2) / 0.10,
  HR_PFS = 0.8,
  HR_OS = 0.8,
  alpha_spending_PFS = "OF",
  alpha_spending_OS = "Pocock",
  integration_seed = 777,
  simulation_seed = 777,
  n_sim = 5000
)

scenarios <- list(
  theoretical = base_args,
  simulation = modifyList(
    base_args,
    list(simulation = TRUE)
  ),
  futility_hsd = modifyList(
    base_args,
    list(
      alpha_spending_PFS = "HSD",
      alpha_spending_gamma_PFS = -4,
      simulation = TRUE,
      efficacy_looks = list(PFS = c(2, 3), OS = c(1, 3)),
      futility_looks = list(PFS = 1, OS = c(1, 2)),
      futility_HR = list(PFS = 1.2, OS = c(1.30, 1.15))
    )
  )
)

for (scenario_name in names(scenarios)) {
  state <- do.call(
    CorrSurvGSD::run_pipeline,
    scenarios[[scenario_name]]
  )
  print_pipeline_summary(scenario_name, state)
}
