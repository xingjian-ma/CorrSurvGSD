# run_example.R — execute an example group sequential design.

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("../R/utils.R")
source("../R/per_subject_moments.R")
source("../R/calendar_cutoff.R")
source("../R/joint_cor_matrix.R")
source("../R/closed_gsd_os_and_pfs.R")
source("../R/simulation.R")
source("../R/pipeline.R")

print_pipeline_results <- function(state) {
  cat("\nExpected calendar cutoff times (A_vec):\n")
  print(state$design$A_vec)

  cat("\nTheoretical results:\n")
  print(state$theoretical_results)

  if (state$options$simulation) {
    cat("\nEmpirical results:\n")
    print(state$empirical_results)
  }
}

state <- run_pipeline(
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
  n_sim = 5000,
  seed = 777
)
print_pipeline_results(state)

state <- run_pipeline(
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
  simulation = TRUE,
  n_sim = 5000,
  seed = 777
)
print_pipeline_results(state)

state <- run_pipeline(
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
  simulation = TRUE,
  efficacy_start = c(PFS = 2, OS = 1),
  n_sim = 5000,
  seed = 777
)
print_pipeline_results(state)
