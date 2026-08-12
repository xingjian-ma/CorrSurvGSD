# construct_joint_correlation_matrix.R — Module 3: (2L)×(2L) joint correlation matrix
#
# §3.1  Delta method: marginal variance Var(θ̂_ℓ^k) and cross-endpoint
#       covariance Cov(θ̂_ℓ₁^PFS, θ̂_ℓ₂^OS), summed over T and C groups.
# §3.2  Assemble correlation matrix:
#         - Within-endpoint: canonical G-S form √(d_min / d_max)
#         - Cross-endpoint:  Cov / √(Var_PFS · Var_OS)
#       d_OS_vec is computed from moments and written into design.
#
# Dependencies: compute_per_subject_moments.R (moments list)

# -----------------------------------------------------------------
 # §3.1.1  marginal_variances — summed over T and C groups via delta
 #         method.  ∇g = (1/m_δ, -1/m_T̃)^T, Σ from single_look_cov.
 #
 # Returns: log_HR_var array, dims (endpoint, look)

marginal_variances <- function(state) {
  n_T <- state$design$n_T
  n_C <- state$design$n_C
  L   <- state$design$L
  sm  <- state$moments$single_look_means
  sc  <- state$moments$single_look_covariances

  groups <- c("T", "C")
  eps    <- c("PFS", "OS")

  log_HR_var <- array(0, dim = c(2, L),
    dimnames = list(eps, NULL))

  for (g in groups) {
    n_g <- if (g == "T") n_T else n_C
    for (ep in eps) {
      for (ell in seq_len(L)) {
        grad  <- c(1, -1) / sm[g, , ep, ell]
        Sigma <- sc[g, , , ep, ell]
        log_HR_var[ep, ell] <- log_HR_var[ep, ell] +
          as.numeric(t(grad) %*% Sigma %*% grad) / n_g
      }
    }
  }

  log_HR_var
}

# -----------------------------------------------------------------
# §3.1.2  cross_endpoint_covariances — summed over T and C groups
#         via delta method.  ∇g_PFS^T Σ_cross ∇g_OS.
#
# Returns: log_HR_cov matrix, dims (PFS_look, OS_look)

cross_endpoint_covariances <- function(state) {
  n_T <- state$design$n_T
  n_C <- state$design$n_C
  L   <- state$design$L
  sm  <- state$moments$single_look_means
  cc  <- state$moments$cross_look_covariances

  groups <- c("T", "C")

  log_HR_cov <- matrix(0, L, L)

  for (g in groups) {
    n_g <- if (g == "T") n_T else n_C
    for (ell1 in seq_len(L)) {
      for (ell2 in seq_len(L)) {
        grad_PFS <- c(1, -1) / sm[g, , "PFS", ell1]
        grad_OS  <- c(1, -1) / sm[g, , "OS",  ell2]
        Sigma_cross <- cc[g, , , ell1, ell2]
        log_HR_cov[ell1, ell2] <- log_HR_cov[ell1, ell2] +
          as.numeric(t(grad_PFS) %*% Sigma_cross %*% grad_OS) / n_g
      }
    }
  }

  log_HR_cov
}

 # -----------------------------------------------------------------
# §3.2  construct_joint_correlation_matrix — assemble (2L)×(2L) correlation matrix
#
# Args:
#   design  — trial design list from utils → calendar_cutoff pipeline
#   moments — per-subject moments list from compute_per_subject_moments()
#
# Returns:
 #   (2L)×(2L) correlation matrix, dimnames PFS_1…PFS_L, OS_1…OS_L

construct_joint_correlation_matrix <- function(state) {
  L          <- state$design$L
  d_PFS_vec  <- state$design$d_PFS_vec
  d_OS_vec   <- state$design$d_OS_vec

  # --- §3.1.1  marginal variances (summed over groups) --------------
  log_HR_var <- marginal_variances(state)

  # --- §3.1.2  cross-endpoint covariances (summed over groups) ------
  log_HR_cov <- cross_endpoint_covariances(state)

  # --- §3.2  assemble (2L)×(2L) correlation matrix ------------------
  R_mat <- diag(1, 2 * L)
  nm <- c(paste0("PFS_", seq_len(L)), paste0("OS_", seq_len(L)))
  rownames(R_mat) <- colnames(R_mat) <- nm

  # §3.2.1  within-endpoint: canonical G-S form
  for (ell1 in seq_len(L)) {
    for (ell2 in seq_len(L)) {
      if (ell1 > ell2) next

      # PFS block (upper-left)
      rho_pfs <- sqrt(d_PFS_vec[ell1] / d_PFS_vec[ell2])
      R_mat[ell1, ell2] <- rho_pfs
      R_mat[ell2, ell1] <- rho_pfs

      # OS block (lower-right)
      rho_os <- sqrt(d_OS_vec[ell1] / d_OS_vec[ell2])
      R_mat[L + ell1, L + ell2] <- rho_os
      R_mat[L + ell2, L + ell1] <- rho_os
    }
  }

  # §3.2.2  cross-endpoint: delta-method correlation
  for (ell1 in seq_len(L)) {
    for (ell2 in seq_len(L)) {
      cov_val <- log_HR_cov[ell1, ell2]
      var_pfs <- log_HR_var["PFS", ell1]
      var_os  <- log_HR_var["OS",  ell2]
      rho     <- cov_val / sqrt(var_pfs * var_os)

      R_mat[ell1,     L + ell2] <- rho   # upper-right
      R_mat[L + ell2, ell1]     <- rho   # lower-left (symmetric)
    }
  }

  state$theoretical_results$joint_correlation_matrix <- R_mat
  state
}
