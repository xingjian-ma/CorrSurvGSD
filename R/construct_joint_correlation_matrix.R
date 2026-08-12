# construct_joint_correlation_matrix.R — joint correlation matrix.
#
# Dependencies: compute_per_subject_moments.R.

# Calculate marginal log-HR variances.

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

# Calculate PFS/OS cross-endpoint log-HR covariances.

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

# Assemble the joint correlation matrix.

construct_joint_correlation_matrix <- function(state) {
  L          <- state$design$L
  d_PFS_vec  <- state$design$d_PFS_vec
  d_OS_vec   <- state$design$d_OS_vec

  log_HR_var <- marginal_variances(state)

  log_HR_cov <- cross_endpoint_covariances(state)

  R_mat <- diag(1, 2 * L)
  nm <- c(paste0("PFS_", seq_len(L)), paste0("OS_", seq_len(L)))
  rownames(R_mat) <- colnames(R_mat) <- nm

  for (ell1 in seq_len(L)) {
    for (ell2 in seq_len(L)) {
      if (ell1 > ell2) next

      rho_pfs <- sqrt(d_PFS_vec[ell1] / d_PFS_vec[ell2])
      R_mat[ell1, ell2] <- rho_pfs
      R_mat[ell2, ell1] <- rho_pfs

      rho_os <- sqrt(d_OS_vec[ell1] / d_OS_vec[ell2])
      R_mat[L + ell1, L + ell2] <- rho_os
      R_mat[L + ell2, L + ell1] <- rho_os
    }
  }

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
