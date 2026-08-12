# calculate_calendar_cutoffs.R — Module 1: calendar cutoffs via root-finding
#
# For each look, solve n_T·m_δ^PFS(T, A_ℓ) + n_C·m_δ^PFS(C, A_ℓ)
# = d_PFS_ℓ for A_ℓ.  The left-hand side is strictly increasing
# in A_ℓ, so uniroot converges reliably.
# m_δ is computed via compute_mean("delta", ...) from Module 2.
#
# Dependencies: compute_per_subject_moments.R (compute_mean, find_segment)

# -----------------------------------------------------------------
# pfs_mean — per-subject PFS event probability for one group

pfs_mean <- function(A, R, t_vec, p_vec, Lambda) {
  B <- min(A, R)
  b <- find_segment(B, t_vec)
  compute_mean("delta", A, B, b, Lambda, t_vec, p_vec)
}

# -----------------------------------------------------------------
# pfs_event_diff — target function f(A) = LHS − d_target for uniroot

pfs_event_diff <- function(A, n_T, n_C, Lambda_T, Lambda_C,
                           R, t_vec, p_vec, d_target) {
  lhs <- n_T * pfs_mean(A, R, t_vec, p_vec, Lambda_T) +
         n_C * pfs_mean(A, R, t_vec, p_vec, Lambda_C)
  lhs - d_target
}

# -----------------------------------------------------------------
# calculate_calendar_cutoff — solve one A_ℓ

calculate_calendar_cutoff <- function(n_T, n_C, Lambda_T, Lambda_C,
                                   R, t_vec, p_vec, d_target,
                                   A_prev = 0, tol) {
  lower <- A_prev
  upper <- A_prev + 2 * R

  f_upper <- pfs_event_diff(upper, n_T, n_C, Lambda_T, Lambda_C,
                            R, t_vec, p_vec, d_target)
  while (f_upper < 0) {
    upper   <- upper + R
    f_upper <- pfs_event_diff(upper, n_T, n_C, Lambda_T, Lambda_C,
                              R, t_vec, p_vec, d_target)
  }

  result <- stats::uniroot(pfs_event_diff, interval = c(lower, upper),
                    n_T = n_T, n_C = n_C,
                    Lambda_T = Lambda_T, Lambda_C = Lambda_C,
                    R = R, t_vec = t_vec, p_vec = p_vec,
                    d_target = d_target, tol = tol)
  result$root
}

# -----------------------------------------------------------------
# calculate_calendar_cutoffs — module entry point

calculate_calendar_cutoffs <- function(state) {
  design    <- state$design
  n_T       <- design$n_T
  n_C       <- design$n_C
  R         <- design$R
  L         <- design$L
  t_vec     <- design$t_vec
  p_vec     <- design$p_vec
  Lambda_T  <- design$Lambda_T
  Lambda_C  <- design$Lambda_C
  d_PFS_vec <- design$d_PFS_vec
  tol       <- state$options$tol

  A_vec  <- numeric(L)
  A_prev <- 0

  for (ell in seq_len(L)) {
    A_ell <- calculate_calendar_cutoff(
      n_T, n_C, Lambda_T, Lambda_C,
      R, t_vec, p_vec,
      d_target = d_PFS_vec[ell],
      A_prev   = A_prev,
      tol      = tol
    )
    A_vec[ell] <- A_ell
    A_prev     <- A_ell
  }

  state$design$A_vec <- A_vec
  state
}
