# per_subject_moments.R — Module 2: per-subject means (piecewise uniform)
#
# Piecewise-uniform accrual.  Functions implement the closed-form
# unconditional means from Blueprint §2.2 and §2.3.
#
# Dependencies: none (pure arithmetic)

# =================================================================
# 2.1  Segment integrals (i₁…i₅) and cumulative sums (I₁…I₅)
# =================================================================

i_1 <- function(tau1, tau2, A, mu) {
  tau2 - tau1
}

i_2 <- function(tau1, tau2, A, mu) {
  (tau2 - tau1) * (2 * A - tau1 - tau2) / 2
}

i_3 <- function(tau1, tau2, A, mu) {
  (tau2 - tau1) * (3 * A^2 - 3 * A * (tau1 + tau2) +
                   tau1^2 + tau1 * tau2 + tau2^2) / 3
}

i_4 <- function(tau1, tau2, A, mu) {
  (exp(-mu * (A - tau2)) - exp(-mu * (A - tau1))) / mu
}

i_5 <- function(tau1, tau2, A, mu) {
  (exp(-mu * (A - tau2)) * (A - tau2 + 1 / mu) -
   exp(-mu * (A - tau1)) * (A - tau1 + 1 / mu)) / mu
}

find_segment <- function(B, t_vec) {
  b <- which(B <= t_vec)[1]
  if (is.na(b)) length(t_vec) else b
}

I <- function(i_fn, A, B, b, mu, t_vec, p_vec) {
  total <- 0
  tau_start <- 0
  for (i in seq_len(b - 1)) {
    total     <- total + p_vec[i] * i_fn(tau_start, t_vec[i], A, mu)
    tau_start <- t_vec[i]
  }
  total + p_vec[b] * i_fn(tau_start, B, A, mu)
}

I_1 <- function(B, b, t_vec, p_vec)           I(i_1, 0, B, b, 0, t_vec, p_vec)
I_2 <- function(A, B, b, t_vec, p_vec)        I(i_2, A, B, b, 0, t_vec, p_vec)
I_3 <- function(A, B, b, t_vec, p_vec)        I(i_3, A, B, b, 0, t_vec, p_vec)
I_4 <- function(A, B, b, mu, t_vec, p_vec)    I(i_4, A, B, b, mu, t_vec, p_vec)
I_5 <- function(A, B, b, mu, t_vec, p_vec)    I(i_5, A, B, b, mu, t_vec, p_vec)

# =================================================================
# 2.2  Single-endpoint mean
# =================================================================

compute_mean <- function(variable, A, B, b, mu, t_vec, p_vec) {
  if (variable == "delta") {
    I_1(B, b, t_vec, p_vec) - I_4(A, B, b, mu, t_vec, p_vec)
  } else if (variable == "time") {
    (I_1(B, b, t_vec, p_vec) - I_4(A, B, b, mu, t_vec, p_vec)) / mu
  }

}

# =================================================================
# 2.2  Single-look product mean
# =================================================================

compute_single_product_mean <- function(variable1, variable2, A, B, b, mu, t_vec, p_vec) {
  if (variable1 == "delta" && variable2 == "delta") {
    I_1(B, b, t_vec, p_vec) - I_4(A, B, b, mu, t_vec, p_vec)
  } else if (variable1 == "time" && variable2 == "time") {
    2 / mu^2 * (I_1(B, b, t_vec, p_vec) - I_4(A, B, b, mu, t_vec, p_vec) -
            mu * I_5(A, B, b, mu, t_vec, p_vec))
  } else if ((variable1 == "delta" && variable2 == "time") ||
           (variable1 == "time" && variable2 == "delta")) {
    (I_1(B, b, t_vec, p_vec) - I_4(A, B, b, mu, t_vec, p_vec)) / mu -
      I_5(A, B, b, mu, t_vec, p_vec)
  }
}

# =================================================================
# 2.3  Cross-look product mean
# =================================================================

compute_cross_product_mean <- function(variable1, variable2, A1, A2, B, b,
                                      A_star, Lambda, lambda_1, lambda_2,
                                       t_vec, p_vec) {
  if (A1 > A2) {
    if (variable1 == "delta" && variable2 == "delta") {
      I_1(B, b, t_vec, p_vec) -
      I_4(A2, B, b, lambda_2, t_vec, p_vec)
    } else if (variable1 == "delta" && variable2 == "time") {
      (I_1(B, b, t_vec, p_vec) -
       I_4(A2, B, b, lambda_2, t_vec, p_vec)) / lambda_2 -
      exp(-Lambda * (A1 - A2)) *
      I_5(A2, B, b, Lambda, t_vec, p_vec)
    } else if (variable1 == "time" && variable2 == "delta") {
      (I_1(B, b, t_vec, p_vec) -
       I_4(A2, B, b, Lambda, t_vec, p_vec)) / Lambda -
      (I_4(A2, B, b, lambda_2, t_vec, p_vec) -
       I_4(A2, B, b, Lambda, t_vec, p_vec)) / lambda_1
    } else if (variable1 == "time" && variable2 == "time") {
      (I_1(B, b, t_vec, p_vec) -
       I_4(A2, B, b, Lambda, t_vec, p_vec)) / Lambda^2 -
      exp(-Lambda * (A1 - A2)) *
      I_5(A2, B, b, Lambda, t_vec, p_vec) / Lambda +
      (I_1(B, b, t_vec, p_vec) -
       I_4(A2, B, b, Lambda, t_vec, p_vec)) / (lambda_2 * Lambda) -
      (I_4(A2, B, b, lambda_2, t_vec, p_vec) -
       I_4(A2, B, b, Lambda, t_vec, p_vec)) / (lambda_1 * lambda_2)
    }
  } else {
    if (variable1 == "delta" && variable2 == "delta") {
      I_1(B, b, t_vec, p_vec) -
      I_4(A2, B, b, lambda_2, t_vec, p_vec) -
      I_4(A1, B, b, Lambda,   t_vec, p_vec) +
      I_4(A_star, B, b, Lambda, t_vec, p_vec)
    } else if (variable1 == "delta" && variable2 == "time") {
      (I_1(B, b, t_vec, p_vec) - I_4(A2, B, b, lambda_2, t_vec, p_vec)) /
        lambda_2 -
        I_5(A1, B, b, Lambda, t_vec, p_vec) -
        (I_4(A1, B, b, Lambda, t_vec, p_vec) -
         I_4(A_star, B, b, Lambda, t_vec, p_vec)) / lambda_2
    } else if (variable1 == "time" && variable2 == "delta") {
      (I_1(B, b, t_vec, p_vec) - I_4(A1, B, b, Lambda, t_vec, p_vec)) /
        Lambda -
        (I_4(A2, B, b, lambda_2, t_vec, p_vec) -
         I_4(A_star, B, b, Lambda, t_vec, p_vec)) / lambda_1
    } else if (variable1 == "time" && variable2 == "time") {
      I_1(B, b, t_vec, p_vec) / Lambda^2 -
      I_5(A1, B, b, Lambda, t_vec, p_vec) / Lambda -
      I_4(A1, B, b, Lambda, t_vec, p_vec) / Lambda^2 +
      (I_1(B, b, t_vec, p_vec) - I_4(A1, B, b, Lambda, t_vec, p_vec)) /
        (lambda_2 * Lambda) -
      I_4(A2, B, b, lambda_2, t_vec, p_vec) / (lambda_1 * lambda_2) +
      I_4(A_star, B, b, Lambda, t_vec, p_vec) / (lambda_1 * lambda_2)
    }
  }
}

# =================================================================
# 2.4  Array computation
# =================================================================

per_subject_moments <- function(state) {
  design <- state$design
  R     <- design$R
  t_vec <- design$t_vec
  p_vec <- design$p_vec
  L     <- design$L
  A_vec <- design$A_vec

  dn1 <- c("T", "C")
  dn2 <- c("delta", "time")
  dn3 <- c("PFS", "OS")
  single_look_means <- array(0, dim = c(2, 2, 2, L),
    dimnames = list(dn1, dn2, dn3, NULL))
  single_look_covariances <- array(0, dim = c(2, 2, 2, 2, L),
    dimnames = list(dn1, dn2, dn2, dn3, NULL))
  cross_look_covariances <- array(0, dim = c(2, 2, 2, L, L),
    dimnames = list(dn1, dn2, dn2, NULL, NULL))

  groups <- list(
    T = list(lam1 = design$lambda_1_T, lam2 = design$lambda_2_T,
             Lam  = design$Lambda_T),
    C = list(lam1 = design$lambda_1_C, lam2 = design$lambda_2_C,
             Lam  = design$Lambda_C)
  )

  vars <- c("delta", "time")

  for (g in c("T", "C")) {
    lam1 <- groups[[g]]$lam1
    lam2 <- groups[[g]]$lam2
    Lam  <- groups[[g]]$Lam

    for (ell in seq_len(L)) {
      A <- A_vec[ell]
      B <- min(A, R)
      b <- find_segment(B, t_vec)
      for (ep in c("PFS", "OS")) {
        mu <- switch(ep, PFS = Lam, OS = lam2)
        for (type in vars)
          single_look_means[g, type, ep, ell] <-
            compute_mean(type, A, B, b, mu, t_vec, p_vec)
        for (v1 in vars) for (v2 in vars)
          single_look_covariances[g, v1, v2, ep, ell] <-
            compute_single_product_mean(v1, v2, A, B, b, mu, t_vec, p_vec) -
            single_look_means[g, v1, ep, ell] *
            single_look_means[g, v2, ep, ell]
      }
    }

    if (L > 1)
      for (ell1 in seq_len(L))
        for (ell2 in seq_len(L)) {
          A1 <- A_vec[ell1]
          A2 <- A_vec[ell2]
          B  <- min(A1, A2, R)
          b  <- find_segment(B, t_vec)
          A_star <- (lam1 * A1 + lam2 * A2) / Lam
          for (v1 in vars) for (v2 in vars)
            cross_look_covariances[g, v1, v2, ell1, ell2] <-
              compute_cross_product_mean(v1, v2, A1, A2, B, b, A_star,
                                         Lam, lam1, lam2, t_vec, p_vec) -
              single_look_means[g, v1, "PFS", ell1] *
              single_look_means[g, v2, "OS",  ell2]
        }
  }

  state$moments <- list(
    single_look_means       = single_look_means,
    single_look_covariances = single_look_covariances,
    cross_look_covariances  = cross_look_covariances
  )

  d_OS_vec <- numeric(L)
  for (ell in seq_len(L))
    d_OS_vec[ell] <- design$n_T * single_look_means["T", "delta", "OS", ell] +
                     design$n_C * single_look_means["C", "delta", "OS", ell]
  state$design$d_OS_vec <- d_OS_vec

  state
}
