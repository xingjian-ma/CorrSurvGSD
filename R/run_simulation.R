# run_simulation.R — Module 5: Monte Carlo simulation
#
# Simulates the Fleischer model under piecewise-uniform accrual with
# multiple interim looks. The helper functions use explicit inputs.
# run_simulation() is the Module 5 state handler.
#
# Dependencies: completed outputs of Modules 0 to 4

# -----------------------------------------------------------------
# simulate_trial_statistics — run Monte Carlo replicates
#
# Args:
#   design      — completed trial design list
#   n_sim       — number of valid simulation replicates
#   seed        — random seed for reproducibility
#   max_attempts — maximum total replicate-generation attempts
#
# Returns an n_sim by 2L matrix of valid Z-score statistics.

simulate_trial_statistics <- function(design, n_sim, seed,
                                      max_attempts = 10 * n_sim) {
  set.seed(seed)

  n_T   <- design$n_T
  n_C   <- design$n_C
  n     <- design$n
  L     <- design$L
  R_vec <- design$R_vec
  t_vec <- design$t_vec
  v_vec <- design$v_vec
  s     <- length(v_vec)

  column_names <- c(
    paste0("PFS_", seq_len(L)),
    paste0("OS_", seq_len(L))
  )
  z_mat <- matrix(
    NA_real_,
    nrow = n_sim,
    ncol = 2 * L,
    dimnames = list(NULL, column_names)
  )

  group <- c(rep("T", n_T), rep("C", n_C))
  endpoints <- c("PFS", "OS")
  d_PFS <- design$d_PFS_vec
  n_valid <- 0
  attempts <- 0

  while (n_valid < n_sim && attempts < max_attempts) {
    attempts <- attempts + 1

    # --- enrollment times: piecewise uniform -------------------------
    seg_probs <- v_vec * R_vec / n
    seg <- sample(seq_len(s), n, replace = TRUE, prob = seg_probs)
    U <- numeric(n)
    for (i in seq_len(s)) {
      idx <- which(seg == i)
      if (length(idx) == 0) {
        next
      }
      t_start <- if (i == 1) 0 else t_vec[i - 1]
      t_end   <- t_vec[i]
      U[idx] <- stats::runif(length(idx), t_start, t_end)
    }

    # --- latent failure times: Fleischer model -----------------------
    t1 <- numeric(n)
    t2 <- numeric(n)
    t1[group == "T"] <- stats::rexp(n_T, rate = design$lambda_1_T)
    t1[group == "C"] <- stats::rexp(n_C, rate = design$lambda_1_C)
    t2[group == "T"] <- stats::rexp(n_T, rate = design$lambda_2_T)
    t2[group == "C"] <- stats::rexp(n_C, rate = design$lambda_2_C)

    pfs_time <- pmin(t1, t2)
    os_time  <- t2

    # --- event-driven calendar cutoffs -------------------------------
    cal_pfs <- U + pfs_time
    A_sim   <- sort(cal_pfs)[d_PFS]

    # --- per-look Z-scores -------------------------------------------
    z <- numeric(2 * L)
    valid <- TRUE

    for (ell in seq_len(L)) {
      follow_up <- pmax(0, A_sim[ell] - U)

      for (ep in endpoints) {
        if (ep == "PFS") {
          delta   <- as.numeric(pfs_time <= follow_up)
          T_tilde <- pmin(pfs_time, follow_up)
        } else {
          delta   <- as.numeric(os_time <= follow_up)
          T_tilde <- pmin(os_time, follow_up)
        }

        d_T <- sum(delta[group == "T"])
        t_T <- sum(T_tilde[group == "T"])
        d_C <- sum(delta[group == "C"])
        t_C <- sum(T_tilde[group == "C"])

        if (d_T == 0 || t_T == 0 || d_C == 0 || t_C == 0) {
          valid <- FALSE
          break
        }

        idx <- if (ep == "PFS") ell else L + ell
        theta_hat <- log(d_T / t_T) - log(d_C / t_C)
        z[idx] <- -theta_hat / sqrt(1 / d_T + 1 / d_C)
      }
      if (!valid) {
        break
      }
    }

    if (valid) {
      n_valid <- n_valid + 1
      z_mat[n_valid, ] <- z
    }
  }

  if (n_valid < n_sim) {
    stop(
      "Unable to generate n_sim valid replicates within max_attempts; ",
      "check the trial design."
    )
  }

  z_mat
}

# -----------------------------------------------------------------
# first_efficacy_crossing — first allowed efficacy crossing

first_efficacy_crossing <- function(z_vec, boundary, gate_start = 1) {
  for (ell in seq_along(z_vec)) {
    if (z_vec[ell] <= boundary[ell, "futility"]) {
      return(NA_integer_)
    }
    if (ell >= gate_start &&
        z_vec[ell] >= boundary[ell, "efficacy"]) {
      return(ell)
    }
  }

  NA_integer_
}

# -----------------------------------------------------------------
# summarize_simulation_power — marginal and gatekeeping joint power

summarize_simulation_power <- function(z, boundary, hierarchy_order) {
  L <- ncol(z) / 2
  n_sim <- nrow(z)

  endpoint_indices <- list(
    PFS = seq_len(L),
    OS = L + seq_len(L)
  )
  endpoint_boundary <- list(
    PFS = boundary[endpoint_indices$PFS, , drop = FALSE],
    OS = boundary[endpoint_indices$OS, , drop = FALSE]
  )
  endpoint_incremental <- list(PFS = numeric(L), OS = numeric(L))
  primary_endpoint <- hierarchy_order[["primary"]]
  secondary_endpoint <- hierarchy_order[["secondary"]]
  joint_counts <- matrix(0, nrow = L, ncol = L)

  for (sim in seq_len(n_sim)) {
    endpoint_z <- list(
      PFS = z[sim, endpoint_indices$PFS],
      OS = z[sim, endpoint_indices$OS]
    )
    endpoint_look <- list()
    for (endpoint in c("PFS", "OS")) {
      endpoint_look[[endpoint]] <- first_efficacy_crossing(
        z_vec = endpoint_z[[endpoint]],
        boundary = endpoint_boundary[[endpoint]]
      )
      if (!is.na(endpoint_look[[endpoint]])) {
        endpoint_incremental[[endpoint]][endpoint_look[[endpoint]]] <-
          endpoint_incremental[[endpoint]][endpoint_look[[endpoint]]] + 1
      }
    }

    primary_look <- endpoint_look[[primary_endpoint]]
    if (!is.na(primary_look)) {
      secondary_look <- first_efficacy_crossing(
        z_vec = endpoint_z[[secondary_endpoint]],
        boundary = endpoint_boundary[[secondary_endpoint]],
        gate_start = primary_look
      )
      if (!is.na(secondary_look)) {
        joint_counts[primary_look, secondary_look] <-
          joint_counts[primary_look, secondary_look] + 1
      }
    }
  }

  marginal_power <- rbind(
    cbind(
      endpoint_incremental$PFS / n_sim,
      cumsum(endpoint_incremental$PFS / n_sim)
    ),
    cbind(
      endpoint_incremental$OS / n_sim,
      cumsum(endpoint_incremental$OS / n_sim)
    )
  )
  colnames(marginal_power) <- c("incremental", "cumulative")
  rownames(marginal_power) <- c(
    paste0("PFS_", seq_len(L)),
    paste0("OS_", seq_len(L))
  )

  joint_power_matrix <- joint_counts / n_sim
  joint_power_matrix[lower.tri(joint_power_matrix)] <- NA_real_
  rownames(joint_power_matrix) <- paste0(primary_endpoint, "_", seq_len(L))
  colnames(joint_power_matrix) <-
    paste0(secondary_endpoint, "_", seq_len(L))

  list(
    marginal_power = marginal_power,
    joint_power_matrix = joint_power_matrix,
    joint_power = sum(joint_power_matrix, na.rm = TRUE)
  )
}

# -----------------------------------------------------------------
# run_simulation — Module 5 state handler

run_simulation <- function(state) {
  design <- state$design
  n_sim <- state$options$n_sim
  simulation_seed <- state$options$simulation_seed

  z <- simulate_trial_statistics(
    design = design,
    n_sim = n_sim,
    seed = simulation_seed
  )
  state$empirical_results <- c(
    list(
      joint_mean_vector = colMeans(z),
      joint_correlation_matrix = stats::cor(z)
    ),
    summarize_simulation_power(
      z = z,
      boundary = design$boundary,
      hierarchy_order = design$hierarchy_order
    )
  )
  state
}
