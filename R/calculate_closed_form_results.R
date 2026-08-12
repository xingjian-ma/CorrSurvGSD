# calculate_closed_form_results.R — Module 4: power calculations
#
# The helper functions in this file use explicit numerical inputs.
# calculate_closed_form_results() is the Module 4 state handler and writes the
# resulting boundaries, mean vector, marginal power, and joint power.
#
# Dependencies: gsDesign, mvtnorm

# -----------------------------------------------------------------
# mean_drift — mean Z-score drift for one endpoint

mean_drift <- function(HR, r, d_vec) {
  -log(HR) * sqrt(d_vec * r / (1 + r)^2)
}

# -----------------------------------------------------------------
# build_boundary — futility and efficacy boundaries for one endpoint

build_boundary <- function(d_vec, L, alpha, alpha_spending,
                           alpha_spending_gamma, efficacy_looks,
                           futility_looks, futility_HR, r) {
  boundary <- matrix(
    c(-Inf, Inf),
    nrow = L,
    ncol = 2,
    byrow = TRUE,
    dimnames = list(NULL, c("futility", "efficacy"))
  )

  active_d_vec <- d_vec[efficacy_looks]
  if (length(efficacy_looks) == 1) {
    boundary[efficacy_looks, "efficacy"] <- stats::qnorm(1 - alpha)
  } else {
    gs_design_arguments <- list(
      k = length(efficacy_looks),
      alpha = alpha,
      timing = active_d_vec / d_vec[L],
      test.type = 1
    )
    if (identical(alpha_spending, "HSD")) {
      gs_design_arguments$sfu <- gsDesign::sfHSD
      gs_design_arguments$sfupar <- alpha_spending_gamma
    } else {
      gs_design_arguments$sfu <- alpha_spending
    }
    efficacy_design <- do.call(gsDesign::gsDesign, gs_design_arguments)
    boundary[efficacy_looks, "efficacy"] <- efficacy_design$upper$bound
  }

  if (length(futility_looks) > 0) {
    boundary[futility_looks, "futility"] <- mean_drift(
      HR = futility_HR,
      r = r,
      d_vec = d_vec[futility_looks]
    )
  }

  boundary
}

# -----------------------------------------------------------------
# convert_boundary_scales — convert finite Z boundaries to HR and p scales

convert_boundary_scales <- function(boundary, d_vec, r) {

  information_scale <- sqrt(d_vec * r / (1 + r)^2)
  boundary_HR <- matrix(
    NA_real_,
    nrow = nrow(boundary),
    ncol = ncol(boundary),
    dimnames = dimnames(boundary)
  )
  boundary_p <- matrix(
    NA_real_,
    nrow = nrow(boundary),
    ncol = ncol(boundary),
    dimnames = dimnames(boundary)
  )

  for (boundary_type in colnames(boundary)) {
    finite_boundary <- is.finite(boundary[, boundary_type])
    boundary_HR[finite_boundary, boundary_type] <- exp(
      -boundary[finite_boundary, boundary_type] /
        information_scale[finite_boundary]
    )
    boundary_p[finite_boundary, boundary_type] <- stats::pnorm(
      boundary[finite_boundary, boundary_type],
      lower.tail = FALSE
    )
  }

  list(HR = boundary_HR, p = boundary_p)
}

# -----------------------------------------------------------------
# first_crossing_bounds — bounds for first efficacy crossing at look

first_crossing_bounds <- function(boundary, look, gate_start = 1) {
  lower <- rep(-Inf, look)
  upper <- rep(Inf, look)

  if (look > 1) {
    previous_looks <- seq_len(look - 1)
    lower[previous_looks] <- boundary[previous_looks, "futility"]
    if (gate_start < look) {
      efficacy_indices <- seq.int(gate_start, look - 1)
      upper[efficacy_indices] <- boundary[efficacy_indices, "efficacy"]
    }
  }
  lower[look] <- boundary[look, "efficacy"]

  list(lower = lower, upper = upper)
}

# -----------------------------------------------------------------
# calculate_marginal_power — incremental and cumulative endpoint power

calculate_marginal_power <- function(mean_vector,
                                     correlation_matrix,
                                     boundary,
                                     L,
                                     seed = NULL) {
  power <- matrix(0, nrow = L, ncol = 2)
  colnames(power) <- c("incremental", "cumulative")

  for (ell in seq_len(L)) {
    if (!is.finite(boundary[ell, "efficacy"])) {
      next
    }
    bounds <- first_crossing_bounds(
      boundary = boundary,
      look = ell
    )
    indices <- seq_len(ell)
    if (ell == 1) {
      power[ell, "incremental"] <- stats::pnorm(
        bounds$upper,
        mean = mean_vector[indices]
      ) - stats::pnorm(
        bounds$lower,
        mean = mean_vector[indices]
      )
    } else {
      power[ell, "incremental"] <- as.numeric(mvtnorm::pmvnorm(
        lower = bounds$lower,
        upper = bounds$upper,
        mean = mean_vector[indices],
        corr = correlation_matrix[indices, indices, drop = FALSE],
        seed = seed
      ))
    }
  }

  power[, "cumulative"] <- cumsum(power[, "incremental"])

  power
}

# -----------------------------------------------------------------
# calculate_joint_power — generic primary-to-secondary joint power
#
# Returns an L by L upper-triangular matrix. Rows identify the primary
# endpoint efficacy look and columns identify the secondary endpoint
# efficacy look. Entries below the diagonal are NA.

calculate_joint_power <- function(joint_mean_vector,
                                  joint_correlation_matrix,
                                  boundary,
                                  hierarchy_order,
                                  L,
                                  seed = NULL) {
  endpoint_indices <- list(
    PFS = seq_len(L),
    OS = L + seq_len(L)
  )
  primary_endpoint <- hierarchy_order[["primary"]]
  secondary_endpoint <- hierarchy_order[["secondary"]]
  primary_indices <- endpoint_indices[[primary_endpoint]]
  secondary_indices <- endpoint_indices[[secondary_endpoint]]
  primary_boundary <- boundary[primary_indices, , drop = FALSE]
  secondary_boundary <- boundary[secondary_indices, , drop = FALSE]

  joint_power_matrix <- matrix(NA_real_, nrow = L, ncol = L)
  rownames(joint_power_matrix) <- paste0(primary_endpoint, "_", seq_len(L))
  colnames(joint_power_matrix) <- paste0(secondary_endpoint, "_", seq_len(L))

  for (primary_look in seq_len(L)) {
    for (secondary_look in seq.int(primary_look, L)) {
      if (!is.finite(primary_boundary[primary_look, "efficacy"]) ||
          !is.finite(secondary_boundary[secondary_look, "efficacy"])) {
        joint_power_matrix[primary_look, secondary_look] <- 0
        next
      }

      primary_bounds <- first_crossing_bounds(
        boundary = primary_boundary,
        look = primary_look
      )
      secondary_bounds <- first_crossing_bounds(
        boundary = secondary_boundary,
        look = secondary_look,
        gate_start = primary_look
      )
      joint_indices <- c(
        primary_indices[seq_len(primary_look)],
        secondary_indices[seq_len(secondary_look)]
      )
      lower <- c(primary_bounds$lower, secondary_bounds$lower)
      upper <- c(primary_bounds$upper, secondary_bounds$upper)

      joint_power_matrix[primary_look, secondary_look] <- as.numeric(
        mvtnorm::pmvnorm(
          lower = lower,
          upper = upper,
          mean = joint_mean_vector[joint_indices],
          corr = joint_correlation_matrix[joint_indices, joint_indices, drop = FALSE],
          seed = seed
        )
      )
    }
  }

  joint_power_matrix
}

# -----------------------------------------------------------------
# calculate_closed_form_results — Module 4 state handler
#
# Expects the completed outputs of Modules 1 to 3. It calculates only
# power-related quantities. Boundaries are appended to state$design,
# while power results are appended to state$theoretical_results.

calculate_closed_form_results <- function(state) {
  design <- state$design
  correlation_matrix <-
    state$theoretical_results$joint_correlation_matrix
  integration_seed <- state$options$integration_seed
  L <- design$L

  PFS_boundary <- build_boundary(
    d_vec = design$d_PFS_vec,
    L = L,
    alpha = design$alpha,
    alpha_spending = design$alpha_spending_PFS,
    alpha_spending_gamma = design$alpha_spending_gamma_PFS,
    efficacy_looks = design$efficacy_looks[["PFS"]],
    futility_looks = design$futility_looks[["PFS"]],
    futility_HR = design$futility_HR[["PFS"]],
    r = design$r
  )
  OS_boundary <- build_boundary(
    d_vec = design$d_OS_vec,
    L = L,
    alpha = design$alpha,
    alpha_spending = design$alpha_spending_OS,
    alpha_spending_gamma = design$alpha_spending_gamma_OS,
    efficacy_looks = design$efficacy_looks[["OS"]],
    futility_looks = design$futility_looks[["OS"]],
    futility_HR = design$futility_HR[["OS"]],
    r = design$r
  )
  PFS_boundary_scales <- convert_boundary_scales(
    boundary = PFS_boundary,
    d_vec = design$d_PFS_vec,
    r = design$r
  )
  OS_boundary_scales <- convert_boundary_scales(
    boundary = OS_boundary,
    d_vec = design$d_OS_vec,
    r = design$r
  )
  boundary_rownames <- c(
    paste0("PFS_", seq_len(L)),
    paste0("OS_", seq_len(L))
  )
  state$design$boundary <- rbind(PFS_boundary, OS_boundary)
  state$design$boundary_HR <- rbind(
    PFS_boundary_scales$HR,
    OS_boundary_scales$HR
  )
  state$design$boundary_p <- rbind(
    PFS_boundary_scales$p,
    OS_boundary_scales$p
  )
  rownames(state$design$boundary) <- boundary_rownames
  rownames(state$design$boundary_HR) <- boundary_rownames
  rownames(state$design$boundary_p) <- boundary_rownames
  boundary <- state$design$boundary

  state$theoretical_results$joint_mean_vector <- c(
    mean_drift(
      HR = design$HR_PFS,
      r = design$r,
      d_vec = design$d_PFS_vec
    ),
    mean_drift(
      HR = design$HR_OS,
      r = design$r,
      d_vec = design$d_OS_vec
    )
  )
  mean_vector <- state$theoretical_results$joint_mean_vector
  PFS_indices <- seq_len(L)
  OS_indices <- L + seq_len(L)
  PFS_marginal_power <- calculate_marginal_power(
    mean_vector = mean_vector[PFS_indices],
    correlation_matrix = correlation_matrix[PFS_indices, PFS_indices],
    boundary = boundary[PFS_indices, , drop = FALSE],
    L = L,
    seed = integration_seed
  )
  OS_marginal_power <- calculate_marginal_power(
    mean_vector = mean_vector[OS_indices],
    correlation_matrix = correlation_matrix[OS_indices, OS_indices],
    boundary = boundary[OS_indices, , drop = FALSE],
    L = L,
    seed = integration_seed
  )
  state$theoretical_results$marginal_power <- rbind(
    PFS_marginal_power,
    OS_marginal_power
  )
  rownames(state$theoretical_results$marginal_power) <- c(
    paste0("PFS_", PFS_indices),
    paste0("OS_", seq_len(L))
  )

  joint_power_matrix <- calculate_joint_power(
    joint_mean_vector = mean_vector,
    joint_correlation_matrix = correlation_matrix,
    boundary = boundary,
    hierarchy_order = design$hierarchy_order,
    L = L,
    seed = integration_seed
  )

  state$theoretical_results$joint_power_matrix <- joint_power_matrix
  state$theoretical_results$joint_power <-
    sum(joint_power_matrix, na.rm = TRUE)

  state
}
