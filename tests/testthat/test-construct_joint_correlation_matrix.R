# test-construct_joint_correlation_matrix.R — correlation matrix tests.
#
# Run with `devtools::test()` from the package root.
#
# =================================================================
# 1. Dimensions and structure
# =================================================================

test_that("output is (2L)×(2L)", {
  local_correlation_fixture()
  expect_equal(dim(R_mat), c(2 * L, 2 * L))
})

test_that("dimnames are PFS_i and OS_i", {
  local_correlation_fixture()
  nm <- c(paste0("PFS_", seq_len(L)), paste0("OS_", seq_len(L)))
  expect_equal(rownames(R_mat), nm)
  expect_equal(colnames(R_mat), nm)
})

test_that("matrix is symmetric", {
  local_correlation_fixture()
  expect_true(isSymmetric(R_mat))
})

test_that("diagonal is all ones", {
  local_correlation_fixture()
  expect_equal(as.vector(diag(R_mat)), rep(1, 2 * L))
})

# =================================================================
# 2. d_OS_vec
# =================================================================

test_that("d_OS_vec is positive and increasing", {
  local_correlation_fixture()
  expect_true(all(d_OS > 0))
  expect_gt(d_OS[2], d_OS[1])
})

test_that("d_OS_ell <= n (cannot exceed sample size)", {
  local_correlation_fixture()
  expect_true(all(d_OS <= design$n))
})

# =================================================================
# 3. Within-endpoint PFS block
# =================================================================

test_that("PFS block matches canonical G-S √(d_min/d_max)", {
  local_correlation_fixture()
  for (ell1 in seq_len(L))
    for (ell2 in seq_len(L)) {
      expected <- sqrt(d_PFS[min(ell1, ell2)] / d_PFS[max(ell1, ell2)])
      expect_equal(R_mat[ell1, ell2], expected, tolerance = 1e-8)
    }
})

# =================================================================
# 4. Within-endpoint OS block
# =================================================================

test_that("OS block matches canonical G-S √(d_min/d_max)", {
  local_correlation_fixture()
  for (ell1 in seq_len(L))
    for (ell2 in seq_len(L)) {
      expected <- sqrt(d_OS[min(ell1, ell2)] / d_OS[max(ell1, ell2)])
      expect_equal(R_mat[L + ell1, L + ell2], expected, tolerance = 1e-7)
    }
})

# =================================================================
# 5. Cross-endpoint block
# =================================================================

test_that("cross-endpoint block is finite and symmetric through transpose", {
  local_correlation_fixture()
  cross <- R_mat[seq_len(L), L + seq_len(L), drop = FALSE]
  expect_true(all(is.finite(cross)))
  expect_equal(cross, t(R_mat[L + seq_len(L), seq_len(L), drop = FALSE]))
  expect_true(any(cross != 0))
})

# =================================================================
# 6. Cross-endpoint — manual delta-method verification
# =================================================================

test_that("cross-endpoint ρ(PFS₁, OS₂) matches manual delta method", {
  local_correlation_fixture()
  ell1 <- 1
  ell2 <- 2

  for (g in groups) {
    n_g <- if (g == "T") n_T else n_C

    # marginal var PFS₁
    m_d1  <- sm[g, "delta", "PFS", ell1]
    m_T1  <- sm[g, "time",  "PFS", ell1]
    s_dd1 <- sc[g, "delta", "delta", "PFS", ell1]
    s_dT1 <- sc[g, "delta", "time",  "PFS", ell1]
    s_TT1 <- sc[g, "time",  "time",  "PFS", ell1]
    var1  <- (s_dd1 / m_d1^2 - 2 * s_dT1 / (m_d1 * m_T1) +
              s_TT1 / m_T1^2) / n_g

    # marginal var OS₂
    m_d2  <- sm[g, "delta", "OS", ell2]
    m_T2  <- sm[g, "time",  "OS", ell2]
    s_dd2 <- sc[g, "delta", "delta", "OS", ell2]
    s_dT2 <- sc[g, "delta", "time",  "OS", ell2]
    s_TT2 <- sc[g, "time",  "time",  "OS", ell2]
    var2  <- (s_dd2 / m_d2^2 - 2 * s_dT2 / (m_d2 * m_T2) +
              s_TT2 / m_T2^2) / n_g

    # cross cov PFS₁ × OS₂
    m_dP  <- sm[g, "delta", "PFS", ell1]
    m_TP  <- sm[g, "time",  "PFS", ell1]
    m_dO  <- sm[g, "delta", "OS",  ell2]
    m_TO  <- sm[g, "time",  "OS",  ell2]
    s_dd  <- cc[g, "delta", "delta", ell1, ell2]
    s_dT  <- cc[g, "delta", "time",  ell1, ell2]
    s_Td  <- cc[g, "time",  "delta", ell1, ell2]
    s_TT  <- cc[g, "time",  "time",  ell1, ell2]
    cov_g <- (s_dd / (m_dP * m_dO) - s_dT / (m_dP * m_TO) -
              s_Td / (m_TP * m_dO) + s_TT / (m_TP * m_TO)) / n_g

    var1_total <- if (g == "T") var1 else var1_total + var1
    var2_total <- if (g == "T") var2 else var2_total + var2
    cov_total  <- if (g == "T") cov_g else cov_total + cov_g
  }

  expected <- cov_total / sqrt(var1_total * var2_total)

  expect_equal(R_mat[ell1, L + ell2], expected, tolerance = 1e-8)
})

# =================================================================
# 7. Positive definiteness
# =================================================================

test_that("correlation matrix is positive definite", {
  local_correlation_fixture()
  ev <- eigen(R_mat, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(ev > 0))
})

# =================================================================
# 8. All correlations in (0, 1]
# =================================================================

test_that("all off-diagonal entries in [0, 1]", {
  local_correlation_fixture()
  off <- R_mat[upper.tri(R_mat)]
  expect_true(all(off >= 0))
  expect_true(all(off <= 1))
})

# =================================================================
# 9. L=1 edge case
# =================================================================

test_that("L=1 returns 2×2 correlation matrix", {
  state1 <- make_test_state(
    n_T = 80,
    n_C = 120,
    d_PFS_vec = c(30),
    t_vec = c(6),
    v_vec = c(200 / 12, 200 / 12),
    median_PFS_C = log(2) / 0.25,
    median_OS_C = log(2) / 0.10,
    median_PFS_T = log(2) / 0.08,
    median_OS_T = log(2) / 0.03,
    alpha_spending_PFS = "OF",
    alpha_spending_OS = "Pocock"
  )
  R1 <- state1$theoretical_results$joint_correlation_matrix

  expect_equal(dim(R1), c(2, 2))
  expect_equal(as.vector(diag(R1)), c(1, 1))
  expect_true(isSymmetric(R1))

  # within-endpoint blocks degenerate to 1 (only one look each)
  expect_equal(R1[1, 1], 1)
  expect_equal(R1[2, 2], 1)

  expect_true(is.finite(R1[1, 2]))
  expect_equal(R1[1, 2], R1[2, 1])
})
