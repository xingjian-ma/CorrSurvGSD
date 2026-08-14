# test-construct_joint_correlation_matrix.R — correlation matrix tests.
#
# Run with \`devtools::test()\` from the package root.

test_that("correlation matrix preserves structure and names", {
  fixture <- local_correlation_fixture()
  R_mat <- fixture$R_mat
  expected_names <- c(
    paste0("PFS_", seq_len(fixture$L)),
    paste0("OS_", seq_len(fixture$L))
  )

  expect_equal(dim(R_mat), c(2 * fixture$L, 2 * fixture$L))
  expect_equal(dimnames(R_mat), list(expected_names, expected_names))
  expect_true(isSymmetric(R_mat))
  expect_equal(unname(diag(R_mat)), rep(1, 2 * fixture$L))
})

test_that("event counts and within-endpoint blocks are valid", {
  fixture <- local_correlation_fixture()
  R_mat <- fixture$R_mat
  L <- fixture$L

  expect_true(all(fixture$d_OS > 0))
  expect_true(all(diff(fixture$d_OS) > 0))
  expect_true(all(fixture$d_OS <= fixture$design$n))

  for (ell1 in seq_len(L)) {
    for (ell2 in seq_len(L)) {
      pfs_expected <- sqrt(
        fixture$d_PFS[min(ell1, ell2)] /
          fixture$d_PFS[max(ell1, ell2)]
      )
      os_expected <- sqrt(
        fixture$d_OS[min(ell1, ell2)] /
          fixture$d_OS[max(ell1, ell2)]
      )
      expect_equal(R_mat[ell1, ell2], pfs_expected, tolerance = 1e-8)
      expect_equal(
        R_mat[L + ell1, L + ell2],
        os_expected,
        tolerance = 1e-7
      )
    }
  }
})

test_that("cross-endpoint correlations match the delta-method reference", {
  fixture <- local_correlation_fixture()
  R_mat <- fixture$R_mat
  L <- fixture$L
  cross <- R_mat[seq_len(L), L + seq_len(L), drop = FALSE]

  expect_true(all(is.finite(cross)))
  expect_equal(
    cross,
    t(R_mat[L + seq_len(L), seq_len(L), drop = FALSE])
  )
  expect_true(any(cross != 0))

  ell1 <- 1
  ell2 <- 2
  var1_total <- 0
  var2_total <- 0
  cov_total <- 0
  for (g in fixture$groups) {
    n_g <- if (g == "T") fixture$n_T else fixture$n_C
    m_d1 <- fixture$sm[g, "delta", "PFS", ell1]
    m_T1 <- fixture$sm[g, "time", "PFS", ell1]
    s_dd1 <- fixture$sc[g, "delta", "delta", "PFS", ell1]
    s_dT1 <- fixture$sc[g, "delta", "time", "PFS", ell1]
    s_TT1 <- fixture$sc[g, "time", "time", "PFS", ell1]
    var1 <- (
      s_dd1 / m_d1^2 -
        2 * s_dT1 / (m_d1 * m_T1) +
        s_TT1 / m_T1^2
    ) / n_g

    m_d2 <- fixture$sm[g, "delta", "OS", ell2]
    m_T2 <- fixture$sm[g, "time", "OS", ell2]
    s_dd2 <- fixture$sc[g, "delta", "delta", "OS", ell2]
    s_dT2 <- fixture$sc[g, "delta", "time", "OS", ell2]
    s_TT2 <- fixture$sc[g, "time", "time", "OS", ell2]
    var2 <- (
      s_dd2 / m_d2^2 -
        2 * s_dT2 / (m_d2 * m_T2) +
        s_TT2 / m_T2^2
    ) / n_g

    m_dP <- fixture$sm[g, "delta", "PFS", ell1]
    m_TP <- fixture$sm[g, "time", "PFS", ell1]
    m_dO <- fixture$sm[g, "delta", "OS", ell2]
    m_TO <- fixture$sm[g, "time", "OS", ell2]
    s_dd <- fixture$cc[g, "delta", "delta", ell1, ell2]
    s_dT <- fixture$cc[g, "delta", "time", ell1, ell2]
    s_Td <- fixture$cc[g, "time", "delta", ell1, ell2]
    s_TT <- fixture$cc[g, "time", "time", ell1, ell2]
    cov_g <- (
      s_dd / (m_dP * m_dO) -
        s_dT / (m_dP * m_TO) -
        s_Td / (m_TP * m_dO) +
        s_TT / (m_TP * m_TO)
    ) / n_g

    var1_total <- var1_total + var1
    var2_total <- var2_total + var2
    cov_total <- cov_total + cov_g
  }

  expect_equal(
    R_mat[ell1, L + ell2],
    cov_total / sqrt(var1_total * var2_total),
    tolerance = 1e-8
  )
})

test_that("correlation matrix remains valid at ordinary and single-look designs", {
  fixture <- local_correlation_fixture()
  eigenvalues <- eigen(fixture$R_mat, symmetric = TRUE, only.values = TRUE)$values
  off_diagonal <- fixture$R_mat[upper.tri(fixture$R_mat)]

  expect_true(all(eigenvalues > 0))
  expect_true(all(off_diagonal >= 0))
  expect_true(all(off_diagonal <= 1))

  state <- make_test_state(
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
  single_look <- state$theoretical_results$joint_correlation_matrix

  expect_equal(dim(single_look), c(2, 2))
  expect_equal(unname(diag(single_look)), c(1, 1))
  expect_true(isSymmetric(single_look))
  expect_true(is.finite(single_look[1, 2]))
  expect_equal(single_look[1, 2], single_look[2, 1])
})
