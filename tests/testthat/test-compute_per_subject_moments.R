# test-compute_per_subject_moments.R — per-subject moment tests.
#
# Run with \`devtools::test()\` from the package root.

test_that("moment arrays preserve structure and look schedule", {
  fixture <- local_moments_fixture()
  design <- fixture$design
  moments <- fixture$moments

  expect_equal(dim(moments$single_look_means), c(2, 2, 2, fixture$L))
  expect_equal(
    dim(moments$single_look_covariances),
    c(2, 2, 2, 2, fixture$L)
  )
  expect_equal(
    dim(moments$cross_look_covariances),
    c(2, 2, 2, fixture$L, fixture$L)
  )
  expect_equal(
    dimnames(moments$single_look_means)[1:3],
    list(fixture$groups, fixture$types, fixture$eps)
  )
  expect_equal(design$d_PFS_vec, c(30, 60))
  expect_true(all(diff(fixture$A_vec) > 0))
})

test_that("find_segment resolves boundaries and single-segment accrual", {
  expect_equal(
    vapply(
      c(3, 6, 9, 12, 15),
      find_segment,
      integer(1),
      t_vec = c(6, 12, 18)
    ),
    c(1L, 1L, 2L, 2L, 3L)
  )
  expect_equal(
    vapply(c(3, 10, 20), find_segment, integer(1), t_vec = c(10)),
    rep(1L, 3)
  )
})

test_that("basis integrals agree with closed-form reference values", {
  A <- 10
  mu <- 0.5

  expect_equal(c(i_1(0, 6, 0, 0), i_1(3, 8, 0, 0)), c(6, 5))
  expect_equal(i_2(0, 6, A, 0), 42)

  expected_i3 <- (6 - 0) *
    (3 * A^2 - 3 * A * (0 + 6) + 0^2 + 0 * 6 + 6^2) / 3
  expect_equal(i_3(0, 6, A, 0), expected_i3)

  expected_i4 <- (exp(-mu * (A - 6)) - exp(-mu * A)) / mu
  expect_equal(i_4(0, 6, A, mu), expected_i4)

  i5 <- i_5(0, 6, A, mu)
  expect_gt(i5, 0)
  expect_lt(i5, i_4(0, 6, A, mu) * 6)
})

test_that("cumulative integrals handle piecewise and uniform accrual", {
  t_vec <- c(6, 12)
  p_vec <- c(0.3, 0.7)
  A <- 10
  mu <- 0.5

  expect_equal(
    c(I_1(3, 1, t_vec, p_vec), I_1(9, 2, t_vec, p_vec)),
    c(0.9, 3.9)
  )
  expect_equal(I_2(A, 3, 1, t_vec, p_vec), 7.65)
  expect_gt(I_3(A, 3, 1, t_vec, p_vec), 0)

  expected_i4 <- p_vec[1] *
    (exp(-mu * (A - 3)) - exp(-mu * A)) / mu
  expect_equal(
    I_4(A, 3, 1, mu, t_vec, p_vec),
    expected_i4,
    tolerance = 1e-8
  )

  uniform <- (exp(-mu * (A - 3)) - exp(-mu * A)) / (mu * 6)
  expect_equal(
    I_4(A, 3, 1, mu, c(6), c(1 / 6)),
    uniform,
    tolerance = 1e-8
  )
  expect_gt(I_5(A, 3, 1, mu, t_vec, p_vec), 0)
  expect_equal(
    I_1(6, 1, t_vec, p_vec),
    I(i_1, 0, 6, 1, 0, t_vec, p_vec)
  )
})

pw_t <- c(6, 12)
pw_p <- c(1 / 12, 1 / 12)
A_val <- 10
mu_val <- 0.15

test_that("compute_mean supports event and time means across segments", {
  delta_first <- I_1(6, 1, pw_t, pw_p) -
    I_4(A_val, 6, 1, mu_val, pw_t, pw_p)
  delta_second <- I_1(9, 2, pw_t, pw_p) -
    I_4(A_val, 9, 2, mu_val, pw_t, pw_p)

  expect_equal(
    compute_mean("delta", A_val, 6, 1, mu_val, pw_t, pw_p),
    delta_first,
    tolerance = 1e-8
  )
  expect_equal(
    compute_mean("time", A_val, 6, 1, mu_val, pw_t, pw_p),
    delta_first / mu_val,
    tolerance = 1e-8
  )
  expect_equal(
    compute_mean("delta", A_val, 9, 2, mu_val, pw_t, pw_p),
    delta_second,
    tolerance = 1e-8
  )

  uniform_delta <- I_1(6, 1, c(12), c(1 / 12)) -
    I_4(A_val, 6, 1, mu_val, c(12), c(1 / 12))
  expect_equal(
    compute_mean("delta", A_val, 6, 1, mu_val, pw_t, pw_p),
    uniform_delta
  )
})

test_that("single-product means match all variable combinations", {
  md <- I_1(6, 1, pw_t, pw_p) -
    I_4(A_val, 6, 1, mu_val, pw_t, pw_p)
  i5 <- I_5(A_val, 6, 1, mu_val, pw_t, pw_p)

  expected <- c(
    delta_delta = md,
    time_time = 2 / mu_val^2 * (md - mu_val * i5),
    delta_time = md / mu_val - i5,
    time_delta = md / mu_val - i5
  )
  actual <- c(
    delta_delta = compute_single_product_mean(
      "delta", "delta", A_val, 6, 1, mu_val, pw_t, pw_p
    ),
    time_time = compute_single_product_mean(
      "time", "time", A_val, 6, 1, mu_val, pw_t, pw_p
    ),
    delta_time = compute_single_product_mean(
      "delta", "time", A_val, 6, 1, mu_val, pw_t, pw_p
    ),
    time_delta = compute_single_product_mean(
      "time", "delta", A_val, 6, 1, mu_val, pw_t, pw_p
    )
  )

  expect_equal(actual, expected, tolerance = 1e-8)
})

test_that("cross-product means match all blueprint formulas", {
  A1 <- 6
  A2 <- 12
  B <- 6
  b <- 1
  Lam <- 0.15
  l1 <- 0.1
  l2 <- 0.05
  A_star <- (l1 * A1 + l2 * A2) / Lam
  i1 <- I_1(B, b, pw_t, pw_p)
  i4_a1 <- I_4(A1, B, b, Lam, pw_t, pw_p)
  i4_a2 <- I_4(A2, B, b, l2, pw_t, pw_p)
  i4_star <- I_4(A_star, B, b, Lam, pw_t, pw_p)
  i5_a1 <- I_5(A1, B, b, Lam, pw_t, pw_p)

  expected <- c(
    delta_delta = i1 - i4_a2 - i4_a1 + i4_star,
    delta_time = (i1 - i4_a2) / l2 - i5_a1 -
      (i4_a1 - i4_star) / l2,
    time_delta = (i1 - i4_a1) / Lam -
      (i4_a2 - i4_star) / l1,
    time_time = i1 / Lam^2 - i5_a1 / Lam - i4_a1 / Lam^2 +
      (i1 - i4_a1) / (l2 * Lam) - i4_a2 / (l1 * l2) +
      i4_star / (l1 * l2)
  )
  actual <- c(
    delta_delta = compute_cross_product_mean(
      "delta", "delta", A1, A2, B, b, A_star,
      Lam, l1, l2, pw_t, pw_p
    ),
    delta_time = compute_cross_product_mean(
      "delta", "time", A1, A2, B, b, A_star,
      Lam, l1, l2, pw_t, pw_p
    ),
    time_delta = compute_cross_product_mean(
      "time", "delta", A1, A2, B, b, A_star,
      Lam, l1, l2, pw_t, pw_p
    ),
    time_time = compute_cross_product_mean(
      "time", "time", A1, A2, B, b, A_star,
      Lam, l1, l2, pw_t, pw_p
    )
  )

  expect_equal(actual, expected, tolerance = 1e-8)
})

test_that("computed moments satisfy mean and covariance invariants", {
  fixture <- local_moments_fixture()
  design <- fixture$design
  moments <- fixture$moments

  for (g in fixture$groups) {
    for (ep in fixture$eps) {
      expect_gt(
        moments$single_look_means[g, "delta", ep, 2],
        moments$single_look_means[g, "delta", ep, 1]
      )
      for (ell in seq_len(fixture$L)) {
        md <- moments$single_look_means[g, "delta", ep, ell]
        mtime <- moments$single_look_means[g, "time", ep, ell]
        expect_gte(md, 0)
        expect_lte(md, 1)
        expect_gt(mtime, 0)

        A <- fixture$A_vec[ell]
        B <- min(A, fixture$R)
        b <- find_segment(B, fixture$t_vec)
        mu <- if (g == "T") {
          if (ep == "PFS") design$Lambda_T else design$lambda_2_T
        } else {
          if (ep == "PFS") design$Lambda_C else design$lambda_2_C
        }
        mtime2 <- compute_single_product_mean(
          "time", "time", A, B, b, mu, fixture$t_vec, fixture$p_vec
        )
        expect_equal(
          moments$single_look_covariances[g, "time", "time", ep, ell],
          mtime2 - mtime^2,
          tolerance = 1e-8
        )
        expect_equal(
          moments$single_look_covariances[g, "delta", "delta", ep, ell],
          md * (1 - md),
          tolerance = 1e-8
        )
        expect_equal(
          moments$single_look_covariances[g, "delta", "time", ep, ell],
          moments$single_look_covariances[g, "time", "delta", ep, ell]
        )
      }
    }

    for (ell in seq_len(fixture$L)) {
      A <- fixture$A_vec[ell]
      B <- min(A, fixture$R)
      b <- find_segment(B, fixture$t_vec)
      mu <- if (g == "T") design$Lambda_T else design$Lambda_C
      md <- moments$single_look_means[g, "delta", "PFS", ell]
      mtime <- moments$single_look_means[g, "time", "PFS", ell]
      mv <- compute_single_product_mean(
        "delta", "time", A, B, b, mu, fixture$t_vec, fixture$p_vec
      )
      expect_equal(
        moments$single_look_covariances[g, "delta", "time", "PFS", ell],
        mv - md * mtime,
        tolerance = 1e-8
      )
    }
  }
})

test_that("cross-look covariances use cross moments and means", {
  fixture <- local_moments_fixture()
  design <- fixture$design
  moments <- fixture$moments
  A1 <- fixture$A_vec[1]
  A2 <- fixture$A_vec[2]
  B <- min(A1, fixture$R)
  b <- find_segment(B, fixture$t_vec)

  for (g in fixture$groups) {
    lam1 <- if (g == "T") design$lambda_1_T else design$lambda_1_C
    lam2 <- if (g == "T") design$lambda_2_T else design$lambda_2_C
    Lam <- if (g == "T") design$Lambda_T else design$Lambda_C
    A_star <- (lam1 * A1 + lam2 * A2) / Lam

    for (v1 in fixture$types) {
      for (v2 in fixture$types) {
        cross_mean <- compute_cross_product_mean(
          v1, v2, A1, A2, B, b, A_star, Lam, lam1, lam2,
          fixture$t_vec, fixture$p_vec
        )
        mean_1 <- moments$single_look_means[g, v1, "PFS", 1]
        mean_2 <- moments$single_look_means[g, v2, "OS", 2]
        expect_equal(
          moments$cross_look_covariances[g, v1, v2, 1, 2],
          cross_mean - mean_1 * mean_2,
          tolerance = 1e-8
        )
      }
    }
  }
})
