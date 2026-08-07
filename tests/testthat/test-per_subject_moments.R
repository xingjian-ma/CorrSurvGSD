# test_per_subject_moments.R — unit tests for per_subject_moments.R
#
# Run with `devtools::test()` from the package root.
#
# Setup: piecewise accrual (s=2), multi-look (L=2), asymmetric groups.

# =================================================================
# Shared setup — piecewise + multi-look + asymmetric
# =================================================================

state <- validate_trial_input(
  n_T = 80,
  n_C = 120,
  d_PFS_vec = c(30, 60),
  t_vec = c(6),
  v_vec = c(200 / 12, 200 / 12),
  median_PFS_C = log(2) / 0.25,
  median_OS_C = log(2) / 0.10,
  median_PFS_T = log(2) / 0.08,
  median_OS_T = log(2) / 0.03,
  alpha_spending_PFS = "OF",
  alpha_spending_OS = "Pocock"
)
state <- calendar_cutoff(state)
state <- per_subject_moments(state)
design <- state$design
moments <- state$moments

A_vec  <- design$A_vec
R      <- design$R
t_vec  <- design$t_vec
p_vec  <- design$p_vec
L      <- design$L
groups <- c("T", "C")
types  <- c("delta", "time")
eps    <- c("PFS", "OS")

# =================================================================
# 1. Array dimensions
# =================================================================

test_that("single_look_means dims", {
  expect_equal(dim(moments$single_look_means), c(2, 2, 2, L))
  expect_equal(dimnames(moments$single_look_means)[[1]], groups)
  expect_equal(dimnames(moments$single_look_means)[[2]], types)
  expect_equal(dimnames(moments$single_look_means)[[3]], eps)
})

test_that("single_look_covariances dims", {
  expect_equal(dim(moments$single_look_covariances), c(2, 2, 2, 2, L))
  expect_equal(dimnames(moments$single_look_covariances)[[1]], groups)
  expect_equal(dimnames(moments$single_look_covariances)[[2]], types)
  expect_equal(dimnames(moments$single_look_covariances)[[3]], types)
  expect_equal(dimnames(moments$single_look_covariances)[[4]], eps)
})

test_that("cross_look_covariances dims", {
  expect_equal(dim(moments$cross_look_covariances), c(2, 2, 2, L, L))
  expect_equal(dimnames(moments$cross_look_covariances)[[1]], groups)
  expect_equal(dimnames(moments$cross_look_covariances)[[2]], types)
  expect_equal(dimnames(moments$cross_look_covariances)[[3]], types)
})

# =================================================================
# 2. d_PFS_vec preserved in design
# =================================================================

test_that("d_PFS_vec matches input", {
  expect_equal(design$d_PFS_vec, c(30, 60), tolerance = 1e-6)
})

test_that("A_vec is strictly increasing", {
  expect_gt(A_vec[2], A_vec[1])
})

# =================================================================
# 3. find_segment
# =================================================================

test_that("B inside first segment returns b = 1", {
  expect_equal(find_segment(3, c(6, 12, 18)), 1)
})

test_that("B exactly on boundary returns that segment", {
  expect_equal(find_segment(6, c(6, 12, 18)), 1)
  expect_equal(find_segment(12, c(6, 12, 18)), 2)
})

test_that("B inside later segment returns correct b", {
  expect_equal(find_segment(9, c(6, 12, 18)), 2)
  expect_equal(find_segment(15, c(6, 12, 18)), 3)
})

test_that("single segment (s=1): any B returns b=1", {
  expect_equal(find_segment(3, c(10)), 1)
  expect_equal(find_segment(10, c(10)), 1)
  expect_equal(find_segment(20, c(10)), 1)
})

# =================================================================
# 4. Segment basis integrals i_1 … i_5
# =================================================================

test_that("i_1: integral of 1 du", {
  expect_equal(i_1(0, 6, 0, 0), 6)
  expect_equal(i_1(3, 8, 0, 0), 5)
})

test_that("i_2: integral of (A-u) du", {
  expect_equal(i_2(0, 6, 10, 0), 42)
})

test_that("i_3: integral of (A-u)^2 du", {
  expected <- (6 - 0) * (3 * 100 - 30 * (0 + 6) + 0 + 0 + 36) / 3
  expect_equal(i_3(0, 6, 10, 0), expected)
})

test_that("i_4: integral of exp(-mu(A-u)) du", {
  A <- 10
  mu <- 0.5
  expected <- (exp(-mu * (A - 6)) - exp(-mu * A)) / mu
  expect_equal(i_4(0, 6, A, mu), expected)
})

test_that("i_5: integral of (A-u)exp(-mu(A-u)) du is positive", {
  A <- 10
  mu <- 0.5
  val <- i_5(0, 6, A, mu)
  expect_gt(val, 0)
  expect_lt(val, i_4(0, 6, A, mu) * 6)
})

# =================================================================
# 5. Cumulative sums I_1 … I_5
# =================================================================

test_that("I_1: b=1, single segment", {
  expect_equal(I_1(3, 1, c(6, 12), c(0.3, 0.7)), 0.9)
})

test_that("I_1: b=2, full first + partial second", {
  expect_equal(I_1(9, 2, c(6, 12), c(0.3, 0.7)), 3.9)
})

test_that("I_2: b=1, single segment", {
  expect_equal(I_2(10, 3, 1, c(6, 12), c(0.3, 0.7)), 7.65)
})

test_that("I_3: b=1, single segment", {
  result <- I_3(10, 3, 1, c(6, 12), c(0.3, 0.7))
  expect_gt(result, 0)
})

test_that("I_4: b=1, single segment, matches manual expansion", {
  A <- 10
  mu <- 0.5
  B <- 3
  expected <- 0.3 * (exp(-mu * (A - B)) - exp(-mu * A)) / mu
  expect_equal(I_4(A, B, 1, mu, c(6, 12), c(0.3, 0.7)), expected, tolerance = 1e-8)
})

test_that("I_4: s=1 recovers single-segment uniform formula", {
  A <- 10
  mu <- 0.5
  B <- 3
  uniform <- (exp(-mu * (A - B)) - exp(-mu * A)) / (mu * 6)
  expect_equal(I_4(A, B, 1, mu, c(6), c(1 / 6)), uniform, tolerance = 1e-8)
})

test_that("I_5: returns positive value", {
  A <- 10
  mu <- 0.5
  val <- I_5(A, 3, 1, mu, c(6, 12), c(0.3, 0.7))
  expect_gt(val, 0)
})

test_that("I generic dispatcher matches I_1", {
  B <- 6
  b <- 1
  t <- c(6, 12)
  p <- c(0.3, 0.7)
  expect_equal(I_1(B, b, t, p), I(i_1, 0, B, b, 0, t, p))
})

# =================================================================
# 6. compute_mean — piecewise two-segment
# =================================================================

pw_t <- c(6, 12)
pw_p <- c(1 / 12, 1 / 12)
A_val <- 10
mu_val <- 0.15

test_that("compute_mean delta (b=1, B=6)", {
  expected <- I_1(6, 1, pw_t, pw_p) - I_4(A_val, 6, 1, mu_val, pw_t, pw_p)
  expect_equal(compute_mean("delta", A_val, 6, 1, mu_val, pw_t, pw_p),
               expected, tolerance = 1e-8)
})

test_that("compute_mean T (b=1, B=6)", {
  md <- I_1(6, 1, pw_t, pw_p) - I_4(A_val, 6, 1, mu_val, pw_t, pw_p)
  expect_equal(compute_mean("time", A_val, 6, 1, mu_val, pw_t, pw_p),
               md / mu_val, tolerance = 1e-8)
})

test_that("compute_mean delta (b=2, B=9)", {
  expected <- I_1(9, 2, pw_t, pw_p) - I_4(A_val, 9, 2, mu_val, pw_t, pw_p)
  expect_equal(compute_mean("delta", A_val, 9, 2, mu_val, pw_t, pw_p),
               expected, tolerance = 1e-8)
})

test_that("compute_mean with s=1 equals uniform case", {
  uni_delta <- I_1(6, 1, c(12), c(1 / 12)) - I_4(A_val, 6, 1, mu_val, c(12), c(1/12))
  pw_delta  <- compute_mean("delta", A_val, 6, 1, mu_val, pw_t, pw_p)
  expect_equal(pw_delta, uni_delta)
})

# =================================================================
# 7. compute_single_product_mean
# =================================================================

test_that("single_product_mean (delta,delta) = I_1 - I_4", {
  expected <- I_1(6, 1, pw_t, pw_p) - I_4(A_val, 6, 1, mu_val, pw_t, pw_p)
  expect_equal(compute_single_product_mean("delta", "delta", A_val, 6, 1,
                mu_val, pw_t, pw_p), expected, tolerance = 1e-8)
})

test_that("single_product_mean (T,T) = 2/mu^2 * (md - mu*I_5)", {
  md <- I_1(6, 1, pw_t, pw_p) - I_4(A_val, 6, 1, mu_val, pw_t, pw_p)
  i5 <- I_5(A_val, 6, 1, mu_val, pw_t, pw_p)
  expected <- 2 / mu_val^2 * (md - mu_val * i5)
  expect_equal(compute_single_product_mean("time", "time", A_val, 6, 1,
                mu_val, pw_t, pw_p), expected, tolerance = 1e-8)
})

test_that("single_product_mean (delta,T) = md/mu - I_5", {
  md <- I_1(6, 1, pw_t, pw_p) - I_4(A_val, 6, 1, mu_val, pw_t, pw_p)
  i5 <- I_5(A_val, 6, 1, mu_val, pw_t, pw_p)
  expected <- md / mu_val - i5
  expect_equal(compute_single_product_mean("delta", "time", A_val, 6, 1,
                mu_val, pw_t, pw_p), expected, tolerance = 1e-8)
})

test_that("single_product_mean (T,delta) equals (delta,T)", {
  expect_equal(compute_single_product_mean("time", "delta", A_val, 6, 1,
                mu_val, pw_t, pw_p),
               compute_single_product_mean("delta", "time", A_val, 6, 1,
                mu_val, pw_t, pw_p))
})

# =================================================================
# 8. compute_cross_product_mean
# =================================================================

A1 <- 6
A2 <- 12
B_c <- 6
b_c <- 1
Lam <- 0.15
l1 <- 0.1
l2 <- 0.05
A_star <- (l1 * A1 + l2 * A2) / Lam

test_that("cross_product_mean (delta,delta) matches blueprint", {
  expected <- I_1(B_c, b_c, pw_t, pw_p) -
              I_4(A2, B_c, b_c, l2, pw_t, pw_p) -
              I_4(A1, B_c, b_c, Lam, pw_t, pw_p) +
              I_4(A_star, B_c, b_c, Lam, pw_t, pw_p)
  expect_equal(compute_cross_product_mean("delta", "delta", A1, A2, B_c, b_c,
                A_star, Lam, l1, l2, pw_t, pw_p), expected, tolerance = 1e-8)
})

test_that("cross_product_mean (delta,T) matches blueprint", {
  i1 <- I_1(B_c, b_c, pw_t, pw_p)
  expected <- (i1 - I_4(A2, B_c, b_c, l2, pw_t, pw_p)) / l2 -
              I_5(A1, B_c, b_c, Lam, pw_t, pw_p) -
              (I_4(A1, B_c, b_c, Lam, pw_t, pw_p) -
               I_4(A_star, B_c, b_c, Lam, pw_t, pw_p)) / l2
  expect_equal(compute_cross_product_mean("delta", "time", A1, A2, B_c, b_c,
                A_star, Lam, l1, l2, pw_t, pw_p), expected, tolerance = 1e-8)
})

test_that("cross_product_mean (T,delta) matches blueprint", {
  i1 <- I_1(B_c, b_c, pw_t, pw_p)
  expected <- (i1 - I_4(A1, B_c, b_c, Lam, pw_t, pw_p)) / Lam -
              (I_4(A2, B_c, b_c, l2, pw_t, pw_p) -
               I_4(A_star, B_c, b_c, Lam, pw_t, pw_p)) / l1
  expect_equal(compute_cross_product_mean("time", "delta", A1, A2, B_c, b_c,
                A_star, Lam, l1, l2, pw_t, pw_p), expected, tolerance = 1e-8)
})

test_that("cross_product_mean (T,T) matches blueprint", {
  i1 <- I_1(B_c, b_c, pw_t, pw_p)
  i4_A1 <- I_4(A1, B_c, b_c, Lam, pw_t, pw_p)
  expected <- i1 / Lam^2 -
              I_5(A1, B_c, b_c, Lam, pw_t, pw_p) / Lam -
              i4_A1 / Lam^2 +
              (i1 - i4_A1) / (l2 * Lam) -
              I_4(A2, B_c, b_c, l2, pw_t, pw_p) / (l1 * l2) +
              I_4(A_star, B_c, b_c, Lam, pw_t, pw_p) / (l1 * l2)
  expect_equal(compute_cross_product_mean("time", "time", A1, A2, B_c, b_c,
                A_star, Lam, l1, l2, pw_t, pw_p), expected, tolerance = 1e-8)
})

# =================================================================
# 9. single_look_means — per-look sanity checks
# =================================================================

test_that("all m_delta in [0, 1]", {
  for (g in groups) for (ep in eps) for (ell in seq_len(L)) {
    val <- moments$single_look_means[g, "delta", ep, ell]
    expect_gte(val, 0)
    expect_lte(val, 1)
  }
})

test_that("all m_T are positive", {
  for (g in groups) for (ep in eps) for (ell in seq_len(L)) {
    val <- moments$single_look_means[g, "time", ep, ell]
    expect_gt(val, 0)
  }
})

test_that("m_delta increases across looks", {
  for (g in groups) for (ep in eps) {
    expect_gt(moments$single_look_means[g, "delta", ep, 2],
              moments$single_look_means[g, "delta", ep, 1])
  }
})

# =================================================================
# 10. single_look_covariances — formula verification
# =================================================================

test_that("sigma2_delta = m_delta * (1 - m_delta)", {
  for (g in groups) for (ep in eps) for (ell in seq_len(L)) {
    md  <- moments$single_look_means[g, "delta", ep, ell]
    sig <- moments$single_look_covariances[g, "delta", "delta", ep, ell]
    expect_equal(sig, md * (1 - md), tolerance = 1e-8)
  }
})

test_that("sigma2_T = m_T2 - m_T^2", {
  for (g in groups) for (ep in eps) for (ell in seq_len(L)) {
    A   <- A_vec[ell]
    B   <- min(A, R)
    b   <- find_segment(B, t_vec)
    mu  <- if (g == "T") {
      if (ep == "PFS") design$Lambda_T else design$lambda_2_T
    } else {
      if (ep == "PFS") design$Lambda_C else design$lambda_2_C
    }
    mT  <- moments$single_look_means[g, "time", ep, ell]
    mT2 <- compute_single_product_mean("time", "time", A, B, b, mu, t_vec, p_vec)
    expect_equal(moments$single_look_covariances[g, "time", "time", ep, ell],
                 mT2 - mT^2, tolerance = 1e-8)
  }
})

test_that("single_look_cov is symmetric", {
  for (g in groups) for (ep in eps) for (ell in seq_len(L)) {
    expect_equal(moments$single_look_covariances[g, "delta", "time", ep, ell],
                 moments$single_look_covariances[g, "time", "delta", ep, ell])
  }
})

test_that("cov(delta,T) = mv - m_delta * m_T", {
  for (g in groups) for (ell in seq_len(L)) {
    A  <- A_vec[ell]
    B  <- min(A, R)
    b  <- find_segment(B, t_vec)
    mu <- if (g == "T") design$Lambda_T else design$Lambda_C
    md <- moments$single_look_means[g, "delta", "PFS", ell]
    mT <- moments$single_look_means[g, "time", "PFS", ell]
    mv <- compute_single_product_mean("delta", "time", A, B, b, mu, t_vec, p_vec)
    expect_equal(moments$single_look_covariances[g, "delta", "time", "PFS", ell],
                 mv - md * mT, tolerance = 1e-8)
  }
})

# =================================================================
# 11. cross_look_covariances — structure and values
# =================================================================

test_that("cross_look cov = mv - m_PFS_ell1 * m_OS_ell2", {
  for (g in groups) {
    lam1 <- if (g == "T") design$lambda_1_T else design$lambda_1_C
    lam2 <- if (g == "T") design$lambda_2_T else design$lambda_2_C
    Lam  <- if (g == "T") design$Lambda_T   else design$Lambda_C
    A1   <- A_vec[1]
    A2   <- A_vec[2]
    B    <- min(A1, R)
    b    <- find_segment(B, t_vec)
    A_star <- (lam1 * A1 + lam2 * A2) / Lam
    for (v1 in types) for (v2 in types) {
      mv <- compute_cross_product_mean(v1, v2, A1, A2, B, b, A_star,
                                       Lam, lam1, lam2, t_vec, p_vec)
      m1 <- moments$single_look_means[g, v1, "PFS", 1]
      m2 <- moments$single_look_means[g, v2, "OS", 2]
      expect_equal(moments$cross_look_covariances[g, v1, v2, 1, 2],
                   mv - m1 * m2, tolerance = 1e-8)
    }
  }
})
