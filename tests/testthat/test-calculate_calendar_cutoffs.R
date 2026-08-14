# test-calculate_calendar_cutoffs.R — calendar cutoff tests.
#
# Run with \`devtools::test()\` from the package root.

test_that("pfs_mean has expected boundaries and monotonicity", {
  A_seq <- seq(0, 24, length.out = 10)
  values <- vapply(
    A_seq,
    pfs_mean,
    numeric(1),
    R = 12,
    t_vec = c(6, 12),
    p_vec = c(1 / 12, 1 / 12),
    Lambda = 0.15
  )

  expect_equal(values[1], 0)
  expect_true(all(diff(values) >= 0))
  expect_gt(
    pfs_mean(120, 12, c(6, 12), c(1 / 12, 1 / 12), 0.15),
    0.99
  )

  expected <- 1 - (1 - exp(-1.8)) / 1.8
  expect_equal(
    pfs_mean(12, 12, c(6, 12), c(1 / 12, 1 / 12), 0.15),
    expected,
    tolerance = 1e-6
  )
})

test_that("pfs_event_diff crosses zero monotonically", {
  values <- vapply(
    seq(0, 48, length.out = 20),
    pfs_event_diff,
    numeric(1),
    n_T = 100,
    n_C = 100,
    Lambda_T = 0.15,
    Lambda_C = 0.15,
    R = 12,
    t_vec = c(6, 12),
    p_vec = c(1 / 12, 1 / 12),
    d_target = 50
  )

  expect_lt(values[1], 0)
  expect_gt(values[length(values)], 0)
  expect_true(all(diff(values) > 0))
})

test_that("calendar cutoff root handles guards and bracket extension", {
  cases <- list(
    list(
      n_T = 100,
      n_C = 100,
      Lambda_T = 0.15,
      Lambda_C = 0.15,
      R = 12,
      t_vec = c(6, 12),
      p_vec = c(1 / 12, 1 / 12),
      d_target = 50,
      A_prev = 0,
      tol = 1e-8
    ),
    list(
      n_T = 100,
      n_C = 100,
      Lambda_T = 0.15,
      Lambda_C = 0.15,
      R = 12,
      t_vec = c(6, 12),
      p_vec = c(1 / 12, 1 / 12),
      d_target = 120,
      A_prev = 10,
      tol = 1e-8
    ),
    list(
      n_T = 200,
      n_C = 200,
      Lambda_T = 0.07,
      Lambda_C = 0.07,
      R = 6,
      t_vec = c(3, 6),
      p_vec = c(1 / 6, 1 / 6),
      d_target = 300,
      A_prev = 0,
      tol = 1e-8
    )
  )

  for (case in cases) {
    root <- do.call(calculate_calendar_cutoff, case)
    residual <- do.call(
      pfs_event_diff,
      c(case[names(case) != "A_prev" & names(case) != "tol"],
        list(A = root))
    )
    expect_equal(residual, 0, tolerance = 1e-6)
    expect_gt(root, case$A_prev)
  }
  expect_gt(
    do.call(calculate_calendar_cutoff, cases[[3]]),
    12
  )
})

test_that("calendar cutoff module returns ordered roots", {
  state <- validate_test_input(
    n_T = 150,
    n_C = 150,
    d_PFS_vec = c(50, 100, 150),
    t_vec = c(12),
    v_vec = c(300 / 24, 300 / 24),
    median_PFS_C = log(2) / 0.25,
    median_OS_C = log(2) / 0.10,
    median_PFS_T = log(2) / 0.15,
    median_OS_T = log(2) / 0.05,
    alpha_spending_PFS = "OF",
    alpha_spending_OS = "Pocock"
  )
  state <- calculate_calendar_cutoffs(state)
  design <- state$design

  expect_length(design$A_vec, 3)
  expect_true(all(diff(design$A_vec) >= 0))
  for (ell in seq_along(design$A_vec)) {
    residual <- pfs_event_diff(
      design$A_vec[ell],
      150,
      150,
      0.15,
      0.25,
      24,
      c(12, 24),
      c(1 / 24, 1 / 24),
      design$d_PFS_vec[ell]
    )
    expect_equal(residual, 0, tolerance = 1e-6)
  }
})
