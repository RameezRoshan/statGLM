test_that("SWEEP operator satisfies mathematical properties", {
  A <- matrix(c(4, 2, 2, 8), 2, 2)
  # Sweeping all columns yields the negative inverse
  A_sw <- sweep_matrix(A, 1)
  A_sw <- sweep_matrix(A_sw, 2)
  expect_equal(-A_sw, solve(A), tolerance = 1e-7)
})

test_that("ginv_g2 satisfies reflexive generalized inverse conditions", {
  X <- matrix(c(1, 1, 0,
                1, 1, 0,
                1, 0, 1,
                1, 0, 1), 4, 3, byrow = TRUE)
  XtX <- crossprod(X)
  G <- ginv_g2(XtX)

  # Condition 1: A G A == A
  expect_equal(XtX %*% G %*% XtX, XtX, tolerance = 1e-7)

  # Condition 2: G A G == G (reflexive)
  expect_equal(G %*% XtX %*% G, G, tolerance = 1e-7)

  # Redundant parameter (column 3) is set to zero
  expect_equal(G[3, ], rep(0, 3))
  expect_equal(G[, 3], rep(0, 3))
})

test_that("ginv_svd satisfies all four Moore-Penrose conditions", {
  X <- matrix(c(1, 2, 3,
                2, 4, 6,
                1, 0, 1), 3, 3, byrow = TRUE)
  G <- ginv_svd(X)

  # 1. A A+ A == A
  expect_equal(unname(X %*% G %*% X), unname(X), tolerance = 1e-7)
  # 2. A+ A A+ == A+
  expect_equal(unname(G %*% X %*% G), unname(G), tolerance = 1e-7)
  # 3. (A A+)' == A A+
  expect_equal(unname(t(X %*% G)), unname(X %*% G), tolerance = 1e-7)
  # 4. (A+ A)' == A+ A
  expect_equal(unname(t(G %*% X)), unname(G %*% X), tolerance = 1e-7)
})

test_that("is_estimable correctly identifies estimable and inestimable functions", {
  X <- matrix(c(1, 1, 0,
                1, 1, 0,
                1, 0, 1,
                1, 0, 1), 4, 3, byrow = TRUE)
  XtX <- crossprod(X)

  # Group difference (0, 1, -1) is estimable
  expect_true(is_estimable(c(0, 1, -1), XtX))

  # Group 1 mean (1, 1, 0) is estimable
  expect_true(is_estimable(c(1, 1, 0), XtX))

  # Group 2 mean (1, 0, 1) is estimable
  expect_true(is_estimable(c(1, 0, 1), XtX))

  # Individual cell effect without intercept (0, 1, 0) is NOT estimable in overparameterized model
  expect_false(is_estimable(c(0, 1, 0), XtX))
})

test_that("quadratic form L G L' is invariant to generalized inverse choice", {
  X <- matrix(c(1, 1, 0,
                1, 1, 0,
                1, 0, 1,
                1, 0, 1), 4, 3, byrow = TRUE)
  XtX <- crossprod(X)
  G2 <- ginv_g2(XtX)
  G_svd <- ginv_svd(XtX)

  L <- matrix(c(0, 1, -1), 1, 3)
  val_g2 <- as.numeric(L %*% G2 %*% t(L))
  val_svd <- as.numeric(L %*% G_svd %*% t(L))

  expect_equal(val_g2, val_svd, tolerance = 1e-10)
  expect_equal(val_g2, 1.0, tolerance = 1e-10)
})
