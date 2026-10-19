#' Matrix Utilities for Linear Models and Estimable Functions
#'
#' @description
#' Provides matrix operations fundamental to general linear models,
#' including the Goodnight SWEEP operator, reflexive generalized inverses (g2),
#' Moore-Penrose pseudo-inverses, rank evaluation, and estimability verification.
#'
#' @name utils_matrix
NULL

#' Goodnight Symmetric SWEEP Operator
#'
#' @description
#' Performs the symmetric SWEEP operation on column/row \code{k} of a symmetric
#' positive semi-definite matrix \code{A}.
#'
#' @param A A symmetric numeric matrix.
#' @param k Integer index of the pivot element.
#' @param tol Tolerance for pivot detection. If the pivot element is less than
#'   \code{tol}, the sweep operation is skipped.
#'
#' @return A numeric matrix of the same dimensions after sweeping column/row \code{k}.
#' @references
#' Goodnight, J. H. (1979). A Tutorial on the SWEEP Operator.
#' \emph{The American Statistician}, 33(3), 149-158.
#'
#' @export
sweep_matrix <- function(A, k, tol = 1e-7) {
  if (!is.matrix(A) || nrow(A) != ncol(A)) {
    stop("Input 'A' must be a square numeric matrix.")
  }
  p <- nrow(A)
  pivot <- A[k, k]
  if (is.na(pivot) || abs(pivot) < tol) {
    return(A)
  }

  B <- A
  inv_pivot <- 1 / pivot
  B[k, k] <- -inv_pivot

  # Update row and col k
  col_k <- A[, k]
  col_k[k] <- 0
  scaled_col <- col_k * inv_pivot

  B[, k] <- scaled_col
  B[k, ] <- scaled_col
  B[k, k] <- -inv_pivot

  # Rank-1 outer product update for remaining elements
  idx <- seq_len(p)
  idx_other <- idx[idx != k]

  if (length(idx_other) > 0) {
    outer_update <- tcrossprod(col_k[idx_other]) * inv_pivot
    B[idx_other, idx_other] <- A[idx_other, idx_other] - outer_update
  }

  B
}

#' Reflexive Generalized Inverse (Type g2)
#'
#' @description
#' Computes a reflexive generalized inverse \eqn{G} of a symmetric positive
#' semi-definite matrix \eqn{A} such that \eqn{A G A = A} and \eqn{G A G = G}.
#' Linearly dependent columns/rows are set to zero, matching the parameterization
#' conventions of classical general linear model procedures.
#'
#' @param A A symmetric positive semi-definite numeric matrix.
#' @param tol Numerical tolerance for pivot thresholding. Defaults to \code{1e-7}.
#'
#' @return A symmetric numeric matrix \eqn{G} representing the g2-inverse.
#'
#' @examples
#' X <- matrix(c(1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1), 4, 3, byrow = TRUE)
#' XtX <- crossprod(X)
#' G <- ginv_g2(XtX)
#' # Verifies A %*% G %*% A == A
#' max(abs(XtX %*% G %*% XtX - XtX)) < 1e-7
#'
#' @export
ginv_g2 <- function(A, tol = 1e-7) {
  if (!is.matrix(A) || nrow(A) != ncol(A)) {
    stop("Input 'A' must be a square symmetric matrix.")
  }
  p <- nrow(A)
  G <- A
  swept <- rep(FALSE, p)

  # Track scale of diagonals to apply relative tolerance
  diag_vals <- pmax(abs(diag(A)), 1)

  for (k in seq_len(p)) {
    pivot <- G[k, k]
    if (!is.na(pivot) && pivot > tol * diag_vals[k]) {
      G <- sweep_matrix(G, k, tol = tol * diag_vals[k])
      swept[k] <- TRUE
    } else {
      G[k, ] <- 0
      G[, k] <- 0
    }
  }

  # Negate swept elements to obtain generalized inverse
  G <- -G
  G[!swept, ] <- 0
  G[, !swept] <- 0
  dimnames(G) <- dimnames(A)
  G
}

#' Moore-Penrose Pseudo-Inverse via Singular Value Decomposition
#'
#' @description
#' Computes the unique Moore-Penrose pseudo-inverse \eqn{A^+} using singular
#' value decomposition (SVD).
#'
#' @param A A numeric matrix.
#' @param tol Tolerance for singular values. Defaults to \code{1e-7}.
#'
#' @return A numeric matrix representing the Moore-Penrose pseudo-inverse.
#' @export
ginv_svd <- function(A, tol = 1e-7) {
  if (!is.matrix(A)) {
    A <- as.matrix(A)
  }
  s <- svd(A)
  positive <- s$d > tol * max(s$d[1], 1)
  if (!any(positive)) {
    res <- matrix(0, ncol(A), nrow(A))
  } else {
    d_inv <- 1 / s$d[positive]
    u_pos <- s$u[, positive, drop = FALSE]
    v_pos <- s$v[, positive, drop = FALSE]
    res <- v_pos %*% (d_inv * t(u_pos))
  }
  dimnames(res) <- list(colnames(A), rownames(A))
  res
}

#' Numerical Matrix Rank
#'
#' @description
#' Determines the numerical rank of a matrix using QR decomposition or SVD.
#'
#' @param A A numeric matrix.
#' @param tol Tolerance for determining non-zero singular values or pivots.
#'
#' @importFrom Matrix rankMatrix
#' @export
matrix_rank <- function(A, tol = 1e-7) {
  if (!is.matrix(A) && !inherits(A, "Matrix")) {
    A <- as.matrix(A)
  }
  if (nrow(A) == 0 || ncol(A) == 0) {
    return(0L)
  }
  as.integer(Matrix::rankMatrix(A, tol = tol))
}

#' Estimability Verification for Linear Hypotheses
#'
#' @description
#' Tests whether each row of a hypothesis matrix \eqn{L} represents an
#' estimable linear function of the parameters in a linear model. A linear
#' function \eqn{L \beta} is estimable if and only if \eqn{L (X'X)^- (X'X) = L}.
#'
#' @param L A numeric vector or matrix where each row represents a linear combination.
#' @param XtX The cross-product matrix \eqn{X'X}.
#' @param tol Numerical tolerance for equality check. Defaults to \code{1e-7}.
#' @param ginv_fn Function used to compute the generalized inverse. Defaults to \code{ginv_g2}.
#'
#' @return A logical vector indicating whether each row of \code{L} is estimable.
#'
#' @references
#' Searle, S. R. (1971). \emph{Linear Models}. John Wiley & Sons.
#'
#' @examples
#' X <- matrix(c(1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1), 4, 3, byrow = TRUE)
#' XtX <- crossprod(X)
#' # Difference between groups is estimable
#' is_estimable(c(0, 1, -1), XtX)
#' # Single cell parameter in overparameterized model is not estimable alone
#' is_estimable(c(0, 1, 0), XtX)
#'
#' @export
is_estimable <- function(L, XtX, tol = 1e-7, ginv_fn = ginv_g2) {
  if (is.vector(L)) {
    L <- matrix(L, nrow = 1)
  }
  if (!is.matrix(L)) {
    stop("'L' must be a numeric vector or matrix.")
  }
  if (ncol(L) != ncol(XtX)) {
    stop(sprintf("Dimension mismatch: 'L' has %d columns, but 'XtX' has %d.",
                 ncol(L), ncol(XtX)))
  }

  G <- ginv_fn(XtX, tol = tol)
  # Project: L %*% G %*% XtX
  proj <- L %*% (G %*% XtX)

  # Check row-by-row max absolute difference relative to row scale
  diff_mat <- abs(proj - L)
  row_scales <- pmax(apply(abs(L), 1, max), 1)
  max_diffs <- apply(diff_mat, 1, max)

  max_diffs <= (tol * row_scales)
}

#' Safely Invert Quadratic Form Matrix
#'
#' @description
#' Computes the inverse or generalized inverse of the variance-covariance
#' matrix of a linear contrast \eqn{V = L G L'}.
#'
#' @param V A symmetric positive semi-definite numeric matrix.
#' @param tol Tolerance for pseudo-inversion. Defaults to \code{1e-7}.
#'
#' @return The inverted matrix.
#' @keywords internal
#' @noRd
quad_form_inv <- function(V, tol = 1e-7) {
  if (length(V) == 1) {
    if (abs(V) < tol) return(matrix(0, 1, 1))
    return(matrix(1 / V, 1, 1))
  }
  # Attempt Cholesky or standard solve; fall back to SVD pseudo-inverse if singular
  inv <- tryCatch({
    solve(V)
  }, error = function(e) {
    ginv_svd(V, tol = tol)
  })
  inv
}
