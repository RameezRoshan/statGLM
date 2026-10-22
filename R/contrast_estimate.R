#' General Linear Hypotheses: Contrasts and Estimates
#'
#' @description
#' Implements estimability-validated hypothesis testing for general linear
#' models. \code{contrast_glm} evaluates multi-degree-of-freedom hypotheses
#' \eqn{H_0: L \beta = c} via F tests, and \code{estimate_glm} computes point
#' estimates, standard errors, t-statistics, and confidence intervals for
#' linear functions of the parameters.
#'
#' @name contrast_estimate
NULL

#' Helper to Expand Contrast Specifications to Full Parameter Vectors
#'
#' @param L Numeric vector, matrix, or named list.
#' @param object A \code{"stat_glm"} object.
#' @return A numeric matrix of dimensions \eqn{k \times p}.
#' @keywords internal
#' @noRd
expand_contrast_matrix <- function(L, object) {
  p <- length(object$coefficients)
  col_names <- names(object$coefficients)
  assign_ind <- object$design$assign_ind
  term_labels <- object$design$term_labels

  if (is.matrix(L)) {
    if (ncol(L) != p) {
      stop(sprintf("Matrix 'L' has %d columns, but model has %d parameters.",
                   ncol(L), p))
    }
    return(L)
  }

  if (is.numeric(L) && !is.list(L)) {
    if (length(L) != p) {
      stop(sprintf("Vector 'L' has length %d, but model has %d parameters.",
                   length(L), p))
    }
    return(matrix(L, nrow = 1, dimnames = list(NULL, col_names)))
  }

  if (is.list(L)) {
    # If list of multiple contrast specifications
    if (!is.null(names(L)) && all(names(L) %in% c("Intercept", term_labels))) {
      # Single contrast specified as named list of effect vectors
      vec <- numeric(p)
      names(vec) <- col_names
      for (eff in names(L)) {
        eff_vals <- L[[eff]]
        if (eff == "Intercept") {
          col_idx <- which(assign_ind == 0L)
        } else {
          term_idx <- which(term_labels == eff)
          col_idx <- which(assign_ind == term_idx)
        }
        if (length(eff_vals) != length(col_idx)) {
          stop(sprintf("Effect '%s' has %d columns in model, but %d values were provided in L.",
                       eff, length(col_idx), length(eff_vals)))
        }
        vec[col_idx] <- eff_vals
      }
      return(matrix(vec, nrow = 1, dimnames = list(NULL, col_names)))
    } else {
      # List of multiple contrasts
      rows <- lapply(L, function(sub_l) expand_contrast_matrix(sub_l, object))
      res <- do.call(rbind, rows)
      rownames(res) <- names(L)
      return(res)
    }
  }

  stop("Invalid specification for contrast matrix 'L'.")
}

#' Multi-Degree-of-Freedom Linear Contrast Tests (F-Tests)
#'
#' @description
#' Tests the general linear hypothesis \eqn{H_0: L \beta = c} for arbitrary
#' (single or multiple degree of freedom) contrasts. Validates the estimability
#' of each hypothesis row before calculating the contrast sum of squares,
#' degrees of freedom, mean square, F-statistic, and p-value.
#'
#' @param object An object of class \code{"stat_glm"}.
#' @param L A numeric matrix, vector, or named list of coefficients defining
#'   the linear contrast(s).
#' @param name An optional character string naming the contrast.
#' @param rhs A numeric vector specifying the null hypothesis values \eqn{c}.
#'   Defaults to \code{0}.
#' @param tol Numerical tolerance for estimability verification. Defaults to \code{1e-7}.
#'
#' @return An object of class \code{"contrast_stat_glm"} containing:
#' \itemize{
#'   \item \code{name}: Contrast name.
#'   \item \code{L}: The expanded \eqn{L} matrix.
#'   \item \code{df}: Numerator degrees of freedom (\eqn{\mathrm{rank}(L)}).
#'   \item \code{ss}: Contrast sum of squares.
#'   \item \code{ms}: Contrast mean square.
#'   \item \code{f_value}: F-statistic.
#'   \item \code{p_value}: Significance p-value.
#'   \item \code{estimable}: Logical vector indicating estimability of each row.
#' }
#'
#' @examples
#' data(harv)
#' fit <- stat_glm(bw ~ bw0 + stock + culturetype + sex, data = harv)
#'
#' # Test difference between Stock 1 and Stock 2
#' contrast_glm(fit, list(stock = c(1, -1)), name = "Stock 1 vs 2")
#'
#' @export
contrast_glm <- function(object, L, name = NULL, rhs = 0, tol = 1e-7) {
  if (!inherits(object, "stat_glm")) {
    stop("'object' must be an object of class 'stat_glm'.")
  }

  L_mat <- expand_contrast_matrix(L, object)
  k <- nrow(L_mat)

  if (length(rhs) == 1 && k > 1) {
    rhs <- rep(rhs, k)
  }
  if (length(rhs) != k) {
    stop(sprintf("'rhs' length (%d) must match number of contrast rows (%d).",
                 length(rhs), k))
  }

  # Verify estimability
  ginv_fn <- if (object$ginv_method == "g2") ginv_g2 else ginv_svd
  est_flags <- is_estimable(L_mat, object$XtX, tol = tol, ginv_fn = ginv_fn)

  if (!all(est_flags)) {
    inest_idx <- which(!est_flags)
    stop(sprintf("Hypothesis is non-estimable: row(s) %s violate L = L(X'X)^-(X'X).",
                 paste(inest_idx, collapse = ", ")))
  }

  b <- object$coefficients
  G <- object$G
  mse <- object$mse
  dfe <- object$df.residual

  diff_vec <- as.vector(L_mat %*% b - rhs)
  V <- L_mat %*% G %*% t(L_mat)
  V_inv <- quad_form_inv(V, tol = tol)

  ss_h <- as.numeric(t(diff_vec) %*% V_inv %*% diff_vec)
  df_h <- matrix_rank(V, tol = tol)
  ms_h <- if (df_h > 0) ss_h / df_h else NA_real_
  f_val <- if (!is.na(ms_h) && !is.na(mse) && mse > 0) ms_h / mse else NA_real_
  p_val <- if (!is.na(f_val) && df_h > 0 && dfe > 0) {
    stats::pf(f_val, df_h, dfe, lower.tail = FALSE)
  } else {
    NA_real_
  }

  if (is.null(name)) {
    name <- "Linear Contrast"
  }

  res <- list(
    name = name,
    L = L_mat,
    df = df_h,
    ss = ss_h,
    ms = ms_h,
    f_value = f_val,
    p_value = p_val,
    estimable = est_flags,
    df.residual = dfe,
    mse = mse
  )
  class(res) <- "contrast_stat_glm"
  res
}

#' @export
print.contrast_stat_glm <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nContrast Test:", x$name, "\n\n")
  tab <- data.frame(
    DF = x$df,
    `Contrast SS` = x$ss,
    `Mean Square` = x$ms,
    `F Value` = x$f_value,
    `Pr > F` = x$p_value,
    check.names = FALSE,
    row.names = x$name
  )
  stats::printCoefmat(as.matrix(tab), digits = digits, cs.ind = 1L, tst.ind = 4L,
                      zap.ps = TRUE, P.values = TRUE, has.Pvalue = TRUE, na.print = "")
  invisible(x)
}

#' Linear Function Parameter Estimates (t-Tests and Confidence Intervals)
#'
#' @description
#' Calculates point estimates, standard errors, t-statistics, and confidence
#' intervals for estimable linear functions of parameters \eqn{L \beta}.
#'
#' @param object An object of class \code{"stat_glm"}.
#' @param L A numeric matrix, vector, or named list specifying linear combinations.
#' @param name An optional character string naming the estimate.
#' @param rhs Hypothesized value \eqn{c} under \eqn{H_0: L \beta = c}. Defaults to \code{0}.
#' @param alpha Significance level for confidence intervals. Defaults to \code{0.05}.
#' @param tol Numerical tolerance for estimability verification. Defaults to \code{1e-7}.
#'
#' @return An object of class \code{"estimate_stat_glm"} containing estimates,
#'   standard errors, t-values, p-values, and confidence bounds.
#'
#' @examples
#' data(harv)
#' fit <- stat_glm(bw ~ bw0 + stock + culturetype + sex, data = harv)
#'
#' # Estimate Stock 1 - Stock 2 difference
#' estimate_glm(fit, list(stock = c(1, -1)), name = "Stock 1 - Stock 2")
#'
#' @export
estimate_glm <- function(object, L, name = NULL, rhs = 0, alpha = 0.05, tol = 1e-7) {
  if (!inherits(object, "stat_glm")) {
    stop("'object' must be an object of class 'stat_glm'.")
  }

  L_mat <- expand_contrast_matrix(L, object)
  k <- nrow(L_mat)

  if (length(rhs) == 1 && k > 1) {
    rhs <- rep(rhs, k)
  }
  if (length(rhs) != k) {
    stop(sprintf("'rhs' length (%d) must match number of estimate rows (%d).",
                 length(rhs), k))
  }

  # Estimability check
  ginv_fn <- if (object$ginv_method == "g2") ginv_g2 else ginv_svd
  est_flags <- is_estimable(L_mat, object$XtX, tol = tol, ginv_fn = ginv_fn)

  b <- object$coefficients
  G <- object$G
  mse <- object$mse
  dfe <- object$df.residual

  est_vals <- as.vector(L_mat %*% b)
  var_vals <- diag(L_mat %*% G %*% t(L_mat)) * mse
  se_vals <- sqrt(pmax(0, var_vals))

  # Inestimable rows receive NA
  est_vals[!est_flags] <- NA_real_
  se_vals[!est_flags] <- NA_real_

  t_vals <- (est_vals - rhs) / se_vals
  p_vals <- 2 * stats::pt(abs(t_vals), df = dfe, lower.tail = FALSE)

  t_crit <- stats::qt(1 - alpha / 2, df = dfe)
  ci_lower <- est_vals - t_crit * se_vals
  ci_upper <- est_vals + t_crit * se_vals

  row_labels <- rownames(L_mat)
  if (is.null(row_labels)) {
    if (!is.null(name) && k == 1) {
      row_labels <- name
    } else {
      row_labels <- paste("Row", seq_len(k))
    }
  }

  conf_pct <- sprintf("%g%%", (1 - alpha) * 100)
  tab <- data.frame(
    Parameter = row_labels,
    Estimate = est_vals,
    `Std Error` = se_vals,
    `t Value` = t_vals,
    `Pr > |t|` = p_vals,
    `Lower CL` = ci_lower,
    `Upper CL` = ci_upper,
    check.names = FALSE,
    row.names = row_labels
  )
  names(tab)[6:7] <- c(paste("Lower", conf_pct), paste("Upper", conf_pct))

  res <- list(
    table = tab,
    L = L_mat,
    estimable = est_flags,
    alpha = alpha,
    df.residual = dfe
  )
  class(res) <- "estimate_stat_glm"
  res
}

#' @export
print.estimate_stat_glm <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nParameter Estimates and Tests:\n\n")
  tab <- x$table[, -1, drop = FALSE]
  print(tab, digits = digits, na.print = "Non-est")
  invisible(x)
}
