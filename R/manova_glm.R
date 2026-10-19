#' Multivariate Analysis of Variance (MANOVA)
#'
#' @description
#' Conducts multivariate hypothesis testing for multivariate general linear models,
#' computing the four classical multivariate test statistics: Wilks' Lambda,
#' Pillai's Trace, Hotelling-Lawley Trace, and Roy's Greatest Root, along with
#' their exact and approximate F distributions and significance levels.
#'
#' @param object An object of class \code{"stat_mlm"} or \code{"stat_glm"} with a
#'   multivariate response.
#' @param type Numeric vector specifying which sums of squares and cross-products
#'   type to evaluate: \code{1}, \code{2}, or \code{3}. Defaults to \code{3}.
#' @param h Optional custom hypothesis SSCP matrix or contrast matrix.
#' @param e Optional custom error SSCP matrix.
#' @param tol Numerical tolerance for eigenvalues and singularity. Defaults to \code{1e-7}.
#' @param ... Additional arguments.
#'
#' @return An object of class \code{"manova_stat_glm"} containing:
#' \itemize{
#'   \item \code{results}: List of multivariate test tables for each evaluated effect.
#'   \item \code{E}: The error SSCP matrix.
#'   \item \code{df.residual}: Error degrees of freedom.
#'   \item \code{response_names}: Names of the response variables.
#' }
#'
#' @references
#' Rao, C. R. (1973). \emph{Linear Statistical Inference and Its Applications}. John Wiley & Sons.
#'
#' Pillai, K. C. S. (1960). \emph{Statistical Tables for Tests of Multivariate Hypotheses}.
#' The Statistical Center, University of the Philippines.
#'
#' @examples
#' data(iris)
#' fit <- stat_glm(cbind(Sepal.Length, Sepal.Width) ~ Species, data = iris)
#' manova_glm(fit)
#'
#' @export
manova_glm <- function(object, type = 3, h = NULL, e = NULL, tol = 1e-7, ...) {
  if (!inherits(object, "stat_mlm") && (!inherits(object, "stat_glm") || !is.matrix(object$design$y))) {
    stop("MANOVA requires a model fitted with a multivariate response matrix.")
  }

  type <- type[1]
  d <- object$design
  term_labels <- d$term_labels
  n_terms <- length(term_labels)
  p_vars <- ncol(d$y)
  ne <- object$df.residual

  # Determine Error SSCP matrix E
  if (is.null(e)) {
    E <- object$E
  } else {
    E <- as.matrix(e)
    if (nrow(E) != p_vars || ncol(E) != p_vars) {
      stop(sprintf("Error matrix 'e' dimensions (%d x %d) do not match response dimensions (%d x %d).",
                   nrow(E), ncol(E), p_vars, p_vars))
    }
  }

  E_inv <- tryCatch(solve(E), error = function(err) ginv_svd(E, tol = tol))

  # Helper to compute 4 test statistics given H, nh, and E, ne
  compute_mstat <- function(H, nh) {
    if (nh <= 0) {
      return(NULL)
    }

    # Generalized eigenvalues of E^-1 H
    HE <- E_inv %*% H
    eig <- eigen(HE, only.values = TRUE)$values
    # Filter real, positive eigenvalues
    real_eig <- Re(eig[abs(Im(eig)) < tol])
    pos_eig <- sort(real_eig[real_eig > tol], decreasing = TRUE)
    if (length(pos_eig) == 0) {
      pos_eig <- 0
    }

    s <- min(p_vars, nh)
    m <- (abs(p_vars - nh) - 1) / 2
    n_param <- (ne - p_vars - 1) / 2

    # 1. Wilks' Lambda
    wilks_val <- prod(1 / (1 + pos_eig))
    # Rao's F approximation
    r_val <- ne - (p_vars - nh + 1) / 2
    p2_nh2 <- p_vars^2 + nh^2
    u_val <- if (p2_nh2 == 5) 1 else sqrt(max(0, (p_vars^2 * nh^2 - 4) / (p2_nh2 - 5)))
    if (is.nan(u_val) || u_val == 0) u_val <- 1

    df1_w <- p_vars * nh
    df2_w <- r_val * u_val - (p_vars * nh - 2) / 2
    f_w <- if (wilks_val > 0 && wilks_val < 1) {
      ((1 - wilks_val^(1 / u_val)) / (wilks_val^(1 / u_val))) * (df2_w / df1_w)
    } else {
      0
    }
    p_w <- if (df2_w > 0 && f_w >= 0) stats::pf(f_w, df1_w, df2_w, lower.tail = FALSE) else NA_real_

    # 2. Pillai's Trace
    pillai_val <- sum(pos_eig / (1 + pos_eig))
    df1_p <- s * (2 * m + s + 1)
    df2_p <- s * (2 * n_param + s + 1)
    f_p <- if (s > pillai_val) (pillai_val / (s - pillai_val)) * (df2_p / df1_p) else NA_real_
    p_p <- if (!is.na(f_p) && df2_p > 0 && f_p >= 0) stats::pf(f_p, df1_p, df2_p, lower.tail = FALSE) else NA_real_

    # 3. Hotelling-Lawley Trace
    hl_val <- sum(pos_eig)
    df1_h <- s * (2 * m + s + 1)
    df2_h <- 2 * (s * n_param + 1)
    f_h <- (hl_val / s) * (df2_h / df1_h)
    p_h <- if (df2_h > 0 && f_h >= 0) stats::pf(f_h, df1_h, df2_h, lower.tail = FALSE) else NA_real_

    # 4. Roy's Greatest Root
    roy_val <- max(pos_eig)
    k_roy <- max(p_vars, nh)
    df1_r <- k_roy
    df2_r <- ne - k_roy + nh
    f_r <- roy_val * (df2_r / df1_r)
    p_r <- if (df2_r > 0 && f_r >= 0) stats::pf(f_r, df1_r, df2_r, lower.tail = FALSE) else NA_real_

    data.frame(
      Statistic = c("Wilks' Lambda", "Pillai's Trace", "Hotelling-Lawley Trace", "Roy's Greatest Root"),
      Value = c(wilks_val, pillai_val, hl_val, roy_val),
      `F Value` = c(f_w, f_p, f_h, f_r),
      `Num DF` = c(df1_w, df1_p, df1_h, df1_r),
      `Den DF` = c(df2_w, df2_p, df2_h, df2_r),
      `Pr > F` = c(p_w, p_p, p_h, p_r),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }

  results <- list()

  if (!is.null(h)) {
    # Custom hypothesis matrix
    H <- as.matrix(h)
    nh <- matrix_rank(H, tol = tol)
    results[["Custom Contrast"]] <- compute_mstat(H, nh)
  } else {
    # Evaluate for each effect
    # Type I, II, or III SSCP matrix H
    X_s <- d$X_sum
    XtX_s <- crossprod(X_s)
    XtX_inv <- tryCatch(solve(XtX_s), error = function(e) ginv_svd(XtX_s))
    B_sum <- XtX_inv %*% crossprod(X_s, d$y)

    for (k in seq_len(n_terms)) {
      t_lbl <- term_labels[k]
      cols_k <- which(d$assign_sum == k)
      if (length(cols_k) > 0) {
        B_k <- B_sum[cols_k, , drop = FALSE]
        V_k <- XtX_inv[cols_k, cols_k, drop = FALSE]
        V_k_inv <- quad_form_inv(V_k, tol = tol)
        # H = B_k' V_k^-1 B_k
        H_k <- crossprod(B_k, V_k_inv %*% B_k)
        nh_k <- matrix_rank(V_k, tol = tol)
        results[[t_lbl]] <- compute_mstat(H_k, nh_k)
      }
    }
  }

  res <- list(
    results = results,
    E = E,
    df.residual = ne,
    response_names = d$response_names,
    type = type
  )
  class(res) <- "manova_stat_glm"
  res
}

#' @export
print.manova_stat_glm <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nMultivariate Analysis of Variance (MANOVA)\n")
  cat("Dependent Variables:", paste(x$response_names, collapse = ", "), "\n")
  cat("Error Degrees of Freedom:", x$df.residual, "\n\n")

  for (eff in names(x$results)) {
    cat("Effect:", eff, "\n")
    tab <- x$results[[eff]]
    if (!is.null(tab)) {
      stats::printCoefmat(as.matrix(tab[, -1]), digits = digits, cs.ind = 1L, tst.ind = 2L,
                          zap.ps = TRUE, P.values = TRUE, has.Pvalue = TRUE, na.print = "")
    }
    cat("\n")
  }
  invisible(x)
}
