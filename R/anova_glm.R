#' Analysis of Variance (ANOVA) for General Linear Models
#'
#' @description
#' Computes Type I (sequential), Type II (hierarchical/containment), and
#' Type III (orthogonal/partial) sums of squares for linear models, along with
#' comprehensive effect size measures (partial and semipartial eta-squared,
#' omega-squared, noncentrality parameters, and noncentral confidence intervals).
#'
#' @param object An object of class \code{"stat_glm"}.
#' @param type Numeric vector specifying which sums of squares to compute.
#'   Any combination of \code{1}, \code{2}, and \code{3}. Defaults to \code{c(1, 3)}.
#' @param effectsize Logical; if \code{TRUE}, adds effect size estimates and
#'   confidence intervals to the ANOVA table. Defaults to \code{FALSE}.
#' @param alpha Significance level for effect size confidence intervals. Defaults to \code{0.05}
#'   (producing 100 * (1 - alpha) percent confidence limits, e.g., 90 percent when \code{alpha = 0.10}).
#' @param error_term Optional character string naming an effect to use as the denominator
#'   error term in hypothesis tests (equivalent to the SAS GLM \code{TEST} statement).
#' @param ... Additional arguments passed to methods.
#'
#' @return An object of class \code{"anova_stat_glm"} containing ANOVA tables for
#'   each requested sum of squares type.
#'
#' @references
#' Searle, S. R. (1971). \emph{Linear Models}. John Wiley & Sons.
#'
#' Searle, S. R. (1987). \emph{Linear Models for Unbalanced Data}. John Wiley & Sons.
#'
#' Goodnight, J. H. (1978). Tests of Hypotheses in Fixed Effects Linear Models.
#' \emph{Communications in Statistics - Theory and Methods}, 7(9), 793-801.
#'
#' @examples
#' a <- factor(c(1,1,1, 2,2,2, 3,3,3,3))
#' b <- factor(c(1,1,2, 1,2,2, 1,1,2,2))
#' y <- c(23.5, 23.7, 28.7, 8.9, 5.6, 8.9, 10.3, 12.5, 13.6, 14.6)
#' df <- data.frame(a, b, y)
#' fit <- stat_glm(y ~ a * b, data = df)
#'
#' # Type I, II, and III Sums of Squares
#' anova_glm(fit, type = 1)
#' anova_glm(fit, type = 2)
#' anova_glm(fit, type = 3, effectsize = TRUE)
#'
#' @export
anova_glm <- function(object, type = c(1, 3), effectsize = FALSE,
                      alpha = 0.05, error_term = NULL, ...) {
  if (!inherits(object, "stat_glm")) {
    stop("'object' must be an object of class 'stat_glm'.")
  }

  if (inherits(object, "stat_mlm")) {
    return(manova_glm(object, type = type, ...))
  }

  valid_types <- c(1, 2, 3)
  type <- intersect(type, valid_types)
  if (length(type) == 0) {
    stop("'type' must contain at least one of 1, 2, or 3.")
  }

  d <- object$design
  term_labels <- d$term_labels
  n_terms <- length(term_labels)
  X_ind <- d$X_ind
  X_sum <- d$X_sum
  y <- d$y
  weights <- object$weights
  n <- d$n_obs
  dfe <- object$df.residual
  mse <- object$mse
  sst <- object$sst
  sse <- object$sse

  w_sqrt <- sqrt(weights)
  y_w <- y * w_sqrt

  # Determine denominator MSE and DFE (standard or custom error term)
  custom_error <- FALSE
  err_ms <- mse
  err_df <- dfe
  err_name <- "Error"

  # Helper function to fit a weighted submodel on columns
  fit_submodel <- function(X_cols) {
    X_sub <- X_cols * w_sqrt
    rk_sub <- matrix_rank(X_sub)
    if (rk_sub == 0) {
      y_hat <- rep(0, n)
    } else {
      G_sub <- ginv_g2(crossprod(X_sub))
      b_sub <- G_sub %*% crossprod(X_sub, y_w)
      y_hat <- X_cols %*% b_sub
    }
    sse_sub <- sum(weights * (y - as.vector(y_hat))^2)
    list(sse = sse_sub, rank = rk_sub)
  }

  tables <- list()

  # ---------------------------------------------------------
  # 1. Type I Sums of Squares (Sequential)
  # ---------------------------------------------------------
  if (1 %in% type) {
    ss1_vals <- numeric(n_terms)
    df1_vals <- integer(n_terms)

    # Base model: Intercept only
    curr_cols <- if (d$has_intercept) which(d$assign_ind == 0L) else integer(0)
    if (length(curr_cols) > 0) {
      fit_prev <- fit_submodel(X_ind[, curr_cols, drop = FALSE])
    } else {
      fit_prev <- list(sse = sum(weights * y^2), rank = 0L)
    }

    for (k in seq_len(n_terms)) {
      term_cols <- which(d$assign_ind == k)
      curr_cols <- c(curr_cols, term_cols)
      fit_curr <- fit_submodel(X_ind[, curr_cols, drop = FALSE])

      ss1_vals[k] <- max(0, fit_prev$sse - fit_curr$sse)
      df1_vals[k] <- max(0L, fit_curr$rank - fit_prev$rank)
      fit_prev <- fit_curr
    }

    ms1_vals <- ifelse(df1_vals > 0, ss1_vals / df1_vals, NA_real_)
    tables[["Type I"]] <- data.frame(
      DF = df1_vals,
      `Sum of Squares` = ss1_vals,
      `Mean Square` = ms1_vals,
      row.names = term_labels,
      check.names = FALSE
    )
  }

  # ---------------------------------------------------------
  # 2. Type II Sums of Squares (Hierarchical / Containment)
  # ---------------------------------------------------------
  if (2 %in% type) {
    ss2_vals <- numeric(n_terms)
    df2_vals <- integer(n_terms)

    for (k in seq_len(n_terms)) {
      t_lbl <- term_labels[k]
      t_vars <- d$term_vars[[t_lbl]]

      # Identify terms containing t_lbl
      contained_in_others <- logical(n_terms)
      for (other_k in seq_len(n_terms)) {
        if (other_k != k) {
          other_vars <- d$term_vars[[term_labels[other_k]]]
          contained_in_others[other_k] <- is_contained(t_vars, other_vars, d$var_types)
        }
      }

      # Non-containing terms (excluding k)
      non_containing_idx <- which(!contained_in_others & (seq_len(n_terms) != k))

      # Model 0: Intercept + non-containing terms
      cols_0 <- if (d$has_intercept) which(d$assign_ind == 0L) else integer(0)
      for (nc in non_containing_idx) {
        cols_0 <- c(cols_0, which(d$assign_ind == nc))
      }

      if (length(cols_0) > 0) {
        fit_0 <- fit_submodel(X_ind[, cols_0, drop = FALSE])
      } else {
        fit_0 <- list(sse = sum(weights * y^2), rank = 0L)
      }

      # Model 1: Model 0 + term k
      cols_1 <- c(cols_0, which(d$assign_ind == k))
      fit_1 <- fit_submodel(X_ind[, cols_1, drop = FALSE])

      ss2_vals[k] <- max(0, fit_0$sse - fit_1$sse)
      df2_vals[k] <- max(0L, fit_1$rank - fit_0$rank)
    }

    ms2_vals <- ifelse(df2_vals > 0, ss2_vals / df2_vals, NA_real_)
    tables[["Type II"]] <- data.frame(
      DF = df2_vals,
      `Sum of Squares` = ss2_vals,
      `Mean Square` = ms2_vals,
      row.names = term_labels,
      check.names = FALSE
    )
  }

  # ---------------------------------------------------------
  # 3. Type III Sums of Squares (Orthogonal / Partial)
  # ---------------------------------------------------------
  if (3 %in% type) {
    ss3_vals <- numeric(n_terms)
    df3_vals <- integer(n_terms)

    # Full-rank sum-contrast quadratic form
    X_s_w <- X_sum * w_sqrt
    XtX_sum <- crossprod(X_s_w)
    XtX_inv <- tryCatch(solve(XtX_sum), error = function(e) ginv_svd(XtX_sum))
    b_sum <- as.vector(XtX_inv %*% crossprod(X_s_w, y_w))

    for (k in seq_len(n_terms)) {
      cols_k <- which(d$assign_sum == k)
      if (length(cols_k) > 0) {
        b_k <- b_sum[cols_k]
        V_k <- XtX_inv[cols_k, cols_k, drop = FALSE]
        V_k_inv <- quad_form_inv(V_k)
        ss_val <- as.numeric(crossprod(b_k, V_k_inv %*% b_k))
        df_val <- matrix_rank(V_k)

        ss3_vals[k] <- max(0, ss_val)
        df3_vals[k] <- df_val
      } else {
        ss3_vals[k] <- 0
        df3_vals[k] <- 0L
      }
    }

    ms3_vals <- ifelse(df3_vals > 0, ss3_vals / df3_vals, NA_real_)
    tables[["Type III"]] <- data.frame(
      DF = df3_vals,
      `Sum of Squares` = ss3_vals,
      `Mean Square` = ms3_vals,
      row.names = term_labels,
      check.names = FALSE
    )
  }

  # Handle custom error term (TEST statement)
  if (!is.null(error_term)) {
    if (!error_term %in% term_labels) {
      stop(sprintf("Error term '%s' not found among model terms: %s",
                   error_term, paste(term_labels, collapse = ", ")))
    }
    # Use specified error term from the highest type requested
    last_tab <- tables[[length(tables)]]
    err_ms <- last_tab[error_term, "Mean Square"]
    err_df <- last_tab[error_term, "DF"]
    err_name <- error_term
    custom_error <- TRUE
  }

  # Add F Values, p-values, and optional effect sizes to each table
  for (t_name in names(tables)) {
    tab <- tables[[t_name]]
    ms_vec <- tab[["Mean Square"]]
    df_vec <- tab[["DF"]]
    ss_vec <- tab[["Sum of Squares"]]

    f_vec <- ifelse(!is.na(ms_vec) & !is.na(err_ms) & err_ms > 0 & df_vec > 0,
                    ms_vec / err_ms, NA_real_)
    p_vec <- ifelse(!is.na(f_vec) & err_df > 0,
                    stats::pf(f_vec, df_vec, err_df, lower.tail = FALSE), NA_real_)

    tab$`F Value` <- f_vec
    tab$`Pr > F` <- p_vec

    if (effectsize) {
      # Partial and Semipartial Effect Sizes
      part_eta <- ifelse(ss_vec + sse > 0, ss_vec / (ss_vec + sse), NA_real_)
      semipart_eta <- ifelse(sst > 0, ss_vec / sst, NA_real_)

      part_omega <- ifelse(ss_vec + (n - df_vec) * mse > 0,
                           (ss_vec - df_vec * mse) / (ss_vec + (n - df_vec) * mse), NA_real_)
      part_omega <- pmax(0, part_omega)

      semipart_omega <- ifelse(sst + mse > 0,
                               (ss_vec - df_vec * mse) / (sst + mse), NA_real_)
      semipart_omega <- pmax(0, semipart_omega)

      # Noncentrality parameter estimates
      ncp_unbiased <- ifelse(!is.na(f_vec) & err_df > 2,
                             df_vec * ((err_df - 2) / err_df) * f_vec - df_vec, NA_real_)
      ncp_low_mse <- ifelse(!is.na(f_vec) & err_df > 4,
                            df_vec * ((err_df - 4) / err_df) * f_vec - df_vec * ((err_df - 4) / (err_df - 2)), NA_real_)

      # Confidence intervals for NCP and Partial Eta2
      ci_lower_ncp <- numeric(n_terms)
      ci_upper_ncp <- numeric(n_terms)
      ci_lower_eta <- numeric(n_terms)
      ci_upper_eta <- numeric(n_terms)

      conf_level <- 1 - alpha
      for (i in seq_len(n_terms)) {
        f_i <- f_vec[i]
        df_i <- df_vec[i]
        if (!is.na(f_i) && f_i > 0 && df_i > 0 && err_df > 0) {
          # Lower NCP limit
          p_zero <- stats::pf(f_i, df_i, err_df, ncp = 0)
          if (p_zero <= 1 - alpha / 2) {
            ncp_l <- 0
          } else {
            f_root_l <- function(lam) stats::pf(f_i, df_i, err_df, ncp = lam) - (1 - alpha / 2)
            up_l <- max(100, f_i * df_i * 2)
            while (f_root_l(up_l) > 0 && up_l < 1e5) up_l <- up_l * 2
            ncp_l <- tryCatch(stats::uniroot(f_root_l, c(0, up_l))$root, error = function(e) 0)
          }

          # Upper NCP limit
          f_root_u <- function(lam) stats::pf(f_i, df_i, err_df, ncp = lam) - (alpha / 2)
          up_u <- max(100, f_i * df_i * 4)
          while (f_root_u(up_u) > 0 && up_u < 1e6) up_u <- up_u * 2
          ncp_u <- tryCatch(stats::uniroot(f_root_u, c(0, up_u))$root, error = function(e) NA_real_)

          ci_lower_ncp[i] <- ncp_l
          ci_upper_ncp[i] <- ncp_u
          ci_lower_eta[i] <- ncp_l / (ncp_l + n)
          ci_upper_eta[i] <- if (!is.na(ncp_u)) ncp_u / (ncp_u + n) else NA_real_
        } else {
          ci_lower_ncp[i] <- NA_real_
          ci_upper_ncp[i] <- NA_real_
          ci_lower_eta[i] <- NA_real_
          ci_upper_eta[i] <- NA_real_
        }
      }

      tab$`Partial Eta2` <- part_eta
      tab$`Partial Omega2` <- part_omega
      tab$`Semipartial Eta2` <- semipart_eta
      tab$`Semipartial Omega2` <- semipart_omega
      tab$`NCP Est` <- ncp_unbiased
      tab[[sprintf("NCP Lower %g%%", conf_level * 100)]] <- ci_lower_ncp
      tab[[sprintf("NCP Upper %g%%", conf_level * 100)]] <- ci_upper_ncp
      tab[[sprintf("Eta2 Lower %g%%", conf_level * 100)]] <- ci_lower_eta
      tab[[sprintf("Eta2 Upper %g%%", conf_level * 100)]] <- ci_upper_eta
    }

    tables[[t_name]] <- tab
  }

  res <- list(
    tables = tables,
    error_term = err_name,
    error_ms = err_ms,
    error_df = err_df,
    custom_error = custom_error,
    dependent_variable = d$response_names[1],
    effectsize = effectsize,
    alpha = alpha,
    call = match.call()
  )
  class(res) <- "anova_stat_glm"
  res
}

#' @export
print.anova_stat_glm <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nANOVA Tables for Dependent Variable:", x$dependent_variable, "\n")
  if (x$custom_error) {
    cat(sprintf("Tests evaluated with denominator: %s (DF = %g, MS = %g)\n",
                x$error_term, x$error_df, x$error_ms))
  }

  for (t_name in names(x$tables)) {
    cat("\n", t_name, " Sums of Squares:\n", sep = "")
    tab <- x$tables[[t_name]]
    stats::printCoefmat(as.matrix(tab), digits = digits, cs.ind = 1L, tst.ind = 4L,
                        zap.ps = TRUE, P.values = TRUE, has.Pvalue = TRUE, na.print = "")
  }
  cat("\n")
  invisible(x)
}
