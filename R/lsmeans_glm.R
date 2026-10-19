#' Least-Squares Means (Population Marginal Means) and Group Comparisons
#'
#' @description
#' Computes least-squares means (LS-means, also known as population marginal
#' means) and pairwise differences with multiple comparison adjustments
#' (Tukey-Kramer, Bonferroni, Scheffe, Sidak). Constructs estimable linear
#' combinations according to standard general linear model definitions.
#'
#' @param object An object of class \code{"stat_glm"}.
#' @param effect Character string specifying the classification effect for which
#'   LS-means are to be computed (e.g., \code{"A"} or \code{"A:B"}).
#' @param by Optional character vector of conditioning variables.
#' @param at Optional named list specifying custom values at which covariates
#'   should be evaluated. Defaults to the overall sample mean.
#' @param adjust Multiple comparison adjustment method: \code{"none"}, \code{"tukey"}
#'   (Tukey-Kramer, default), \code{"bonferroni"}, \code{"scheffe"}, or \code{"sidak"}.
#' @param alpha Significance level for confidence intervals. Defaults to \code{0.05}.
#' @param tol Numerical tolerance for estimability checks. Defaults to \code{1e-7}.
#'
#' @return An object of class \code{"lsmeans_stat_glm"} containing:
#' \itemize{
#'   \item \code{lsmeans}: Data frame of LS-means, standard errors, t-values, and confidence intervals.
#'   \item \code{diffs}: Data frame of pairwise differences, standard errors, t-values, and adjusted p-values.
#'   \item \code{effect}: The evaluated effect.
#'   \item \code{adjust}: The adjustment method used.
#'   \item \code{L}: The matrix of linear combinations defining the LS-means.
#' }
#'
#' @references
#' Goodnight, J. H., & Harvey, W. R. (1978). Least-Squares Means in the Fixed Effects
#' General Linear Models. \emph{SAS Technical Report R-103}.
#'
#' Searle, S. R., Speed, F. M., & Milliken, G. A. (1980). Population Marginal Means
#' in the Linear Model: An Alternative to Least Squares Means.
#' \emph{The American Statistician}, 34(4), 216-221.
#'
#' @examples
#' data(warpbreaks)
#' fit <- stat_glm(breaks ~ wool * tension, data = warpbreaks)
#'
#' # LS-means for main effect
#' lsmeans_glm(fit, effect = "wool")
#'
#' # LS-means for interaction with Tukey adjustment
#' lsmeans_glm(fit, effect = "wool:tension", adjust = "tukey")
#'
#' @export
lsmeans_glm <- function(object, effect, by = NULL, at = list(),
                        adjust = c("tukey", "none", "bonferroni", "scheffe", "sidak"),
                        alpha = 0.05, tol = 1e-7) {
  if (!inherits(object, "stat_glm")) {
    stop("'object' must be an object of class 'stat_glm'.")
  }
  adjust <- match.arg(adjust)

  d <- object$design
  term_labels <- d$term_labels
  var_types <- d$var_types
  var_levels <- d$var_levels
  mf <- d$model_frame
  assign_ind <- d$assign_ind
  p <- ncol(d$X_ind)

  # Validate effect: match either "A" or "A*B" or "A:B"
  clean_name <- function(x) gsub("[*:]", ":", x)
  norm_term_labels <- clean_name(term_labels)
  norm_effect <- clean_name(effect)

  effect_idx <- match(norm_effect, norm_term_labels)
  if (is.na(effect_idx)) {
    # Check if effect is a single variable or interaction of variables
    eff_vars <- unlist(strsplit(norm_effect, ":"))
    if (!all(eff_vars %in% names(var_types))) {
      stop(sprintf("Effect '%s' not recognized in model terms: %s",
                   effect, paste(term_labels, collapse = ", ")))
    }
  } else {
    eff_vars <- d$term_vars[[term_labels[effect_idx]]]
  }

  # Ensure all variables in effect are classification factors
  if (any(var_types[eff_vars] != "factor")) {
    stop("LS-means can only be computed for classification effects (factors).")
  }

  # Compute covariate reference values (at user values or sample means)
  num_vars <- names(var_types)[var_types == "numeric"]
  cov_means <- list()
  for (nv in num_vars) {
    if (nv %in% names(at)) {
      cov_means[[nv]] <- as.numeric(at[[nv]])
    } else {
      cov_means[[nv]] <- mean(mf[[nv]], na.rm = TRUE)
    }
  }

  # Expand grid of levels for the evaluated effect
  effect_grid <- expand.grid(var_levels[eff_vars], stringsAsFactors = FALSE)
  n_levels <- nrow(effect_grid)

  # Labels for each level of the effect
  level_names <- character(n_levels)
  for (i in seq_len(n_levels)) {
    parts <- character(length(eff_vars))
    for (j in seq_along(eff_vars)) {
      parts[j] <- paste(eff_vars[j], effect_grid[i, j])
    }
    level_names[i] <- paste(parts, collapse = " ")
  }

  # Construct L matrix for LS-means: dimension n_levels x p
  L_mat <- matrix(0, nrow = n_levels, ncol = p)
  colnames(L_mat) <- colnames(d$X_ind)
  rownames(L_mat) <- level_names

  # Intercept
  if (d$has_intercept) {
    int_col <- which(assign_ind == 0L)
    L_mat[, int_col] <- 1
  }

  # For each term in model, assign weights according to LSMEANS rules
  for (k in seq_along(term_labels)) {
    t_lbl <- term_labels[k]
    t_vars <- d$term_vars[[t_lbl]]
    cols_k <- which(assign_ind == k)
    cnames_k <- colnames(d$X_ind)[cols_k]

    f_in_term <- t_vars[var_types[t_vars] == "factor"]
    num_in_term <- t_vars[var_types[t_vars] == "numeric"]

    # Continuous scaling factor
    cov_scale <- 1
    if (length(num_in_term) > 0) {
      for (nv in num_in_term) {
        cov_scale <- cov_scale * cov_means[[nv]]
      }
    }

    if (length(f_in_term) == 0) {
      # Pure continuous covariate: set to cov_scale for all rows
      L_mat[, cols_k] <- cov_scale
    } else {
      # Term contains factors
      grid_term <- expand.grid(var_levels[f_in_term], stringsAsFactors = FALSE)
      n_term_combs <- nrow(grid_term)

      # Determine relationship with evaluated effect:
      # Which factors are shared?
      shared_f <- intersect(f_in_term, eff_vars)
      other_f <- setdiff(f_in_term, eff_vars)

      for (row_i in seq_len(n_levels)) {
        target_eff_vals <- effect_grid[row_i, , drop = FALSE]

        if (length(shared_f) == 0) {
          # Unrelated factor: equal weights 1 / n_levels_term
          weight <- cov_scale / n_term_combs
          L_mat[row_i, cols_k] <- weight
        } else {
          # Shares factors with effect
          # Check which columns in term match target_eff_vals on the shared factors
          match_term <- rep(TRUE, n_term_combs)
          for (sf in shared_f) {
            match_term <- match_term & (grid_term[[sf]] == target_eff_vals[[sf]])
          }

          n_matches <- sum(match_term)
          if (n_matches > 0) {
            weight <- cov_scale / n_matches
            matched_cols <- cols_k[match_term]
            L_mat[row_i, matched_cols] <- weight
          }
        }
      }
    }
  }

  # Verify Estimability
  ginv_fn <- if (object$ginv_method == "g2") ginv_g2 else ginv_svd
  est_flags <- is_estimable(L_mat, object$XtX, tol = tol, ginv_fn = ginv_fn)

  b <- object$coefficients
  G <- object$G
  mse <- object$mse
  dfe <- object$df.residual

  est_lsm <- as.vector(L_mat %*% b)
  cov_lsm <- as.matrix((L_mat %*% G %*% t(L_mat)) * mse)
  var_lsm <- diag(cov_lsm)
  se_lsm <- sqrt(pmax(0, var_lsm))

  est_lsm[!est_flags] <- NA_real_
  se_lsm[!est_flags] <- NA_real_

  t_lsm <- est_lsm / se_lsm
  p_lsm <- 2 * stats::pt(abs(t_lsm), df = dfe, lower.tail = FALSE)

  t_crit <- stats::qt(1 - alpha / 2, df = dfe)
  ci_lower <- est_lsm - t_crit * se_lsm
  ci_upper <- est_lsm + t_crit * se_lsm

  conf_pct <- sprintf("%g%%", (1 - alpha) * 100)
  lsm_df <- data.frame(
    Effect = level_names,
    Estimate = est_lsm,
    `Std Error` = se_lsm,
    DF = rep(dfe, n_levels),
    `t Value` = t_lsm,
    `Pr > |t|` = p_lsm,
    `Lower CL` = ci_lower,
    `Upper CL` = ci_upper,
    check.names = FALSE,
    row.names = seq_len(n_levels)
  )
  names(lsm_df)[7:8] <- c(paste("Lower", conf_pct), paste("Upper", conf_pct))

  # ---------------------------------------------------------
  # Pairwise Differences (Diffs)
  # ---------------------------------------------------------
  diff_rows <- list()
  diff_names <- character(0)
  pair_indices <- list()

  if (n_levels > 1) {
    for (i in 1:(n_levels - 1)) {
      for (j in (i + 1):n_levels) {
        d_vec <- L_mat[i, ] - L_mat[j, ]
        diff_rows[[length(diff_rows) + 1]] <- d_vec
        diff_names <- c(diff_names, paste(level_names[i], "-", level_names[j]))
        pair_indices[[length(pair_indices) + 1]] <- c(i, j)
      }
    }
  }

  if (length(diff_rows) > 0) {
    L_diff <- do.call(rbind, diff_rows)
    est_diff_flags <- is_estimable(L_diff, object$XtX, tol = tol, ginv_fn = ginv_fn)

    diff_est <- as.vector(L_diff %*% b)
    diff_var <- diag(L_diff %*% G %*% t(L_diff)) * mse
    diff_se <- sqrt(pmax(0, diff_var))

    diff_est[!est_diff_flags] <- NA_real_
    diff_se[!est_diff_flags] <- NA_real_

    diff_t <- diff_est / diff_se
    raw_p <- 2 * stats::pt(abs(diff_t), df = dfe, lower.tail = FALSE)

    # Multiplicity adjustments
    m_diffs <- length(diff_rows)
    adj_p <- raw_p

    if (adjust == "tukey") {
      # Tukey-Kramer adjustment via Studentized range distribution
      q_val <- abs(diff_t) * sqrt(2)
      adj_p <- stats::ptukey(q_val, nmeans = n_levels, df = dfe, lower.tail = FALSE)
    } else if (adjust == "bonferroni") {
      adj_p <- pmin(1, raw_p * m_diffs)
    } else if (adjust == "scheffe") {
      f_scheffe <- (diff_t^2) / (n_levels - 1)
      adj_p <- stats::pf(f_scheffe, df1 = n_levels - 1, df2 = dfe, lower.tail = FALSE)
    } else if (adjust == "sidak") {
      adj_p <- 1 - (1 - raw_p)^m_diffs
    }

    # Confidence limits for differences
    if (adjust == "tukey") {
      q_crit <- stats::qtukey(1 - alpha, nmeans = n_levels, df = dfe) / sqrt(2)
      diff_lower <- diff_est - q_crit * diff_se
      diff_upper <- diff_est + q_crit * diff_se
    } else if (adjust == "bonferroni") {
      t_bonf <- stats::qt(1 - (alpha / (2 * m_diffs)), df = dfe)
      diff_lower <- diff_est - t_bonf * diff_se
      diff_upper <- diff_est + t_bonf * diff_se
    } else {
      diff_lower <- diff_est - t_crit * diff_se
      diff_upper <- diff_est + t_crit * diff_se
    }

    diff_df <- data.frame(
      Comparison = diff_names,
      Estimate = diff_est,
      `Std Error` = diff_se,
      DF = rep(dfe, m_diffs),
      `t Value` = diff_t,
      `Pr > |t|` = raw_p,
      `Adj P` = adj_p,
      `Lower CL` = diff_lower,
      `Upper CL` = diff_upper,
      check.names = FALSE,
      row.names = seq_len(m_diffs)
    )
    names(diff_df)[8:9] <- c(paste("Lower", conf_pct), paste("Upper", conf_pct))
  } else {
    diff_df <- data.frame()
  }

  res <- list(
    lsmeans = lsm_df,
    diffs = diff_df,
    effect = effect,
    adjust = adjust,
    alpha = alpha,
    L = L_mat,
    vcov = cov_lsm,
    grid = effect_grid,
    df.residual = dfe
  )
  class(res) <- "lsmeans_stat_glm"
  res
}

#' @export
print.lsmeans_stat_glm <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nLeast-Squares Means for Effect:", x$effect, "\n\n")
  tab_lsm <- x$lsmeans[, -1, drop = FALSE]
  rownames(tab_lsm) <- x$lsmeans$Effect
  stats::printCoefmat(as.matrix(tab_lsm), digits = digits, cs.ind = 1:2, tst.ind = 4L,
                      zap.ps = TRUE, P.values = TRUE, has.Pvalue = TRUE, na.print = "Non-est")

  if (nrow(x$diffs) > 0) {
    cat(sprintf("\nDifferences of Least-Squares Means (Adjustment: %s):\n\n", toupper(x$adjust)))
    tab_diff <- x$diffs[, -1, drop = FALSE]
    rownames(tab_diff) <- x$diffs$Comparison
    stats::printCoefmat(as.matrix(tab_diff), digits = digits, cs.ind = 1:2, tst.ind = 4L,
                        zap.ps = TRUE, P.values = TRUE, has.Pvalue = TRUE, na.print = "Non-est")
  }
  cat("\n")
  invisible(x)
}
