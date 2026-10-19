#' Mixed Model Analysis of Variance and Expected Mean Squares
#'
#' @description
#' Implements the SAS \code{PROC GLM} \code{RANDOM} statement. Computes the
#' Expected Mean Squares (EMS) coefficients for model terms using the foundational
#' Goodnight & Speed (1978) formulation, and under \code{test = TRUE}, performs
#' mixed model hypothesis tests using appropriate error terms (exact single-term
#' denominators or Satterthwaite pseudo-F linear combinations of mean squares).
#'
#' @param object An object of class \code{"stat_glm"}.
#' @param effects An optional character vector, one-sided formula (e.g.
#'   \code{~ person + machine:person}), or comma-separated string declaring which
#'   model effects are random. If omitted, uses \code{object$random_effects}.
#' @param test Logical; if \code{TRUE} (default, matching the SAS GLM \code{TEST}
#'   option), performs hypothesis tests for each model effect using appropriate
#'   error terms determined from expected mean squares.
#' @param q Logical; if \code{TRUE} (matching the SAS GLM \code{Q} option),
#'   computes and includes the quadratic forms in fixed effects. Defaults to \code{FALSE}.
#' @param type Integer specifying the sum of squares type to use for expected
#'   mean squares and tests (1 or 3). Defaults to \code{3} (Type III), matching SAS PROC GLM.
#' @param tol Tolerance for rank, zero-coefficient detection, and linear system solves.
#'   Defaults to \code{1e-7}.
#' @param ... Additional arguments (currently unused).
#'
#' @return An object of class \code{"random_stat_glm"} containing:
#'   \item{ems_matrix}{Matrix of Expected Mean Squares coefficients for Error and random effects.}
#'   \item{ems_table}{Data frame summarizing the symbolic Expected Mean Squares for each effect.}
#'   \item{tests}{Data frame of mixed model ANOVA hypothesis tests with numerator and denominator
#'     degrees of freedom, sums of squares, mean squares, F-values, p-values, and error term definitions (if \code{test = TRUE}).}
#'   \item{q_matrices}{List of numeric quadratic form matrices for fixed effects (if \code{q = TRUE}).}
#'   \item{random_effects}{Character vector of effects declared random.}
#'   \item{fixed_effects}{Character vector of effects remaining fixed.}
#'   \item{type}{Sum of squares type used (1 or 3).}
#'   \item{notes}{Character vector of diagnostic notes or model caveats.}
#'
#' @references
#' Goodnight, J. H., and Speed, F. M. (1978). Computing Expected Mean Squares.
#' \emph{Communications in Statistics - Theory and Methods}, 7(13), 1231-1245.
#'
#' Milliken, G. A., and Johnson, D. E. (1984). \emph{Analysis of Messy Data, Volume 1: Designed Experiments}.
#' Van Nostrand Reinhold, New York.
#'
#' Satterthwaite, F. E. (1946). An Approximate Distribution of Estimates of Variance Components.
#' \emph{Biometrics Bulletin}, 2(6), 110-114.
#'
#' @examples
#' # Milliken & Johnson (1984) Unbalanced Mixed Model (Machine-Person Experiment)
#' m <- c(
#'   1, 1, 52.0, 1, 2, 51.8, 1, 2, 52.8, 1, 3, 60.0, 1, 4, 51.1, 1, 4, 52.3, 1, 5, 50.9,
#'   1, 5, 51.8, 1, 5, 51.4, 1, 6, 46.4, 1, 6, 44.8, 1, 6, 49.2, 2, 1, 64.0, 2, 2, 59.7,
#'   2, 2, 60.0, 2, 2, 59.0, 2, 3, 68.6, 2, 3, 65.8, 2, 4, 63.2, 2, 4, 62.8, 2, 4, 62.2,
#'   2, 5, 64.8, 2, 5, 65.0, 2, 6, 43.7, 2, 6, 44.2, 2, 6, 43.0, 3, 1, 67.5, 3, 1, 67.2,
#'   3, 1, 66.9, 3, 2, 61.5, 3, 2, 61.7, 3, 2, 62.3, 3, 3, 70.8, 3, 3, 70.6, 3, 3, 71.0,
#'   3, 4, 64.1, 3, 4, 66.2, 3, 4, 64.0, 3, 5, 72.1, 3, 5, 72.0, 3, 5, 71.1, 3, 6, 62.0,
#'   3, 6, 61.4, 3, 6, 60.5
#' )
#' mat <- matrix(m, ncol = 3, byrow = TRUE)
#' dat <- data.frame(
#'   machine = factor(mat[, 1]),
#'   person = factor(mat[, 2]),
#'   rating = mat[, 3]
#' )
#'
#' # Fit model and declare random effects
#' fit <- stat_glm(rating ~ machine * person, data = dat)
#' rnd <- random_glm(fit, effects = ~ person + machine:person, test = TRUE)
#' print(rnd)
#'
#' @export
random_glm <- function(object, effects = NULL, test = TRUE, q = FALSE,
                       type = 3, tol = 1e-7, ...) {
  if (!inherits(object, "stat_glm")) {
    stop("'object' must be an object of class 'stat_glm'.")
  }
  if (inherits(object, "stat_mlm")) {
    stop("Random-effects analysis is currently supported only for univariate models.")
  }

  d <- object$design
  term_labels <- d$term_labels
  n_terms <- length(term_labels)
  if (n_terms == 0) {
    stop("The fitted model does not contain any independent predictor terms.")
  }

  # Resolve effects argument
  if (is.null(effects)) {
    effects <- object$random_effects
  }
  if (is.null(effects)) {
    stop("No random effects specified. Please provide 'effects' (e.g. effects = ~ person + machine:person).")
  }

  random_effects <- parse_random_effects(effects, term_labels)
  fixed_effects <- setdiff(term_labels, random_effects)

  if (!type %in% c(1, 3)) {
    stop("'type' must be either 1 or 3.")
  }

  # Extract design matrices and weights
  X_ind <- d$X_ind
  X_sum <- d$X_sum
  weights <- object$weights
  w_sqrt <- sqrt(weights)
  assign_ind <- d$assign_ind
  assign_sum <- d$assign_sum

  # Compute standard ANOVA table for requested type (pass random = FALSE to prevent recursion)
  an_res <- anova_glm(object, type = type, random = FALSE)
  type_label <- if (type == 1) "Type I" else "Type III"
  tab_ss <- an_res$tables[[type_label]]

  # Goodnight & Speed (1978) formulation
  # 1. Full-rank sum-contrast quadratic form components
  X_s_w <- X_sum * w_sqrt
  XtX_sum <- crossprod(X_s_w)
  XtX_inv <- tryCatch(solve(XtX_sum), error = function(e) ginv_svd(XtX_sum, tol = tol))

  X_i_w <- X_ind * w_sqrt
  XtX_ind_cross <- crossprod(X_s_w, X_i_w)

  # EMS coefficient matrix: Rows = model terms; Cols = Error, random_effects
  ems_cols <- c("Error", random_effects)
  ems_matrix <- matrix(0, nrow = n_terms, ncol = length(ems_cols),
                      dimnames = list(term_labels, ems_cols))
  ems_matrix[, "Error"] <- 1.0 # Every mean square expectation includes 1.0 * Var(Error)

  # Full coefficients across all terms (including fixed) for Q detection
  all_coeffs <- matrix(0, nrow = n_terms, ncol = n_terms,
                       dimnames = list(term_labels, term_labels))

  q_matrices <- list()

  for (k in seq_len(n_terms)) {
    cols_k <- which(assign_sum == k)
    rk_k <- length(cols_k)

    if (rk_k > 0) {
      L_k <- XtX_inv[cols_k, , drop = FALSE] %*% XtX_ind_cross
      V_k <- XtX_inv[cols_k, cols_k, drop = FALSE]
      V_k_inv <- quad_form_inv(V_k, tol = tol)

      # For each term j in the model:
      for (j in seq_len(n_terms)) {
        cols_j <- which(assign_ind == j)
        if (length(cols_j) > 0) {
          L_kj <- L_k[, cols_j, drop = FALSE]
          # trace(V_k_inv %*% L_kj %*% t(L_kj)) = sum((t(L_kj) %*% V_k_inv) * t(L_kj))
          val <- sum((t(L_kj) %*% V_k_inv) * t(L_kj)) / rk_k
          all_coeffs[k, j] <- if (abs(val) > tol) val else 0
          if (term_labels[j] %in% random_effects) {
            ems_matrix[k, term_labels[j]] <- if (abs(val) > tol) val else 0
          }
        }
      }

      if (q) {
        # Quadratic form matrix on beta_ind: (1 / rk_k) * L_k' V_k_inv L_k
        Q_k <- (crossprod(L_k, V_k_inv %*% L_k)) / rk_k
        colnames(Q_k) <- colnames(X_ind)
        rownames(Q_k) <- colnames(X_ind)
        q_matrices[[term_labels[k]]] <- Q_k
      }
    }
  }

  # Format symbolic Expected Mean Squares
  ems_expressions <- character(n_terms)
  names(ems_expressions) <- term_labels

  for (k in seq_len(n_terms)) {
    t_name <- term_labels[k]
    expr_parts <- "Var(Error)"

    # Add random effects
    for (r_name in random_effects) {
      cf <- ems_matrix[t_name, r_name]
      if (cf > tol) {
        if (abs(cf - 1.0) < tol) {
          expr_parts <- paste0(expr_parts, sprintf(" + Var(%s)", r_name))
        } else {
          expr_parts <- paste0(expr_parts, sprintf(" + %0.4g Var(%s)", cf, r_name))
        }
      }
    }

    # Add fixed effects appearing in Q
    fixed_in_q <- character(0)
    for (f_name in fixed_effects) {
      if (all_coeffs[t_name, f_name] > tol) {
        fixed_in_q <- c(fixed_in_q, f_name)
      }
    }
    if (length(fixed_in_q) > 0) {
      expr_parts <- paste0(expr_parts, sprintf(" + Q(%s)", paste(fixed_in_q, collapse = ", ")))
    }

    ems_expressions[k] <- expr_parts
  }

  ems_table <- data.frame(
    Source = term_labels,
    `Expected Mean Square` = ems_expressions,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  # ---------------------------------------------------------
  # Mixed Model ANOVA Tests (TEST Option)
  # ---------------------------------------------------------
  tests_df <- NULL
  notes <- character(0)

  if (test) {
    # Full EMS rows including Error
    ems_all <- rbind(
      ems_matrix,
      Error = c(1.0, rep(0, length(random_effects)))
    )

    ms_vec <- c(tab_ss[["Mean Square"]], object$mse)
    df_vec <- c(tab_ss[["DF"]], object$df.residual)
    names(ms_vec) <- c(term_labels, "Error")
    names(df_vec) <- c(term_labels, "Error")

    num_df_vals <- integer(n_terms)
    num_ss_vals <- numeric(n_terms)
    num_ms_vals <- numeric(n_terms)
    f_vals <- numeric(n_terms)
    p_vals <- numeric(n_terms)
    denom_ms_vals <- numeric(n_terms)
    denom_df_vals <- numeric(n_terms)
    denom_ss_vals <- numeric(n_terms)
    denom_labels <- character(n_terms)
    has_satterthwaite <- FALSE
    has_fixed_q_assumption <- FALSE

    for (k in seq_len(n_terms)) {
      term_k <- term_labels[k]
      num_df <- df_vec[term_k]
      num_ss <- tab_ss[term_k, "Sum of Squares"]
      num_ms <- ms_vec[term_k]

      num_df_vals[k] <- num_df
      num_ss_vals[k] <- num_ss
      num_ms_vals[k] <- num_ms

      # Check if other fixed effects are in Q
      other_fixed_in_q <- setdiff(fixed_effects[all_coeffs[term_k, fixed_effects] > tol], term_k)
      if (length(other_fixed_in_q) > 0) {
        has_fixed_q_assumption <- TRUE
      }

      # Target expectation under H0
      target <- ems_all[term_k, ]
      if (term_k %in% random_effects) {
        target[term_k] <- 0
      }

      # Candidate terms to form denominator: random effects (excluding term_k) and Error
      cand_terms <- c(random_effects, "Error")
      cand_terms <- setdiff(cand_terms, term_k)

      # Candidate EMS matrix (rows: cand_terms, cols: Error, random_effects)
      E_cand <- ems_all[cand_terms, , drop = FALSE]

      # Solve linear system: t(E_cand) %*% a = target
      a <- tryCatch({
        as.vector(solve(t(E_cand), target))
      }, error = function(e) {
        as.vector(ginv_svd(t(E_cand), tol = tol) %*% target)
      })
      names(a) <- cand_terms

      # Verify solution quality
      residual_norm <- sum(abs(as.vector(t(E_cand) %*% a) - target))
      if (residual_norm > 1e-4) {
        # No exact error term found
        denom_labels[k] <- "No suitable error term"
        denom_df_vals[k] <- NA_real_
        denom_ms_vals[k] <- NA_real_
        denom_ss_vals[k] <- NA_real_
        f_vals[k] <- NA_real_
        p_vals[k] <- NA_real_
      } else {
        # Extract non-zero weights
        nz <- abs(a) > tol
        a_nz <- a[nz]
        cand_nz <- cand_terms[nz]

        if (length(a_nz) == 0) {
          denom_labels[k] <- "None"
          denom_df_vals[k] <- NA_real_
          denom_ms_vals[k] <- NA_real_
          denom_ss_vals[k] <- NA_real_
          f_vals[k] <- NA_real_
          p_vals[k] <- NA_real_
        } else if (length(a_nz) == 1 && abs(a_nz[1] - 1.0) < tol) {
          # Exact single effect denominator
          s_name <- cand_nz[1]
          d_df <- df_vec[s_name]
          d_ms <- ms_vec[s_name]
          d_ss <- d_df * d_ms

          denom_labels[k] <- paste0("MS(", s_name, ")")
          denom_df_vals[k] <- d_df
          denom_ms_vals[k] <- d_ms
          denom_ss_vals[k] <- d_ss

          f_val <- if (!is.na(d_ms) && d_ms > 0) num_ms / d_ms else NA_real_
          p_val <- if (!is.na(f_val) && d_df > 0) stats::pf(f_val, num_df, d_df, lower.tail = FALSE) else NA_real_

          f_vals[k] <- f_val
          p_vals[k] <- p_val
        } else {
          # Satterthwaite linear combination of mean squares
          has_satterthwaite <- TRUE
          d_ms <- sum(a_nz * ms_vec[cand_nz])
          v_parts <- (a_nz * ms_vec[cand_nz])^2 / df_vec[cand_nz]
          sum_v <- sum(v_parts)
          d_df <- if (sum_v > 0) (d_ms^2) / sum_v else NA_real_
          d_ss <- d_ms * d_df

          # Build formatted label
          terms_lbl <- character(length(a_nz))
          for (i in seq_along(a_nz)) {
            val_a <- a_nz[i]
            val_name <- cand_nz[i]
            if (abs(val_a - 1.0) < tol) {
              terms_lbl[i] <- paste0("MS(", val_name, ")")
            } else if (abs(val_a - (-1.0)) < tol) {
              terms_lbl[i] <- paste0("- MS(", val_name, ")")
            } else {
              terms_lbl[i] <- sprintf("%0.4g*MS(%s)", val_a, val_name)
            }
          }
          denom_labels[k] <- paste(terms_lbl, collapse = " + ")
          denom_labels[k] <- gsub("\\+ -", "- ", denom_labels[k])

          denom_df_vals[k] <- d_df
          denom_ms_vals[k] <- d_ms
          denom_ss_vals[k] <- d_ss

          f_val <- if (!is.na(d_ms) && d_ms > 0) num_ms / d_ms else NA_real_
          p_val <- if (!is.na(f_val) && !is.na(d_df) && d_df > 0) stats::pf(f_val, num_df, d_df, lower.tail = FALSE) else NA_real_

          f_vals[k] <- f_val
          p_vals[k] <- p_val
        }
      }
    }

    ss_col_name <- if (type == 1) "Type I SS" else "Type III SS"
    tests_df <- data.frame(
      Source = term_labels,
      DF = num_df_vals,
      `Sum of Squares` = num_ss_vals,
      `Mean Square` = num_ms_vals,
      `F Value` = f_vals,
      `Pr > F` = p_vals,
      `Error Term` = denom_labels,
      `Denom DF` = denom_df_vals,
      `Denom MS` = denom_ms_vals,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    names(tests_df)[3] <- ss_col_name

    if (has_satterthwaite) {
      notes <- c(notes, "Denominator degrees of freedom computed via Satterthwaite approximation.")
    }
    if (has_fixed_q_assumption) {
      notes <- c(notes, "Tests for some effects assume that other fixed effects involved in Q are zero.")
    }
  }

  res <- list(
    ems_matrix = ems_matrix,
    ems_table = ems_table,
    tests = tests_df,
    q_matrices = if (q) q_matrices else NULL,
    random_effects = random_effects,
    fixed_effects = fixed_effects,
    type = type,
    dependent_variable = d$response_names[1],
    notes = notes,
    call = match.call()
  )
  class(res) <- "random_stat_glm"
  res
}

#' Helper to parse random effects specifications
#' @noRd
parse_random_effects <- function(effects, term_labels) {
  if (inherits(effects, "formula")) {
    ef_vars <- attr(stats::terms(effects), "term.labels")
  } else if (is.character(effects)) {
    if (length(effects) == 1 && grepl(",", effects)) {
      ef_vars <- trimws(strsplit(effects, ",")[[1]])
    } else {
      ef_vars <- trimws(effects)
    }
  } else {
    stop("'effects' must be a character vector or a one-sided formula.")
  }

  # Normalize '*' to ':' for interaction terms
  ef_vars <- gsub("\\*", ":", ef_vars)
  norm_terms <- gsub("\\*", ":", term_labels)

  matched <- character(0)
  for (ev in ef_vars) {
    idx <- which(norm_terms == ev)
    if (length(idx) == 0) {
      # Try matching interaction parts regardless of order
      parts <- strsplit(ev, ":")[[1]]
      for (t_i in seq_along(norm_terms)) {
        t_parts <- strsplit(norm_terms[t_i], ":")[[1]]
        if (setequal(parts, t_parts)) {
          idx <- t_i
          break
        }
      }
    }
    if (length(idx) == 0) {
      stop(sprintf("Random effect '%s' not found among model terms: %s",
                   ev, paste(term_labels, collapse = ", ")))
    }
    matched <- c(matched, term_labels[idx])
  }
  unique(matched)
}

#' @export
print.random_stat_glm <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  type_str <- if (x$type == 1) "Type I" else "Type III"

  cat("\nExpected Values of", type_str, "Mean Squares\n")
  cat("Dependent Variable:", x$dependent_variable, "\n\n")

  # Print EMS table
  print(x$ems_table, row.names = FALSE)

  if (!is.null(x$tests)) {
    cat("\nTests of Hypotheses for Mixed Model Analysis of Variance (", type_str, " SS)\n\n", sep = "")

    # Format each source and its denominator nicely
    for (i in seq_len(nrow(x$tests))) {
      row <- x$tests[i, ]
      cat(sprintf("Source: %s\n", row$Source))

      # Test summary line
      ss_col <- if (x$type == 1) "Type I SS" else "Type III SS"
      f_str <- if (is.na(row$`F Value`)) "NA" else sprintf("%0.2f", row$`F Value`)
      p_str <- if (is.na(row$`Pr > F`)) "NA" else if (row$`Pr > F` < 0.0001) "<.0001" else sprintf("%0.4f", row$`Pr > F`)

      cat(sprintf("  Numerator:   DF = %d, %s = %0.4f, Mean Square = %0.4f\n",
                  row$DF, ss_col, row[[ss_col]], row$`Mean Square`))
      if (!is.na(row$`Denom DF`)) {
        cat(sprintf("  Denominator: DF = %0.3f, Mean Square = %0.4f\n",
                    row$`Denom DF`, row$`Denom MS`))
        cat(sprintf("  Error Term:  %s\n", row$`Error Term`))
      } else {
        cat(sprintf("  Error Term:  %s\n", row$`Error Term`))
      }
      cat(sprintf("  F Value = %s, Pr > F = %s\n\n", f_str, p_str))
    }

    if (length(x$notes) > 0) {
      for (nt in x$notes) {
        cat("Note:", nt, "\n")
      }
      cat("\n")
    }
  }

  if (!is.null(x$q_matrices)) {
    cat("\nQuadratic Forms in Fixed Effects (Q Matrices):\n")
    for (nm in names(x$q_matrices)) {
      cat("\n--- Effect:", nm, "---\n")
      print(round(x$q_matrices[[nm]], digits = digits))
    }
    cat("\n")
  }

  invisible(x)
}

#' @export
as.data.frame.random_stat_glm <- function(x, row.names = NULL, optional = FALSE, ...) {
  if (!is.null(x$tests)) {
    return(as.data.frame(x$tests, row.names = row.names, optional = optional, ...))
  }
  as.data.frame(x$ems_table, row.names = row.names, optional = optional, ...)
}
