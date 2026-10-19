#' Design Matrix and Terms Specification for General Linear Models
#'
#' @description
#' Constructs full indicator (over-parameterized) design matrices,
#' tracks column-to-term associations, identifies factor levels,
#' and evaluates effect containment relationships for hierarchical ANOVA testing.
#'
#' @name design_matrix
NULL

#' Check Containment Relationship Between Model Effects
#'
#' @description
#' Evaluates whether effect \code{term_a} is contained within effect \code{term_b}.
#' According to the general linear model principle of marginality and containment
#' (Searle, 1971):
#' \enumerate{
#'   \item The intercept is contained in every effect.
#'   \item Effect \code{a} is contained in effect \code{b} if both involve the same
#'     continuous variables and all classification variables in \code{a} are
#'     also present in \code{b}, with \code{b} having strictly more variables.
#' }
#'
#' @param vars_a Character vector of variable names in effect \code{a}.
#' @param vars_b Character vector of variable names in effect \code{b}.
#' @param var_types Named character vector indicating variable types ("factor" or "numeric").
#'
#' @return Logical scalar indicating if effect \code{a} is contained in effect \code{b}.
#' @export
is_contained <- function(vars_a, vars_b, var_types) {
  # Intercept (empty vars_a) is contained in all non-intercept effects
  if (length(vars_a) == 0) {
    return(length(vars_b) > 0)
  }
  if (length(vars_b) <= length(vars_a)) {
    return(FALSE)
  }

  # Partition into factors and numeric
  num_a <- vars_a[var_types[vars_a] == "numeric"]
  fac_a <- vars_a[var_types[vars_a] == "factor"]

  num_b <- vars_b[var_types[vars_b] == "numeric"]
  fac_b <- vars_b[var_types[vars_b] == "factor"]

  # Continuous variables must match exactly
  if (length(num_a) != length(num_b) || !all(sort(num_a) == sort(num_b))) {
    return(FALSE)
  }

  # All classification variables in a must be present in b
  all(fac_a %in% fac_b)
}

#' Construct Design Matrices for General Linear Models
#'
#' @description
#' Builds the over-parameterized (cell indicator) design matrix and the
#' orthonormal/sum-contrast design matrix for a formula and data frame.
#' Retains complete mapping between columns, factor levels, and model terms.
#'
#' @param formula A model formula.
#' @param data A data frame containing the variables.
#' @param ref_levels Optional named list of factor levels to enforce for new data.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{X_ind}: The overparameterized indicator design matrix.
#'   \item \code{X_sum}: The full-rank sum-contrast design matrix.
#'   \item \code{y}: The response variable (vector or matrix).
#'   \item \code{response_names}: Names of the response variable(s).
#'   \item \code{terms}: The \code{terms} object.
#'   \item \code{assign_ind}: Column-to-term mapping for \code{X_ind}.
#'   \item \code{assign_sum}: Column-to-term mapping for \code{X_sum}.
#'   \item \code{term_labels}: Character vector of term names.
#'   \item \code{term_vars}: List of variable names for each term.
#'   \item \code{var_levels}: Named list of factor levels.
#'   \item \code{var_types}: Named vector of variable types.
#'   \item \code{n_obs}: Number of complete observations.
#' }
#'
#' @export
build_design_glm <- function(formula, data, ref_levels = NULL) {
  if (!inherits(formula, "formula")) {
    stop("'formula' must be a valid formula object.")
  }
  if (!is.data.frame(data)) {
    data <- as.data.frame(data)
  }

  # If reference levels are provided, apply them before creating model frame
  if (!is.null(ref_levels)) {
    for (nm in names(ref_levels)) {
      if (nm %in% names(data)) {
        data[[nm]] <- factor(data[[nm]], levels = ref_levels[[nm]])
      }
    }
  }

  # Evaluate model frame to handle transformations and missing values
  drop_unused <- is.null(ref_levels)
  mf <- stats::model.frame(formula, data = data, drop.unused.levels = drop_unused)
  mt <- stats::terms(mf)
  has_intercept <- attr(mt, "intercept") == 1
  term_labels <- attr(mt, "term.labels")
  n_terms <- length(term_labels)
  n_obs <- nrow(mf)

  # Extract response
  y <- stats::model.response(mf)
  if (is.matrix(y)) {
    response_names <- colnames(y)
    if (is.null(response_names)) {
      response_names <- paste0("Y", seq_len(ncol(y)))
      colnames(y) <- response_names
    }
  } else {
    response_names <- deparse(formula[[2]])
  }

  # Identify predictors
  pred_names <- all.vars(formula[[3]])
  var_types <- character(length(pred_names))
  names(var_types) <- pred_names
  var_levels <- list()

  for (v in pred_names) {
    val <- mf[[v]]
    if (!is.null(ref_levels) && v %in% names(ref_levels)) {
      val_factor <- factor(val, levels = ref_levels[[v]])
      mf[[v]] <- val_factor
      var_types[v] <- "factor"
      var_levels[[v]] <- ref_levels[[v]]
    } else if (is.character(val) || is.factor(val) || is.logical(val)) {
      val_factor <- as.factor(val)
      mf[[v]] <- val_factor
      var_types[v] <- "factor"
      var_levels[[v]] <- levels(val_factor)
    } else if (is.numeric(val)) {
      var_types[v] <- "numeric"
      var_levels[[v]] <- NULL
    } else {
      # Default fallback to factor
      val_factor <- as.factor(val)
      mf[[v]] <- val_factor
      var_types[v] <- "factor"
      var_levels[[v]] <- levels(val_factor)
    }
  }

  # Parse variables in each term
  term_vars <- list()
  for (i in seq_along(term_labels)) {
    lbl <- term_labels[i]
    vars_in_term <- all.vars(stats::as.formula(paste("~", lbl)))
    term_vars[[lbl]] <- vars_in_term
  }

  # 1. Build Overparameterized Indicator Design Matrix (X_ind)
  col_list <- list()
  col_names_list <- character(0)
  assign_ind <- integer(0)

  if (has_intercept) {
    col_list[[1]] <- matrix(1, n_obs, 1)
    col_names_list <- c(col_names_list, "Intercept")
    assign_ind <- c(assign_ind, 0L)
  }

  for (i in seq_along(term_labels)) {
    t_lbl <- term_labels[i]
    t_vars <- term_vars[[t_lbl]]

    # Compute columns for this term
    # Separate factors and numeric
    f_vars <- t_vars[var_types[t_vars] == "factor"]
    num_vars <- t_vars[var_types[t_vars] == "numeric"]

    # Product of numeric components
    if (length(num_vars) > 0) {
      num_mat <- mf[[num_vars[1]]]
      if (length(num_vars) > 1) {
        for (nv in num_vars[-1]) {
          num_mat <- num_mat * mf[[nv]]
        }
      }
    } else {
      num_mat <- NULL
    }

    if (length(f_vars) == 0) {
      # Pure numeric term
      col_list[[length(col_list) + 1]] <- matrix(num_mat, n_obs, 1)
      col_names_list <- c(col_names_list, t_lbl)
      assign_ind <- c(assign_ind, i)
    } else {
      # Factors present: build indicator columns for all level combinations
      grid_df <- expand.grid(var_levels[f_vars], stringsAsFactors = FALSE)
      n_combs <- nrow(grid_df)

      comb_cols <- matrix(0, n_obs, n_combs)
      comb_names <- character(n_combs)

      for (c_idx in seq_len(n_combs)) {
        match_idx <- rep(TRUE, n_obs)
        c_label_parts <- character(length(f_vars))
        for (f_i in seq_along(f_vars)) {
          fv <- f_vars[f_i]
          target_lvl <- grid_df[c_idx, f_i]
          match_idx <- match_idx & (mf[[fv]] == target_lvl)
          c_label_parts[f_i] <- paste0(fv, " ", target_lvl)
        }
        val_ind <- as.numeric(match_idx)
        if (!is.null(num_mat)) {
          val_ind <- val_ind * num_mat
          comb_names[c_idx] <- paste(c(c_label_parts, num_vars), collapse = "*")
        } else {
          comb_names[c_idx] <- paste(c_label_parts, collapse = "*")
        }
        comb_cols[, c_idx] <- val_ind
      }

      col_list[[length(col_list) + 1]] <- comb_cols
      col_names_list <- c(col_names_list, comb_names)
      assign_ind <- c(assign_ind, rep(i, n_combs))
    }
  }

  X_ind <- do.call(cbind, col_list)
  colnames(X_ind) <- col_names_list
  rownames(X_ind) <- rownames(mf)

  # 2. Build Full-Rank Sum-Contrast Matrix (X_sum) for fast Type III testing
  contrast_list <- list()
  for (v in pred_names) {
    if (var_types[v] == "factor" && nlevels(mf[[v]]) > 1) {
      contrast_list[[v]] <- stats::contr.sum(nlevels(mf[[v]]))
    }
  }

  contr_arg <- if (length(contrast_list) > 0) contrast_list else NULL
  X_sum <- stats::model.matrix(mt, data = mf, contrasts.arg = contr_arg)
  assign_sum <- attr(X_sum, "assign")

  list(
    X_ind = X_ind,
    X_sum = X_sum,
    y = y,
    response_names = response_names,
    terms = mt,
    model_frame = mf,
    assign_ind = assign_ind,
    assign_sum = assign_sum,
    term_labels = term_labels,
    term_vars = term_vars,
    var_levels = var_levels,
    var_types = var_types,
    has_intercept = has_intercept,
    n_obs = n_obs
  )
}
