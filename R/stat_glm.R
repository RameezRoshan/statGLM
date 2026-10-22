#' Fit Statistical General Linear Models
#'
#' @description
#' Fits general linear models using over-parameterized (cell indicator) design
#' matrices, solving the normal equations via generalized inverses (reflexive g2
#' or Moore-Penrose pseudo-inverse). Handles univariate responses, multivariate
#' responses (MANOVA), weighted least squares, and over-parameterized factor models.
#'
#' @param formula An object of class \code{\link[stats]{formula}} specifying the model.
#' @param data A data frame containing the variables in the model.
#' @param weights An optional numeric vector of positive observation weights.
#' @param random An optional character vector, one-sided formula (e.g. \code{~ person + machine:person}),
#'   or comma-separated string declaring which model effects are random (equivalent to
#'   the SAS GLM \code{RANDOM} statement).
#' @param ginv_method Character string specifying the generalized inverse method:
#'   \code{"g2"} (default, reflexive generalized inverse via Goodnight SWEEP with
#'   redundant parameters set to zero) or \code{"svd"} (Moore-Penrose pseudo-inverse).
#' @param tol Numerical tolerance for pivot and singularity checks. Defaults to \code{1e-7}.
#'
#' @return An object of class \code{"stat_glm"} (or \code{"stat_mlm"} for multivariate
#'   responses) containing:
#' \itemize{
#'   \item \code{coefficients}: Parameter estimates \eqn{b = (X'X)^- X'y}.
#'   \item \code{residuals}: Model residuals \eqn{y - \hat{y}}.
#'   \item \code{fitted.values}: Predicted values \eqn{\hat{y} = X b}.
#'   \item \code{rank}: Rank of the design matrix \eqn{X}.
#'   \item \code{df.residual}: Degrees of freedom for error (\eqn{n - \mathrm{rank}(X)}).
#'   \item \code{df.model}: Degrees of freedom for model (\eqn{\mathrm{rank}(X) - 1}).
#'   \item \code{df.total}: Total degrees of freedom.
#'   \item \code{sse}: Error sum of squares.
#'   \item \code{ssm}: Model sum of squares.
#'   \item \code{sst}: Corrected total sum of squares.
#'   \item \code{mse}: Mean square error.
#'   \item \code{msm}: Mean square model.
#'   \item \code{f_value}: Overall model F statistic.
#'   \item \code{p_value}: Overall model p-value.
#'   \item \code{r_square}: Coefficient of determination \eqn{R^2}.
#'   \item \code{adj_r_square}: Adjusted \eqn{R^2}.
#'   \item \code{root_mse}: Root mean square error \eqn{\sqrt{\mathrm{MSE}}}.
#'   \item \code{cv}: Coefficient of variation \eqn{100 \times \mathrm{Root MSE} / \bar{y}}.
#'   \item \code{y_mean}: Mean of the dependent variable.
#'   \item \code{G}: The generalized inverse matrix \eqn{(X'X)^-}.
#'   \item \code{XtX}: The cross-product matrix \eqn{X'X}.
#'   \item \code{design}: The design metadata list from \code{\link{build_design_glm}}.
#'   \item \code{weights}: Observation weights (if provided).
#'   \item \code{call}: The matched call.
#' }
#'
#' @references
#' Goodnight, J. H. (1979). A Tutorial on the SWEEP Operator.
#' \emph{The American Statistician}, 33(3), 149-158.
#'
#' Searle, S. R. (1971). \emph{Linear Models}. John Wiley & Sons.
#'
#' @examples
#' data(harv)
#' fit <- stat_glm(bw ~ bw0 + stock * culturetype + sex, data = harv)
#' summary(fit)
#'
#' @export
stat_glm <- function(formula, data, weights = NULL, random = NULL,
                     ginv_method = c("g2", "svd"), tol = 1e-7) {
  cl <- match.call()
  ginv_method <- match.arg(ginv_method)

  # Build design matrices and extract metadata
  d <- build_design_glm(formula, data)
  X <- d$X_ind
  y <- d$y
  n <- d$n_obs
  p <- ncol(X)

  # Process weights
  if (!is.null(weights)) {
    if (length(weights) != n || any(weights < 0) || any(is.na(weights))) {
      stop("'weights' must be a non-negative numeric vector of length equal to observations.")
    }
    w_sqrt <- sqrt(weights)
    X_w <- X * w_sqrt
    y_w <- if (is.matrix(y)) y * w_sqrt else y * w_sqrt
  } else {
    weights <- rep(1, n)
    X_w <- X
    y_w <- y
  }

  # Normal equations cross products
  XtX <- crossprod(X_w)
  Xty <- crossprod(X_w, y_w)

  # Compute generalized inverse
  if (ginv_method == "g2") {
    G <- ginv_g2(XtX, tol = tol)
  } else {
    G <- ginv_svd(XtX, tol = tol)
  }

  # Parameter estimates
  b <- G %*% Xty
  fitted_vals <- X %*% b
  resids <- y - fitted_vals

  # Matrix rank and degrees of freedom
  rk <- matrix_rank(X_w, tol = tol)
  dfe <- n - rk
  has_intercept <- d$has_intercept
  dfm <- rk - (if (has_intercept) 1L else 0L)
  dft <- n - (if (has_intercept) 1L else 0L)

  is_multivariate <- is.matrix(y) && ncol(y) > 1

  if (!is_multivariate) {
    b_vec <- as.vector(b)
    names(b_vec) <- colnames(X)
    fitted_vec <- as.vector(fitted_vals)
    names(fitted_vec) <- rownames(X)
    resid_vec <- as.vector(resids)
    names(resid_vec) <- rownames(X)

    # Weighted Sums of Squares
    w_total <- sum(weights)
    y_mean <- sum(weights * y) / w_total

    sse <- sum(weights * resid_vec^2)
    sst <- if (has_intercept) sum(weights * (y - y_mean)^2) else sum(weights * y^2)
    ssm <- sst - sse

    mse <- if (dfe > 0) sse / dfe else NA_real_
    msm <- if (dfm > 0) ssm / dfm else NA_real_
    f_val <- if (!is.na(mse) && mse > 0 && dfm > 0) msm / mse else NA_real_
    p_val <- if (!is.na(f_val)) stats::pf(f_val, dfm, dfe, lower.tail = FALSE) else NA_real_

    r_sq <- if (sst > 0) ssm / sst else NA_real_
    adj_r_sq <- if (dft > 0 && dfe > 0 && sst > 0) 1 - (sse / dfe) / (sst / dft) else NA_real_
    root_mse <- if (!is.na(mse) && mse >= 0) sqrt(mse) else NA_real_
    cv <- if (!is.na(root_mse) && abs(y_mean) > 1e-12) 100 * (root_mse / y_mean) else NA_real_

    # Redundant/aliased parameter detection
    aliased <- abs(diag(G)) < tol

    res <- list(
      coefficients = b_vec,
      residuals = resid_vec,
      fitted.values = fitted_vec,
      rank = rk,
      df.residual = dfe,
      df.model = dfm,
      df.total = dft,
      sse = sse,
      ssm = ssm,
      sst = sst,
      mse = mse,
      msm = msm,
      f_value = f_val,
      p_value = p_val,
      r_square = r_sq,
      adj_r_square = adj_r_sq,
      root_mse = root_mse,
      cv = cv,
      y_mean = y_mean,
      aliased = aliased,
      G = G,
      XtX = XtX,
      random_effects = if (!is.null(random)) parse_random_effects(random, d$term_labels) else NULL,
      design = d,
      weights = weights,
      ginv_method = ginv_method,
      call = cl
    )
    class(res) <- "stat_glm"
    return(res)
  } else {
    # Multivariate General Linear Model
    colnames(b) <- d$response_names
    rownames(b) <- colnames(X)

    # Residual SSCP matrix E
    E <- crossprod(resids * sqrt(weights))
    dimnames(E) <- list(d$response_names, d$response_names)

    # Total SSCP matrix
    Y_means <- colSums(y * weights) / sum(weights)
    Y_centered <- sweep(y, 2, Y_means, "-")
    T_sscp <- crossprod(Y_centered * sqrt(weights))

    # Model SSCP matrix H
    H <- T_sscp - E

    res <- list(
      coefficients = b,
      residuals = resids,
      fitted.values = fitted_vals,
      rank = rk,
      df.residual = dfe,
      df.model = dfm,
      df.total = dft,
      E = E,
      H = H,
      T_sscp = T_sscp,
      G = G,
      XtX = XtX,
      design = d,
      weights = weights,
      ginv_method = ginv_method,
      call = cl
    )
    class(res) <- c("stat_mlm", "stat_glm")
    return(res)
  }
}

#' @export
print.stat_glm <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  if (inherits(x, "stat_mlm")) {
    cat("Multivariate General Linear Model\n")
    cat("Number of observations:", x$design$n_obs, "\n")
    cat("Residual degrees of freedom:", x$df.residual, "\n")
    cat("Dependent variables:", paste(x$design$response_names, collapse = ", "), "\n")
  } else {
    cat("General Linear Model Fit:\n")
    cat("  Dependent Variable:", x$design$response_names, "\n")
    cat("  R-Square:          ", format(round(x$r_square, 6), nsmall = 6), "\n")
    cat("  Root MSE:          ", format(round(x$root_mse, 6), nsmall = 6), "\n")
    cat("  Coeff Var (%):     ", format(round(x$cv, 4), nsmall = 4), "\n")
    cat("  Dependent Mean:    ", format(round(x$y_mean, 6), nsmall = 6), "\n\n")

    cat("Overall Model ANOVA:\n")
    tab <- data.frame(
      DF = c(x$df.model, x$df.residual, x$df.total),
      `Sum of Squares` = c(x$ssm, x$sse, x$sst),
      `Mean Square` = c(x$msm, x$mse, NA),
      `F Value` = c(x$f_value, NA, NA),
      `Pr > F` = c(x$p_value, NA, NA),
      check.names = FALSE,
      row.names = c("Model", "Error", "Corrected Total")
    )
    print(as.data.frame(tab), digits = digits, na.print = "")
  }
  invisible(x)
}

#' @export
summary.stat_glm <- function(object, ...) {
  if (inherits(object, "stat_mlm")) {
    return(object)
  }

  se <- sqrt(pmax(0, diag(object$G) * object$mse))
  # Redundant parameters set to zero
  se[object$aliased] <- NA_real_
  t_val <- object$coefficients / se
  p_val <- 2 * stats::pt(abs(t_val), df = object$df.residual, lower.tail = FALSE)

  coef_tab <- data.frame(
    Estimate = object$coefficients,
    `Std. Error` = se,
    `t value` = t_val,
    `Pr(>|t|)` = p_val,
    check.names = FALSE
  )

  res <- list(
    call = object$call,
    fit = object,
    coefficients = coef_tab
  )
  class(res) <- "summary.stat_glm"
  res
}

#' @export
print.summary.stat_glm <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  print.stat_glm(x$fit, digits = digits, ...)
  cat("\nParameter Estimates:\n")
  cat("(Parameters set to 0 denote linearly dependent columns)\n\n")
  stats::printCoefmat(as.matrix(x$coefficients), digits = digits, na.print = "0 (biased)")
  invisible(x)
}

#' @export
coef.stat_glm <- function(object, ...) {
  object$coefficients
}

#' @export
residuals.stat_glm <- function(object, ...) {
  object$residuals
}

#' @export
fitted.stat_glm <- function(object, ...) {
  object$fitted.values
}

#' @export
vcov.stat_glm <- function(object, ...) {
  if (inherits(object, "stat_mlm")) {
    stop("vcov for multivariate models is not directly defined; use manova_glm.")
  }
  object$mse * object$G
}

#' @export
predict.stat_glm <- function(object, newdata = NULL, se.fit = FALSE,
                             interval = c("none", "confidence", "prediction"),
                             level = 0.95, ...) {
  interval <- match.arg(interval)
  if (is.null(newdata)) {
    pred <- object$fitted.values
    X_new <- object$design$X_ind
  } else {
    fake_df <- as.data.frame(newdata)
    for (resp_name in object$design$response_names) {
      fake_df[[resp_name]] <- 0
    }
    form <- stats::formula(object$design$terms)
    d_new <- build_design_glm(form, data = fake_df,
                              ref_levels = object$design$var_levels)
    X_new <- d_new$X_ind
    pred <- as.vector(X_new %*% object$coefficients)
  }

  if (se.fit || interval != "none") {
    mse <- object$mse
    XG <- as.matrix(X_new %*% object$G)
    var_mean <- pmax(0, rowSums(XG * X_new) * mse)
    se <- sqrt(var_mean)

    if (interval != "none") {
      dfe <- max(1L, as.integer(object$df.residual))
      t_crit <- stats::qt(1 - (1 - level) / 2, df = dfe)
      if (interval == "confidence") {
        lwr <- pred - t_crit * se
        upr <- pred + t_crit * se
      } else {
        se_pred <- sqrt(pmax(0, mse + var_mean))
        lwr <- pred - t_crit * se_pred
        upr <- pred + t_crit * se_pred
      }
      fit_mat <- cbind(fit = pred, lwr = lwr, upr = upr)
      if (se.fit) {
        return(list(fit = fit_mat, se.fit = se, df = dfe, residual.scale = sqrt(mse)))
      } else {
        return(fit_mat)
      }
    }

    if (se.fit) {
      return(list(fit = pred, se.fit = se, df = object$df.residual, residual.scale = sqrt(mse)))
    }
  }

  if (is.matrix(pred) && ncol(pred) == 1) as.vector(pred) else pred
}

#' @export
anova.stat_glm <- function(object, ...) {
  anova_glm(object, ...)
}
