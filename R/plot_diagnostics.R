#' Compute Diagnostics for statGLM Models
#'
#' @description
#' Calculates leverage (hat diagonals), internally and externally studentized residuals
#' (RStudent), Cook's distance, and outlier flags for a \code{"stat_glm"} model fit.
#'
#' @param object An object of class \code{"stat_glm"}.
#' @param id_var Optional vector of observation identifiers.
#'
#' @return A data frame containing observation-level diagnostics.
#' @keywords internal
compute_glm_diagnostics <- function(object, id_var = NULL) {
  if (!inherits(object, "stat_glm")) {
    stop("'object' must be an object of class 'stat_glm'.")
  }

  y <- object$design$y
  if (is.matrix(y) && ncol(y) > 1) {
    stop("Diagnostic plots are currently supported for univariate response models.")
  }

  y <- as.vector(y)
  y_hat <- as.vector(object$fitted.values)
  e <- as.vector(object$residuals)
  n <- length(y)
  p <- max(1L, as.integer(object$rank))
  dfe <- max(1L, as.integer(object$df.residual))
  mse <- object$mse
  if (is.na(mse) || mse <= 0) mse <- sum(e^2) / dfe

  w <- if (!is.null(object$weights)) object$weights else rep(1, n)
  X <- object$design$X_ind
  G <- object$G

  # Leverage h_ii = w_i * sum_j (X G)_ij X_ij
  XG <- as.matrix(X %*% G)
  h <- w * rowSums(XG * X)
  h <- pmax(0, pmin(1 - 1e-10, h)) # Numerical safeguard

  # Internally studentized residuals
  se_e <- sqrt(pmax(1e-12, mse * (1 - h) / w))
  r_intern <- e / se_e

  # Externally studentized residuals (RStudent / jackknife)
  denom_rstudent <- dfe - r_intern^2
  r_extern <- ifelse(denom_rstudent > 1e-10,
                     r_intern * sqrt((dfe - 1) / denom_rstudent),
                     r_intern)

  # Cook's Distance
  cooksd <- (r_intern^2 / p) * (h / (1 - h))

  # IDs / Labels
  if (is.null(id_var)) {
    obs_id <- if (!is.null(names(e))) names(e) else as.character(seq_len(n))
  } else {
    obs_id <- as.character(id_var)
  }

  # Outlier and influence cutoffs (matching SAS PROC GLM)
  h_cutoff <- 2 * p / n
  cook_cutoff <- 4 / n
  r_cutoff <- 2

  outlier <- abs(r_extern) > r_cutoff
  high_leverage <- h > h_cutoff
  influential <- cooksd > cook_cutoff

  data.frame(
    obs = seq_len(n),
    id = obs_id,
    y = y,
    fitted = y_hat,
    residual = e,
    leverage = h,
    r_intern = r_intern,
    r_student = r_extern,
    cooksd = cooksd,
    outlier = outlier,
    high_leverage = high_leverage,
    influential = influential,
    label_flag = outlier | high_leverage | influential,
    stringsAsFactors = FALSE
  )
}

#' Diagnostic Plots for statGLM Models
#'
#' @description
#' Generates a comprehensive, modern, publication-ready panel of model diagnostics
#' matching the \code{PLOTS=DIAGNOSTICS} specifications of SAS PROC GLM.
#'
#' Includes:
#' \enumerate{
#'   \item Residuals vs. Predicted Values (with LOESS trend)
#'   \item Studentized Residuals (RStudent) vs. Predicted Values
#'   \item Observed vs. Predicted Values (with 1:1 identity line)
#'   \item Studentized Residuals vs. Leverage (with cutoffs)
#'   \item Cook's Distance by Observation (needle plot with 4/n cutoff)
#'   \item Normal Q-Q Plot of Residuals
#'   \item Residual Histogram with Normal and Kernel Densities
#'   \item Residual-Fit (RF) Spread Plot
#' }
#'
#' @param object An object of class \code{"stat_glm"}.
#' @param label Logical; if \code{TRUE}, labels points exceeding outlier/influence cutoffs. Default is \code{TRUE}.
#' @param unpack Logical; if \code{TRUE}, returns a named list of individual ggplot objects instead of a single paneled grid. Default is \code{FALSE}.
#' @param id_var Optional variable or vector to use for point labels.
#' @param theme_fn Theme function to use. Defaults to \code{theme_statglm}.
#' @param ... Additional arguments.
#'
#' @return A grob/gtable multi-panel arrangement (if \code{unpack = FALSE}) or a list of \code{ggplot} objects (if \code{unpack = TRUE}).
#'
#' @examples
#' data(warpbreaks)
#' fit <- stat_glm(breaks ~ wool * tension, data = warpbreaks)
#' p_diag <- plot_diagnostics(fit)
#'
#' @export
plot_diagnostics <- function(object, label = TRUE, unpack = FALSE, id_var = NULL,
                             theme_fn = theme_statglm, ...) {
  diag_df <- compute_glm_diagnostics(object, id_var = id_var)
  n <- nrow(diag_df)
  p <- max(1L, as.integer(object$rank))
  dfe <- max(1L, as.integer(object$df.residual))
  h_cut <- 2 * p / n
  cook_cut <- 4 / n

  color_point <- "#2b5c8f"
  color_line <- "#d62728"
  color_ref <- "#888888"

  # Helper for labeling
  add_labels <- function(p_obj, subset_cond) {
    if (!label) return(p_obj)
    df_sub <- diag_df[subset_cond, , drop = FALSE]
    if (nrow(df_sub) == 0) return(p_obj)
    p_obj + ggplot2::geom_text(
      data = df_sub,
      ggplot2::aes(label = .data$id),
      size = 2.8, vjust = -0.7, color = "#990000", fontface = "bold"
    )
  }

  # 1. Residuals vs Predicted
  p1 <- ggplot2::ggplot(diag_df, ggplot2::aes(x = .data$fitted, y = .data$residual)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = color_ref, linewidth = 0.5) +
    ggplot2::geom_point(color = color_point, alpha = 0.75, size = 1.8) +
    ggplot2::geom_smooth(method = "loess", formula = y ~ x, se = FALSE,
                         color = color_line, linewidth = 0.8) +
    ggplot2::labs(title = "Residual vs Predicted", x = "Predicted Value", y = "Residual") +
    theme_fn(base_size = 9)

  # 2. Studentized Residuals vs Predicted
  p2 <- ggplot2::ggplot(diag_df, ggplot2::aes(x = .data$fitted, y = .data$r_student)) +
    ggplot2::geom_hline(yintercept = c(-2, 0, 2), linetype = c("dotted", "dashed", "dotted"),
                        color = c(color_line, color_ref, color_line), linewidth = 0.5) +
    ggplot2::geom_point(ggplot2::aes(color = .data$outlier), alpha = 0.8, size = 1.8) +
    ggplot2::scale_color_manual(values = c("FALSE" = color_point, "TRUE" = color_line), guide = "none") +
    ggplot2::labs(title = "RStudent vs Predicted", x = "Predicted Value", y = "Studentized Residual") +
    theme_fn(base_size = 9)
  p2 <- add_labels(p2, diag_df$outlier)

  # 3. Observed vs Predicted
  p3 <- ggplot2::ggplot(diag_df, ggplot2::aes(x = .data$fitted, y = .data$y)) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = color_ref, linewidth = 0.5) +
    ggplot2::geom_point(color = color_point, alpha = 0.75, size = 1.8) +
    ggplot2::labs(title = "Observed vs Predicted", x = "Predicted Value", y = "Observed Value") +
    theme_fn(base_size = 9)

  # 4. Studentized Residuals vs Leverage
  p4 <- ggplot2::ggplot(diag_df, ggplot2::aes(x = .data$leverage, y = .data$r_student)) +
    ggplot2::geom_hline(yintercept = c(-2, 0, 2), linetype = c("dotted", "dashed", "dotted"),
                        color = c(color_line, color_ref, color_line), linewidth = 0.5) +
    ggplot2::geom_vline(xintercept = h_cut, linetype = "dotted", color = color_line, linewidth = 0.5) +
    ggplot2::geom_point(ggplot2::aes(color = .data$high_leverage | .data$outlier), alpha = 0.8, size = 1.8) +
    ggplot2::scale_color_manual(values = c("FALSE" = color_point, "TRUE" = color_line), guide = "none") +
    ggplot2::labs(title = "RStudent vs Leverage", x = paste0("Leverage (cutoff = ", round(h_cut, 3), ")"),
                  y = "Studentized Residual") +
    theme_fn(base_size = 9)
  p4 <- add_labels(p4, diag_df$high_leverage | diag_df$outlier)

  # 5. Cook's Distance Needle Plot
  p5 <- ggplot2::ggplot(diag_df, ggplot2::aes(x = .data$obs, y = .data$cooksd)) +
    ggplot2::geom_hline(yintercept = cook_cut, linetype = "dashed", color = color_line, linewidth = 0.5) +
    ggplot2::geom_segment(ggplot2::aes(x = .data$obs, xend = .data$obs, y = 0, yend = .data$cooksd),
                          color = color_point, linewidth = 0.6) +
    ggplot2::geom_point(ggplot2::aes(color = .data$influential), size = 1.6) +
    ggplot2::scale_color_manual(values = c("FALSE" = color_point, "TRUE" = color_line), guide = "none") +
    ggplot2::labs(title = paste0("Cook's D (cutoff = ", round(cook_cut, 3), ")"),
                  x = "Observation Number", y = "Cook's D") +
    theme_fn(base_size = 9)
  p5 <- add_labels(p5, diag_df$influential)

  # 6. Normal Q-Q Plot
  sorted_r <- sort(diag_df$r_student)
  probs <- (seq_len(n) - 0.375) / (n + 0.25)
  theo_q <- stats::qnorm(probs)
  qq_df <- data.frame(theo = theo_q, sample = sorted_r)
  
  # Quantile reference line passing through 25th and 75th percentiles
  q25_x <- stats::qnorm(0.25); q75_x <- stats::qnorm(0.75)
  q25_y <- stats::quantile(sorted_r, 0.25); q75_y <- stats::quantile(sorted_r, 0.75)
  slope <- (q75_y - q25_y) / (q75_x - q25_x)
  intercept <- q25_y - slope * q25_x

  p6 <- ggplot2::ggplot(qq_df, ggplot2::aes(x = .data$theo, y = .data$sample)) +
    ggplot2::geom_abline(intercept = intercept, slope = slope, linetype = "dashed",
                         color = color_line, linewidth = 0.6) +
    ggplot2::geom_point(color = color_point, alpha = 0.8, size = 1.8) +
    ggplot2::labs(title = "Normal Q-Q Plot", x = "Theoretical Quantiles", y = "Sample Residual Quantiles") +
    theme_fn(base_size = 9)

  # 7. Residual Histogram with Densities
  s_sd <- stats::sd(diag_df$residual)
  p7 <- ggplot2::ggplot(diag_df, ggplot2::aes(x = .data$residual)) +
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)),
                            bins = max(8, min(25, round(sqrt(n)))),
                            fill = "#dce7f5", color = "#4682b4", linewidth = 0.4) +
    ggplot2::geom_density(color = color_line, linewidth = 0.8) +
    ggplot2::stat_function(fun = stats::dnorm, args = list(mean = 0, sd = s_sd),
                           linetype = "dashed", color = color_ref, linewidth = 0.7) +
    ggplot2::labs(title = "Residual Distribution", x = "Residual", y = "Density") +
    theme_fn(base_size = 9)

  # 8. RF (Residual-Fit) Spread Plot
  fit_centered <- sort(diag_df$fitted - mean(diag_df$fitted))
  res_sorted <- sort(diag_df$residual)
  f_vals <- (seq_len(n) - 0.5) / n
  rf_df <- rbind(
    data.frame(f = f_vals, value = fit_centered, component = "Centered Fit"),
    data.frame(f = f_vals, value = res_sorted, component = "Residuals")
  )
  p8 <- ggplot2::ggplot(rf_df, ggplot2::aes(x = .data$f, y = .data$value, color = .data$component)) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::scale_color_manual(values = c("Centered Fit" = color_point, "Residuals" = color_line)) +
    ggplot2::labs(title = "RF-Spread Plot", x = "Fraction of Data", y = "Spread", color = NULL) +
    theme_fn(base_size = 9) +
    ggplot2::theme(legend.position = "top", legend.margin = ggplot2::margin(0, 0, 0, 0))

  plot_list <- list(
    ResidualByPredicted = p1,
    RStudentByPredicted = p2,
    ObservedByPredicted = p3,
    RStudentByLeverage = p4,
    CooksDPlot = p5,
    QQPlot = p6,
    ResidualHistogram = p7,
    RFPlot = p8
  )

  if (unpack) {
    return(plot_list)
  }

  # Multi-panel publication display (2 rows x 4 columns)
  title_grob <- grid::textGrob(
    paste("Fit Diagnostics for", deparse(object$call$formula)),
    gp = grid::gpar(fontface = "bold", fontsize = 12)
  )
  gridExtra::arrangeGrob(
    p1, p2, p3, p4,
    p5, p6, p7, p8,
    ncol = 4, nrow = 2,
    top = title_grob
  )
}

#' Plot Residuals Against Continuous Covariates
#'
#' @description
#' Generates scatter plots of model residuals against each continuous covariate,
#' optionally overlaying a LOESS smoother matching SAS PROC GLM \code{PLOTS=RESIDUALS}.
#'
#' @param object An object of class \code{"stat_glm"}.
#' @param smooth Logical; if \code{TRUE}, adds a LOESS smoother with confidence ribbon. Default is \code{TRUE}.
#' @param unpack Logical; if \code{TRUE}, returns a list of individual ggplot objects. Default is \code{FALSE}.
#' @param theme_fn Theme function to use. Defaults to \code{theme_statglm}.
#' @param ... Additional arguments.
#'
#' @return A grob arrangement or list of \code{ggplot} objects.
#'
#' @export
plot_residuals <- function(object, smooth = TRUE, unpack = FALSE,
                           theme_fn = theme_statglm, ...) {
  if (!inherits(object, "stat_glm")) {
    stop("'object' must be an object of class 'stat_glm'.")
  }

  d <- object$design
  mf <- d$model_frame
  resids <- as.vector(object$residuals)
  num_vars <- names(d$var_types)[d$var_types == "numeric"]

  if (length(num_vars) == 0) {
    message("Model contains no continuous covariates for residual covariate plotting.")
    return(invisible(NULL))
  }

  plots <- list()
  for (nv in num_vars) {
    df_plot <- data.frame(x = mf[[nv]], residual = resids)
    p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = .data$x, y = .data$residual)) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "#888888", linewidth = 0.5) +
      ggplot2::geom_point(color = "#2b5c8f", alpha = 0.75, size = 2) +
      ggplot2::labs(title = paste("Residuals vs", nv), x = nv, y = "Residual") +
      theme_fn()

    if (smooth) {
      p <- p + ggplot2::geom_smooth(method = "loess", formula = y ~ x, se = TRUE,
                                   color = "#d62728", fill = "#f5dada", linewidth = 0.8)
    }
    plots[[nv]] <- p
  }

  if (unpack || length(plots) == 1) {
    if (length(plots) == 1 && !unpack) return(plots[[1]])
    return(plots)
  }

  gridExtra::arrangeGrob(grobs = plots, ncol = min(3, length(plots)))
}
