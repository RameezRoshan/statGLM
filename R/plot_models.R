#' Model Fit Plots for statGLM Models
#'
#' @description
#' Provides modern, high-definition, publication-ready visualizations for general linear
#' models matching SAS PROC GLM specifications:
#' \itemize{
#'   \item \code{plot_fit}: Continuous regression fit with confidence (CLM) and prediction (CLI) limits.
#'   \item \code{plot_contour}: 2D continuous response surface contour plot with observed data overlay.
#'   \item \code{plot_ancova}: Analysis of covariance plots by classification levels.
#'   \item \code{plot_boxplot}: Grouped box plots with jittered observations and mean diamonds.
#'   \item \code{plot_interaction}: Two-way ANOVA interaction plots with confidence limits.
#'   \item \code{plot_nested}: Nested factor box plots.
#'   \item \code{plot.stat_glm}: S3 plot method with intelligent design auto-detection.
#' }
#'
#' @name plot_models
NULL

#' Continuous Predictor Fit Plot
#'
#' @description
#' Plots observed data, fitted response curve, confidence limits for the mean (CLM),
#' and prediction limits for individual observations (CLI) for models with a single
#' continuous predictor (matching SAS PROC GLM \code{PLOTS=FITPLOT}).
#'
#' @param object An object of class \code{"stat_glm"}.
#' @param clm Logical; show confidence limits for the mean (default \code{TRUE}).
#' @param cli Logical; show prediction limits for individuals (default \code{TRUE}).
#' @param limits Logical; if specified, sets both \code{clm} and \code{cli}.
#' @param level Confidence level (default 0.95).
#' @param n_points Number of evaluation points along the predictor grid (default 200).
#' @param theme_fn Theme function (default \code{theme_statglm}).
#' @param ... Additional arguments.
#'
#' @return A \code{ggplot} object.
#' @export
plot_fit <- function(object, clm = TRUE, cli = TRUE, limits = NULL, level = 0.95,
                     n_points = 200, theme_fn = theme_statglm, ...) {
  if (!is.null(limits)) {
    clm <- limits
    cli <- limits
  }

  d <- object$design
  mf <- d$model_frame
  num_vars <- names(d$var_types)[d$var_types == "numeric"]

  if (length(num_vars) == 0) {
    stop("plot_fit requires at least one continuous predictor.")
  }
  x_var <- num_vars[1]
  y_var <- d$response_names[1]

  x_vals <- mf[[x_var]]
  y_vals <- d$y
  x_grid <- seq(min(x_vals, na.rm = TRUE), max(x_vals, na.rm = TRUE), length.out = n_points)

  # Create newdata grid holding other variables at baseline/means
  grid_df <- data.frame(matrix(ncol = length(names(mf)), nrow = n_points))
  names(grid_df) <- names(mf)
  grid_df[[x_var]] <- x_grid

  for (v in names(mf)) {
    if (v == x_var) next
    if (v == y_var) {
      grid_df[[v]] <- 0
    } else if (d$var_types[v] == "numeric") {
      grid_df[[v]] <- mean(mf[[v]], na.rm = TRUE)
    } else {
      grid_df[[v]] <- factor(rep(levels(mf[[v]])[1], n_points), levels = levels(mf[[v]]))
    }
  }

  # Predictions and intervals
  pred_conf <- stats::predict(object, newdata = grid_df, interval = "confidence", level = level)
  pred_pred <- stats::predict(object, newdata = grid_df, interval = "prediction", level = level)

  plot_df <- data.frame(
    x = x_grid,
    fit = pred_conf[, "fit"],
    clm_lwr = pred_conf[, "lwr"],
    clm_upr = pred_conf[, "upr"],
    cli_lwr = pred_pred[, "lwr"],
    cli_upr = pred_pred[, "upr"]
  )

  obs_df <- data.frame(x = x_vals, y = as.vector(y_vals))

  p <- ggplot2::ggplot()

  # Prediction limits (CLI) - dashed boundaries
  if (cli) {
    p <- p +
      ggplot2::geom_ribbon(data = plot_df,
                           ggplot2::aes(x = .data$x, ymin = .data$cli_lwr, ymax = .data$cli_upr),
                           fill = "#e8f0fe", alpha = 0.5) +
      ggplot2::geom_line(data = plot_df, ggplot2::aes(x = .data$x, y = .data$cli_lwr),
                         linetype = "dashed", color = "#6baed6", linewidth = 0.6) +
      ggplot2::geom_line(data = plot_df, ggplot2::aes(x = .data$x, y = .data$cli_upr),
                         linetype = "dashed", color = "#6baed6", linewidth = 0.6)
  }

  # Confidence limits for mean (CLM) - solid fill
  if (clm) {
    p <- p +
      ggplot2::geom_ribbon(data = plot_df,
                           ggplot2::aes(x = .data$x, ymin = .data$clm_lwr, ymax = .data$clm_upr),
                           fill = "#b8d2eb", alpha = 0.75)
  }

  # Fitted response curve
  p <- p +
    ggplot2::geom_line(data = plot_df, ggplot2::aes(x = .data$x, y = .data$fit),
                       color = "#1a5276", linewidth = 1)

  # Observed points
  p <- p +
    ggplot2::geom_point(data = obs_df, ggplot2::aes(x = .data$x, y = .data$y),
                        color = "#2c3e50", alpha = 0.75, size = 2)

  subtitle_text <- paste0(
    "R^2 = ", round(object$r_square, 4),
    " | MSE = ", round(object$mse, 4),
    if (clm) " | Shaded: Mean CL" else "",
    if (cli) " | Dashed: Prediction Limits" else ""
  )

  p <- p +
    ggplot2::labs(
      title = paste("Fit Plot for", y_var),
      subtitle = subtitle_text,
      x = x_var,
      y = y_var
    ) +
    theme_fn()

  p
}

#' Response Surface Contour Fit Plot
#'
#' @description
#' Generates a 2D contour plot and surface heatmap of predicted response values
#' overlaid with observed data for models with two continuous predictors
#' (matching SAS PROC GLM \code{PLOTS=CONTOURFIT}).
#'
#' @param object An object of class \code{"stat_glm"}.
#' @param obs Observation display option: \code{"outlinegradient"} (default), \code{"gradient"}, \code{"outline"}, or \code{"none"}.
#' @param grid_n Grid resolution along each axis (default 50).
#' @param theme_fn Theme function (default \code{theme_statglm}).
#' @param ... Additional arguments.
#'
#' @return A \code{ggplot} object.
#' @export
plot_contour <- function(object, obs = c("outlinegradient", "gradient", "outline", "none"),
                         grid_n = 50, theme_fn = theme_statglm, ...) {
  obs <- match.arg(obs)
  d <- object$design
  mf <- d$model_frame
  num_vars <- names(d$var_types)[d$var_types == "numeric"]

  if (length(num_vars) < 2) {
    stop("plot_contour requires at least two continuous predictors.")
  }

  x1_var <- num_vars[1]
  x2_var <- num_vars[2]
  y_var <- d$response_names[1]

  x1_seq <- seq(min(mf[[x1_var]], na.rm = TRUE), max(mf[[x1_var]], na.rm = TRUE), length.out = grid_n)
  x2_seq <- seq(min(mf[[x2_var]], na.rm = TRUE), max(mf[[x2_var]], na.rm = TRUE), length.out = grid_n)
  grid_df <- expand.grid(x1 = x1_seq, x2 = x2_seq)

  eval_df <- data.frame(matrix(ncol = length(names(mf)), nrow = nrow(grid_df)))
  names(eval_df) <- names(mf)
  eval_df[[x1_var]] <- grid_df$x1
  eval_df[[x2_var]] <- grid_df$x2

  for (v in names(mf)) {
    if (v %in% c(x1_var, x2_var)) next
    if (v == y_var) {
      eval_df[[v]] <- 0
    } else if (d$var_types[v] == "numeric") {
      eval_df[[v]] <- mean(mf[[v]], na.rm = TRUE)
    } else {
      eval_df[[v]] <- factor(rep(levels(mf[[v]])[1], nrow(grid_df)), levels = levels(mf[[v]]))
    }
  }

  grid_df$fit <- stats::predict(object, newdata = eval_df)

  p <- ggplot2::ggplot() +
    # Predicted response surface (raster / tile)
    ggplot2::geom_raster(data = grid_df, ggplot2::aes(x = .data$x1, y = .data$x2, fill = .data$fit), interpolate = TRUE) +
    # Contour lines
    ggplot2::geom_contour(data = grid_df, ggplot2::aes(x = .data$x1, y = .data$x2, z = .data$fit),
                          color = "#ffffff", alpha = 0.55, linewidth = 0.4) +
    ggplot2::scale_fill_viridis_c(option = "mako", name = paste("Predicted", y_var))

  # Observed points overlay
  if (obs != "none") {
    obs_df <- data.frame(
      x1 = mf[[x1_var]],
      x2 = mf[[x2_var]],
      y = as.vector(d$y)
    )

    if (obs == "outlinegradient") {
      p <- p + ggplot2::geom_point(data = obs_df, ggplot2::aes(x = .data$x1, y = .data$x2, fill = .data$y),
                                   shape = 21, color = "#000000", stroke = 0.7, size = 2.5)
    } else if (obs == "gradient") {
      p <- p + ggplot2::geom_point(data = obs_df, ggplot2::aes(x = .data$x1, y = .data$x2, color = .data$y),
                                   size = 2.5) +
        ggplot2::scale_color_viridis_c(option = "mako", guide = "none")
    } else if (obs == "outline") {
      p <- p + ggplot2::geom_point(data = obs_df, ggplot2::aes(x = .data$x1, y = .data$x2),
                                   shape = 21, fill = "transparent", color = "#ffffff", stroke = 0.9, size = 2.5)
    }
  }

  p <- p +
    ggplot2::labs(
      title = paste("Contour Fit Plot for", y_var),
      subtitle = paste("Surface over", x1_var, "and", x2_var),
      x = x1_var,
      y = x2_var
    ) +
    theme_fn()

  p
}

#' Analysis of Covariance (ANCOVA) Plot
#'
#' @description
#' Generates an analysis of covariance plot displaying observed data points and fitted
#' regression lines with confidence bands for each classification level
#' (matching SAS PROC GLM \code{PLOTS=ANCOVAPLOT}).
#'
#' @param object An object of class \code{"stat_glm"}.
#' @param clm Logical; show confidence limits for the mean (default \code{TRUE}).
#' @param cli Logical; show prediction limits for individuals (default \code{FALSE}).
#' @param limits Logical; if specified, sets both \code{clm} and \code{cli}.
#' @param level Confidence level (default 0.95).
#' @param n_points Number of points per factor level for lines (default 100).
#' @param theme_fn Theme function (default \code{theme_statglm}).
#' @param ... Additional arguments.
#'
#' @return A \code{ggplot} object.
#' @export
plot_ancova <- function(object, clm = TRUE, cli = FALSE, limits = NULL, level = 0.95,
                        n_points = 100, theme_fn = theme_statglm, ...) {
  if (!is.null(limits)) {
    clm <- limits
    cli <- limits
  }

  d <- object$design
  mf <- d$model_frame
  fac_vars <- names(d$var_types)[d$var_types == "factor"]
  num_vars <- names(d$var_types)[d$var_types == "numeric"]

  if (length(fac_vars) == 0 || length(num_vars) == 0) {
    stop("plot_ancova requires at least one factor and one continuous covariate.")
  }

  fac_var <- fac_vars[1]
  num_var <- num_vars[1]
  y_var <- d$response_names[1]

  fac_levels <- levels(mf[[fac_var]])
  x_min <- min(mf[[num_var]], na.rm = TRUE)
  x_max <- max(mf[[num_var]], na.rm = TRUE)
  x_seq <- seq(x_min, x_max, length.out = n_points)

  # Generate prediction grid across all levels of the primary factor
  lines_df_list <- list()
  for (lev in fac_levels) {
    eval_df <- data.frame(matrix(ncol = length(names(mf)), nrow = n_points))
    names(eval_df) <- names(mf)
    eval_df[[num_var]] <- x_seq
    eval_df[[fac_var]] <- factor(rep(lev, n_points), levels = fac_levels)

    for (v in names(mf)) {
      if (v %in% c(num_var, fac_var)) next
      if (v == y_var) {
        eval_df[[v]] <- 0
      } else if (d$var_types[v] == "numeric") {
        eval_df[[v]] <- mean(mf[[v]], na.rm = TRUE)
      } else {
        eval_df[[v]] <- factor(rep(levels(mf[[v]])[1], n_points), levels = levels(mf[[v]]))
      }
    }

    pred_res <- stats::predict(object, newdata = eval_df, interval = "confidence", level = level)
    sub_df <- data.frame(
      x = x_seq,
      group = lev,
      fit = pred_res[, "fit"],
      clm_lwr = pred_res[, "lwr"],
      clm_upr = pred_res[, "upr"]
    )

    if (cli) {
      pred_cli <- stats::predict(object, newdata = eval_df, interval = "prediction", level = level)
      sub_df$cli_lwr <- pred_cli[, "lwr"]
      sub_df$cli_upr <- pred_cli[, "upr"]
    }

    lines_df_list[[lev]] <- sub_df
  }

  lines_df <- do.call(rbind, lines_df_list)
  lines_df$group <- factor(lines_df$group, levels = fac_levels)

  obs_df <- data.frame(
    x = mf[[num_var]],
    y = as.vector(d$y),
    group = mf[[fac_var]]
  )

  p <- ggplot2::ggplot()

  if (cli) {
    p <- p + ggplot2::geom_ribbon(data = lines_df,
                                 ggplot2::aes(x = .data$x, ymin = .data$cli_lwr, ymax = .data$cli_upr,
                                              group = .data$group, fill = .data$group),
                                 alpha = 0.12)
  }

  if (clm) {
    p <- p + ggplot2::geom_ribbon(data = lines_df,
                                 ggplot2::aes(x = .data$x, ymin = .data$clm_lwr, ymax = .data$clm_upr,
                                              group = .data$group, fill = .data$group),
                                 alpha = 0.22)
  }

  p <- p +
    ggplot2::geom_line(data = lines_df,
                       ggplot2::aes(x = .data$x, y = .data$fit, color = .data$group, group = .data$group),
                       linewidth = 1) +
    ggplot2::geom_point(data = obs_df,
                        ggplot2::aes(x = .data$x, y = .data$y, color = .data$group, shape = .data$group),
                        size = 2.2, alpha = 0.8) +
    ggplot2::scale_color_brewer(palette = "Set1", name = fac_var) +
    ggplot2::scale_fill_brewer(palette = "Set1", name = fac_var) +
    ggplot2::scale_shape_discrete(name = fac_var) +
    ggplot2::labs(
      title = paste("ANCOVA Plot:", y_var, "by", fac_var),
      subtitle = paste("Fitted regression by group adjusted for", num_var),
      x = num_var,
      y = y_var
    ) +
    theme_fn()

  p
}

#' Grouped Box Plot
#'
#' @description
#' Generates publication-ready box plots of response values across factor levels
#' with jittered observations and mean diamonds (matching SAS PROC GLM \code{PLOTS=BOXPLOT}).
#'
#' @param object An object of class \code{"stat_glm"}.
#' @param effect Optional character string indicating the classification effect.
#' @param jitter Logical; overlay jittered observation points (default \code{TRUE}).
#' @param mean_diamond Logical; overlay diamonds marking group sample means (default \code{TRUE}).
#' @param level Confidence level for mean error bars (default 0.95).
#' @param theme_fn Theme function (default \code{theme_statglm}).
#' @param ... Additional arguments.
#'
#' @return A \code{ggplot} object.
#' @export
plot_boxplot <- function(object, effect = NULL, jitter = TRUE, mean_diamond = TRUE,
                         level = 0.95, theme_fn = theme_statglm, ...) {
  d <- object$design
  mf <- d$model_frame
  fac_vars <- names(d$var_types)[d$var_types == "factor"]

  if (length(fac_vars) == 0) {
    stop("plot_boxplot requires at least one classification factor.")
  }

  fac_var <- if (!is.null(effect) && effect %in% fac_vars) effect else fac_vars[1]
  y_var <- d$response_names[1]

  df_plot <- data.frame(
    group = mf[[fac_var]],
    y = as.vector(d$y)
  )

  # Group summary statistics
  group_means <- stats::aggregate(y ~ group, data = df_plot, FUN = function(x) {
    n <- length(x)
    m <- mean(x)
    se <- stats::sd(x) / sqrt(n)
    tc <- stats::qt(1 - (1 - level) / 2, df = max(1, n - 1))
    c(mean = m, se = se, lwr = m - tc * se, upr = m + tc * se)
  })
  means_df <- data.frame(
    group = group_means$group,
    mean = group_means$y[, "mean"],
    lwr = group_means$y[, "lwr"],
    upr = group_means$y[, "upr"]
  )

  p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = .data$group, y = .data$y, fill = .data$group)) +
    ggplot2::geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.65, color = "#333333", linewidth = 0.6) +
    ggplot2::scale_fill_brewer(palette = "Blues", guide = "none")

  if (jitter) {
    p <- p + ggplot2::geom_jitter(width = 0.15, alpha = 0.55, size = 1.8, color = "#2c3e50")
  }

  if (mean_diamond) {
    p <- p +
      ggplot2::geom_errorbar(data = means_df,
                             ggplot2::aes(x = .data$group, ymin = .data$lwr, ymax = .data$upr),
                             inherit.aes = FALSE, width = 0.18, color = "#d62728", linewidth = 0.8) +
      ggplot2::geom_point(data = means_df,
                          ggplot2::aes(x = .data$group, y = .data$mean),
                          inherit.aes = FALSE, shape = 23, size = 3.5, fill = "#d62728", color = "#800000")
  }

  p <- p +
    ggplot2::labs(
      title = paste("Box Plot of", y_var, "by", fac_var),
      subtitle = if (mean_diamond) "Diamonds indicate group means with 95% CI" else NULL,
      x = fac_var,
      y = y_var
    ) +
    theme_fn()

  p
}

#' Two-Way Interaction Plot
#'
#' @description
#' Generates publication-ready two-way ANOVA interaction plots with cell means,
#' connecting lines, and confidence limits (matching SAS PROC GLM \code{PLOTS=INTPLOT}).
#'
#' @param object An object of class \code{"stat_glm"}.
#' @param effect Optional interaction effect name (e.g. \code{"A*B"} or \code{"A:B"}).
#' @param clm Logical; show confidence limits for the means (default \code{TRUE}).
#' @param level Confidence level (default 0.95).
#' @param theme_fn Theme function (default \code{theme_statglm}).
#' @param ... Additional arguments.
#'
#' @return A \code{ggplot} object.
#' @export
plot_interaction <- function(object, effect = NULL, clm = TRUE, level = 0.95,
                             theme_fn = theme_statglm, ...) {
  d <- object$design
  fac_vars <- names(d$var_types)[d$var_types == "factor"]

  if (length(fac_vars) < 2) {
    stop("plot_interaction requires at least two classification factors.")
  }

  if (!is.null(effect)) {
    eff_clean <- gsub("[*:]", ":", effect)
    eff_parts <- unlist(strsplit(eff_clean, ":"))
    if (length(eff_parts) >= 2 && all(eff_parts[1:2] %in% fac_vars)) {
      fac1 <- eff_parts[1]
      fac2 <- eff_parts[2]
    } else {
      fac1 <- fac_vars[1]
      fac2 <- fac_vars[2]
    }
  } else {
    fac1 <- fac_vars[1]
    fac2 <- fac_vars[2]
  }

  y_var <- d$response_names[1]

  # Compute cell LS-means for fac1:fac2
  eff_name <- paste0(fac1, ":", fac2)
  lsm <- lsmeans_glm(object, effect = eff_name, alpha = 1 - level)
  lsm_df <- lsm$lsmeans

  # Parse the grid factors
  grid <- lsm$grid
  plot_df <- data.frame(
    x = grid[[fac1]],
    group = grid[[fac2]],
    mean = lsm_df$Estimate,
    lwr = lsm_df[[7]], # Lower CL column
    upr = lsm_df[[8]]  # Upper CL column
  )

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$x, y = .data$mean,
                                             color = .data$group, group = .data$group)) +
    ggplot2::geom_line(position = ggplot2::position_dodge(0.15), linewidth = 0.9) +
    ggplot2::geom_point(position = ggplot2::position_dodge(0.15), size = 3)

  if (clm) {
    p <- p + ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data$lwr, ymax = .data$upr),
      position = ggplot2::position_dodge(0.15), width = 0.18, linewidth = 0.7
    )
  }

  p <- p +
    ggplot2::scale_color_brewer(palette = "Set1", name = fac2) +
    ggplot2::labs(
      title = paste("Interaction Plot:", y_var),
      subtitle = paste("LS-Means for", fac1, "by", fac2, if (clm) "with 95% CI" else ""),
      x = fac1,
      y = paste("LS-Mean of", y_var)
    ) +
    theme_fn()

  p
}

#' Nested Classification Box Plot
#'
#' @description
#' Generates clustered box plots for nested classification designs
#' (matching SAS PROC GLM \code{PLOTS=NESTPLOT}).
#'
#' @param object An object of class \code{"stat_glm"}.
#' @param effect Optional effect specification.
#' @param boxwidth Box width (default 0.7).
#' @param clusterwidth Cluster width (default 0.8).
#' @param theme_fn Theme function (default \code{theme_statglm}).
#' @param ... Additional arguments.
#'
#' @return A \code{ggplot} object.
#' @export
plot_nested <- function(object, effect = NULL, boxwidth = 0.7, clusterwidth = 0.8,
                        theme_fn = theme_statglm, ...) {
  d <- object$design
  mf <- d$model_frame
  fac_vars <- names(d$var_types)[d$var_types == "factor"]

  if (length(fac_vars) < 2) {
    stop("plot_nested requires at least two classification factors.")
  }

  parent_fac <- fac_vars[1]
  nested_fac <- fac_vars[2]
  y_var <- d$response_names[1]

  df_plot <- data.frame(
    parent = mf[[parent_fac]],
    nested = mf[[nested_fac]],
    y = as.vector(d$y)
  )

  p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = .data$parent, y = .data$y, fill = .data$nested)) +
    ggplot2::geom_boxplot(position = ggplot2::position_dodge(clusterwidth), width = boxwidth,
                          alpha = 0.75, color = "#333333", linewidth = 0.5) +
    ggplot2::scale_fill_brewer(palette = "Blues", name = nested_fac) +
    ggplot2::labs(
      title = paste("Nested Box Plot of", y_var),
      subtitle = paste(nested_fac, "nested within", parent_fac),
      x = parent_fac,
      y = y_var
    ) +
    theme_fn()

  p
}

#' S3 Plot Method for statGLM Models
#'
#' @description
#' Master plotting interface for \code{"stat_glm"} objects. Automatically detects the
#' design structure of the model and displays the appropriate SAS PROC GLM graphic,
#' or produces specific requested plots via the \code{which} argument.
#'
#' @param x An object of class \code{"stat_glm"}.
#' @param which Specific plot to produce: \code{"auto"} (default), \code{"diagnostics"},
#'   \code{"residuals"}, \code{"fitplot"}, \code{"contourfit"}, \code{"ancovaplot"},
#'   \code{"boxplot"}, \code{"intplot"}, \code{"nestplot"}, or \code{"all"}.
#' @param ... Arguments passed to the specific plotting function.
#'
#' @return A \code{ggplot} object, a grob arrangement, or a list of plots.
#'
#' @examples
#' data(harv)
#' fit <- stat_glm(bw ~ bw0 + stock * culturetype + sex, data = harv)
#' plot(fit) # Automatically produces interaction plot
#' plot(fit, which = "diagnostics")
#'
#' @export
plot.stat_glm <- function(x, which = c("auto", "diagnostics", "residuals", "fitplot",
                                       "contourfit", "ancovaplot", "boxplot",
                                       "intplot", "nestplot", "all"), ...) {
  if (is.null(which)) which <- "auto"
  which <- match.arg(which)

  d <- x$design
  fac_vars <- names(d$var_types)[d$var_types == "factor"]
  num_vars <- names(d$var_types)[d$var_types == "numeric"]
  n_fac <- length(fac_vars)
  n_num <- length(num_vars)

  if (which == "auto") {
    # Check for nested model terms
    has_nested <- any(grepl("[:%in%]", d$term_labels) & !grepl("[*]", d$term_labels)) && n_fac >= 2

    if (n_fac == 1 && n_num == 0) {
      which <- "boxplot"
    } else if (n_fac == 2 && n_num == 0) {
      which <- if (has_nested) "nestplot" else "intplot"
    } else if (n_fac == 0 && n_num == 1) {
      which <- "fitplot"
    } else if (n_fac == 0 && n_num == 2) {
      which <- "contourfit"
    } else if (n_fac >= 1 && n_num >= 1) {
      which <- "ancovaplot"
    } else {
      which <- "diagnostics"
    }
  }

  if (which == "all") {
    res_list <- list()
    res_list$diagnostics <- plot_diagnostics(x, unpack = TRUE, ...)
    if (n_num > 0) res_list$residuals <- plot_residuals(x, unpack = TRUE, ...)
    if (n_num == 1 && n_fac == 0) res_list$fitplot <- plot_fit(x, ...)
    if (n_num == 2 && n_fac == 0) res_list$contourfit <- plot_contour(x, ...)
    if (n_fac >= 1 && n_num >= 1) res_list$ancovaplot <- plot_ancova(x, ...)
    if (n_fac >= 1) res_list$boxplot <- plot_boxplot(x, ...)
    if (n_fac >= 2) res_list$intplot <- plot_interaction(x, ...)
    return(res_list)
  }

  switch(
    which,
    diagnostics = plot_diagnostics(x, ...),
    residuals = plot_residuals(x, ...),
    fitplot = plot_fit(x, ...),
    contourfit = plot_contour(x, ...),
    ancovaplot = plot_ancova(x, ...),
    boxplot = plot_boxplot(x, ...),
    intplot = plot_interaction(x, ...),
    nestplot = plot_nested(x, ...),
    plot_diagnostics(x, ...)
  )
}

#' General Plotting Wrapper for statGLM
#'
#' @description
#' Explicit function interface matching \code{plot.stat_glm}.
#'
#' @param object An object of class \code{"stat_glm"}.
#' @param which Plot type to produce.
#' @param ... Arguments passed to \code{plot.stat_glm}.
#'
#' @export
plot_glm <- function(object, which = "auto", ...) {
  plot.stat_glm(object, which = which, ...)
}
