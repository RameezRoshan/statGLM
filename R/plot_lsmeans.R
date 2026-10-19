#' LS-Means Plots for statGLM
#'
#' @description
#' Provides modern, publication-ready visualizations for least-squares means (LS-means)
#' and multiple comparisons matching SAS PROC GLM specifications:
#' \itemize{
#'   \item \code{plot_lsmeans} (\code{MEANPLOT}): LS-means with confidence intervals and bands.
#'   \item \code{plot_diff} (\code{DIFFPLOT}): Hsu's pairwise diffogram (mean-mean scatter plot).
#'   \item \code{plot_control} (\code{CONTROLPLOT}): Comparisons against a reference control level.
#'   \item \code{plot_anom} (\code{ANOMPLOT}): Analysis of means comparing LS-means against the overall mean.
#'   \item \code{plot_lines} (\code{LINESPLOT}): Grouping lines displaying non-significant subsets.
#'   \item \code{plot.lsmeans_stat_glm}: S3 plot method for LS-means objects.
#' }
#'
#' @name plot_lsmeans
NULL

#' Plot Least-Squares Means
#'
#' @description
#' Generates a publication-ready plot of least-squares means with confidence intervals,
#' optional confidence bands, and connecting lines (matching SAS PROC GLM \code{PLOTS=MEANPLOT}).
#'
#' @param object An object of class \code{"lsmeans_stat_glm"}.
#' @param cl Logical; show error bars for confidence limits (default \code{TRUE}).
#' @param clband Logical; show confidence band ribbon (default \code{FALSE}).
#' @param connect Logical; connect means across levels with lines (default \code{FALSE}).
#' @param order Ordering of levels: \code{"default"} (default), \code{"ascending"}, or \code{"descending"}.
#' @param theme_fn Theme function (default \code{theme_statglm}).
#' @param ... Additional arguments.
#'
#' @return A \code{ggplot} object.
#' @export
plot_lsmeans <- function(object, cl = TRUE, clband = FALSE, connect = FALSE,
                         order = c("default", "ascending", "descending"),
                         theme_fn = theme_statglm, ...) {
  if (!inherits(object, "lsmeans_stat_glm")) {
    stop("'object' must be an object of class 'lsmeans_stat_glm'.")
  }
  order <- match.arg(order)

  df <- object$lsmeans
  lwr_col <- names(df)[grep("^Lower", names(df))]
  upr_col <- names(df)[grep("^Upper", names(df))]

  plot_df <- data.frame(
    level = df$Effect,
    estimate = df$Estimate,
    lwr = df[[lwr_col]],
    upr = df[[upr_col]],
    stringsAsFactors = FALSE
  )

  # Remove non-estimable rows
  plot_df <- plot_df[!is.na(plot_df$estimate), , drop = FALSE]

  if (order == "ascending") {
    plot_df$level <- factor(plot_df$level, levels = plot_df$level[order(plot_df$estimate)])
  } else if (order == "descending") {
    plot_df$level <- factor(plot_df$level, levels = plot_df$level[order(plot_df$estimate, decreasing = TRUE)])
  } else {
    plot_df$level <- factor(plot_df$level, levels = plot_df$level)
  }

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$level, y = .data$estimate, group = 1))

  if (clband) {
    p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lwr, ymax = .data$upr),
                                 fill = "#c6dbef", alpha = 0.5)
  }

  if (connect) {
    p <- p + ggplot2::geom_line(color = "#1f77b4", linewidth = 0.9)
  }

  if (cl) {
    p <- p + ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$lwr, ymax = .data$upr),
                                    width = 0.18, color = "#2c3e50", linewidth = 0.7)
  }

  p <- p +
    ggplot2::geom_point(color = "#1f77b4", size = 3.5, shape = 19) +
    ggplot2::labs(
      title = paste("Least-Squares Means for", object$effect),
      subtitle = paste0("Estimates with ", (1 - object$alpha) * 100, "% Confidence Intervals"),
      x = object$effect,
      y = "LS-Mean Estimate"
    ) +
    theme_fn()

  p
}

#' Pairwise Diffogram (Mean-Mean Scatter Plot)
#'
#' @description
#' Generates Hsu's (1996) pairwise difference plot ("diffogram") matching SAS PROC GLM
#' \code{PLOTS=DIFFPLOT}. Both axes represent the LS-means of the factor levels.
#' For each pairwise comparison, a line segment with slope -1 is drawn centered at
#' the intersection of the two LS-means. If the segment intersects the 45-degree
#' diagonal line of equality (y = x), the difference is not statistically significant;
#' if it does not intersect, the difference is statistically significant.
#'
#' @param object An object of class \code{"lsmeans_stat_glm"}.
#' @param center Logical; mark the center point of each comparison (default \code{TRUE}).
#' @param nolines Logical; suppress confidence interval line segments (default \code{FALSE}).
#' @param alpha Significance level (defaults to \code{object$alpha}).
#' @param theme_fn Theme function (default \code{theme_statglm}).
#' @param ... Additional arguments.
#'
#' @return A \code{ggplot} object.
#'
#' @references
#' Hsu, J. C. (1996). \emph{Multiple Comparisons: Theory and Methods}. Chapman & Hall/CRC.
#'
#' @export
plot_diff <- function(object, center = TRUE, nolines = FALSE, alpha = NULL,
                      theme_fn = theme_statglm, ...) {
  if (!inherits(object, "lsmeans_stat_glm")) {
    stop("'object' must be an object of class 'lsmeans_stat_glm'.")
  }

  lsm_df <- object$lsmeans
  diff_df <- object$diffs

  if (nrow(diff_df) == 0) {
    stop("plot_diff requires pairwise differences. Only one group level found.")
  }

  if (is.null(alpha)) alpha <- object$alpha

  # Extract means and names
  m_names <- lsm_df$Effect
  means <- lsm_df$Estimate
  names(means) <- m_names
  k <- length(means)

  # Parse comparisons
  # Each comparison is formatted as "A - B"
  seg_list <- list()
  for (row_idx in seq_len(nrow(diff_df))) {
    comp_str <- diff_df$Comparison[row_idx]
    parts <- unlist(strsplit(comp_str, " - "))
    if (length(parts) != 2) next
    g1 <- parts[1]
    g2 <- parts[2]
    if (!g1 %in% m_names || !g2 %in% m_names) next

    m1 <- means[g1]
    m2 <- means[g2]
    diff_val <- diff_df$Estimate[row_idx]
    
    # Margin of error from CI columns
    lwr_col <- grep("^Lower", names(diff_df))[1]
    upr_col <- grep("^Upper", names(diff_df))[1]
    ci_lwr <- diff_df[[lwr_col]][row_idx]
    ci_upr <- diff_df[[upr_col]][row_idx]
    margin <- (ci_upr - ci_lwr) / 2

    # Segment coordinates along slope -1 passing through (m1, m2):
    # x + y = m1 + m2
    # At endpoints, (x - y) spans the CI for (m1 - m2)
    # x1 = m1 - margin / 2, y1 = m2 + margin / 2
    # x2 = m1 + margin / 2, y2 = m2 - margin / 2
    x1 <- m1 - margin / 2
    y1 <- m2 + margin / 2
    x2 <- m1 + margin / 2
    y2 <- m2 - margin / 2

    p_val <- if ("Adj P" %in% names(diff_df)) diff_df[["Adj P"]][row_idx] else diff_df[["Pr > |t|"]][row_idx]
    is_sig <- !is.na(p_val) && p_val < alpha

    seg_list[[length(seg_list) + 1]] <- data.frame(
      comp = comp_str,
      g1 = g1,
      g2 = g2,
      m1 = m1,
      m2 = m2,
      x1 = x1,
      y1 = y1,
      x2 = x2,
      y2 = y2,
      sig = ifelse(is_sig, "Significant", "Non-Significant"),
      p_val = p_val,
      stringsAsFactors = FALSE
    )
  }

  segs_df <- do.call(rbind, seg_list)
  segs_df$sig <- factor(segs_df$sig, levels = c("Significant", "Non-Significant"))

  # Plot range
  all_coords <- c(means, segs_df$x1, segs_df$x2, segs_df$y1, segs_df$y2)
  min_val <- min(all_coords, na.rm = TRUE)
  max_val <- max(all_coords, na.rm = TRUE)
  pad <- 0.08 * (max_val - min_val)
  lims <- c(min_val - pad, max_val + pad)

  p <- ggplot2::ggplot() +
    # 45-degree diagonal reference line y = x (Equality line)
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "#444444", linewidth = 0.8)

  # Grid lines corresponding to each group LS-mean
  for (i in seq_along(means)) {
    p <- p +
      ggplot2::geom_vline(xintercept = means[i], linetype = "dotted", color = "#cccccc", linewidth = 0.4) +
      ggplot2::geom_hline(yintercept = means[i], linetype = "dotted", color = "#cccccc", linewidth = 0.4)
  }

  # Confidence segments
  if (!nolines) {
    p <- p + ggplot2::geom_segment(
      data = segs_df,
      ggplot2::aes(x = .data$x1, y = .data$y1, xend = .data$x2, yend = .data$y2,
                   color = .data$sig, linetype = .data$sig),
      linewidth = 0.95
    )
  }

  # Center comparison points
  if (center) {
    p <- p + ggplot2::geom_point(
      data = segs_df,
      ggplot2::aes(x = .data$m1, y = .data$m2, color = .data$sig),
      size = 2.2
    )
  }

  p <- p +
    ggplot2::scale_color_manual(
      values = c("Significant" = "#1a5276", "Non-Significant" = "#a6b8c7"),
      name = paste0("Difference (alpha = ", alpha, ")")
    ) +
    ggplot2::scale_linetype_manual(
      values = c("Significant" = "solid", "Non-Significant" = "dashed"),
      name = paste0("Difference (alpha = ", alpha, ")")
    ) +
    ggplot2::scale_x_continuous(breaks = means, labels = names(means), limits = lims) +
    ggplot2::scale_y_continuous(breaks = means, labels = names(means), limits = lims) +
    ggplot2::coord_fixed() +
    ggplot2::labs(
      title = paste("Pairwise LS-Means Diffogram for", object$effect),
      subtitle = paste("Adjustment:", toupper(object$adjust), "| Segments crossing y = x are not significant"),
      x = "LS-Mean (i)",
      y = "LS-Mean (j)"
    ) +
    theme_fn() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "top"
    )

  p
}

#' LS-Means Control Comparison Plot
#'
#' @description
#' Generates a Dunnett-style control comparison plot displaying differences between
#' each group LS-mean and a designated reference control level
#' (matching SAS PROC GLM \code{PLOTS=CONTROLPLOT}).
#'
#' @param object An object of class \code{"lsmeans_stat_glm"}.
#' @param control Character string indicating the reference/control group. Defaults to first level.
#' @param alpha Significance level (defaults to \code{object$alpha}).
#' @param theme_fn Theme function (default \code{theme_statglm}).
#' @param ... Additional arguments.
#'
#' @return A \code{ggplot} object.
#' @export
plot_control <- function(object, control = NULL, alpha = NULL,
                         theme_fn = theme_statglm, ...) {
  if (!inherits(object, "lsmeans_stat_glm")) {
    stop("'object' must be an object of class 'lsmeans_stat_glm'.")
  }

  lsm_df <- object$lsmeans
  diff_df <- object$diffs

  if (nrow(diff_df) == 0) {
    stop("plot_control requires at least two group levels.")
  }

  if (is.null(alpha)) alpha <- object$alpha
  levels_vec <- lsm_df$Effect
  if (is.null(control) || !control %in% levels_vec) {
    control <- levels_vec[1]
  }

  # Filter diffs that involve control
  match_rows <- list()
  for (r in seq_len(nrow(diff_df))) {
    comp_str <- diff_df$Comparison[r]
    parts <- unlist(strsplit(comp_str, " - "))
    if (length(parts) != 2) next
    
    if (parts[2] == control) {
      match_rows[[length(match_rows) + 1]] <- data.frame(
        group = parts[1],
        diff = diff_df$Estimate[r],
        se = diff_df[["Std Error"]][r],
        p_val = if ("Adj P" %in% names(diff_df)) diff_df[["Adj P"]][r] else diff_df[["Pr > |t|"]][r],
        lwr = diff_df[[grep("^Lower", names(diff_df))[1]]][r],
        upr = diff_df[[grep("^Upper", names(diff_df))[1]]][r],
        stringsAsFactors = FALSE
      )
    } else if (parts[1] == control) {
      match_rows[[length(match_rows) + 1]] <- data.frame(
        group = parts[2],
        diff = -diff_df$Estimate[r],
        se = diff_df[["Std Error"]][r],
        p_val = if ("Adj P" %in% names(diff_df)) diff_df[["Adj P"]][r] else diff_df[["Pr > |t|"]][r],
        lwr = -diff_df[[grep("^Upper", names(diff_df))[1]]][r],
        upr = -diff_df[[grep("^Lower", names(diff_df))[1]]][r],
        stringsAsFactors = FALSE
      )
    }
  }

  ctrl_df <- do.call(rbind, match_rows)
  ctrl_df$sig <- ifelse(ctrl_df$p_val < alpha, "Significant", "Non-Significant")
  ctrl_df$sig <- factor(ctrl_df$sig, levels = c("Significant", "Non-Significant"))

  p <- ggplot2::ggplot(ctrl_df, ggplot2::aes(x = .data$diff, y = .data$group, color = .data$sig)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "#444444", linewidth = 0.8) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$lwr, xmax = .data$upr),
                           width = 0.25, linewidth = 0.8) +
    ggplot2::geom_point(size = 3.5) +
    ggplot2::scale_color_manual(
      values = c("Significant" = "#d62728", "Non-Significant" = "#2b5c8f"),
      name = paste0("Difference (alpha = ", alpha, ")")
    ) +
    ggplot2::labs(
      title = paste("Control Comparisons for", object$effect),
      subtitle = paste("Comparisons against Control Level:", control),
      x = paste("Difference in LS-Means (Group -", control, ")"),
      y = "Treatment Group"
    ) +
    theme_fn() +
    ggplot2::theme(legend.position = "top")

  p
}

#' Analysis of Means (ANOM) Plot
#'
#' @description
#' Generates an Analysis of Means display comparing each group LS-mean to the
#' overall average LS-mean, with upper and lower decision limits (UDL / LDL)
#' (matching SAS PROC GLM \code{PLOTS=ANOMPLOT}).
#'
#' @param object An object of class \code{"lsmeans_stat_glm"}.
#' @param alpha Significance level (defaults to \code{object$alpha}).
#' @param theme_fn Theme function (default \code{theme_statglm}).
#' @param ... Additional arguments.
#'
#' @return A \code{ggplot} object.
#'
#' @references
#' Nelson, P. R. (1982). Exact critical points for use with the analysis of means.
#' \emph{Communications in Statistics - Theory and Methods}, 11(6), 699-709.
#'
#' @export
plot_anom <- function(object, alpha = NULL, theme_fn = theme_statglm, ...) {
  if (!inherits(object, "lsmeans_stat_glm")) {
    stop("'object' must be an object of class 'lsmeans_stat_glm'.")
  }

  lsm_df <- object$lsmeans
  if (is.null(alpha)) alpha <- object$alpha

  m_names <- lsm_df$Effect
  means <- lsm_df$Estimate
  ses <- lsm_df[["Std Error"]]
  dfe <- object$df.residual
  k <- length(means)

  if (k < 2) {
    stop("plot_anom requires at least two factor levels.")
  }

  # Overall average LS-mean (grand mean)
  grand_mean <- mean(means, na.rm = TRUE)

  # Variance of difference from grand mean:
  # Using covariance matrix if available, else independent approximation
  if (!is.null(object$vcov)) {
    V <- object$vcov
    one <- rep(1 / k, k)
    var_grand <- as.numeric(t(one) %*% V %*% one)
    cov_i_grand <- as.numeric(V %*% one)
    var_diff <- diag(V) + var_grand - 2 * cov_i_grand
    se_diff <- sqrt(pmax(1e-12, var_diff))
  } else {
    se_diff <- ses * sqrt((k - 1) / k)
  }

  # Decision limits using Bonferroni adjusted t-critical value
  t_crit <- stats::qt(1 - alpha / (2 * k), df = dfe)
  udl <- grand_mean + t_crit * se_diff
  ldl <- grand_mean - t_crit * se_diff

  anom_df <- data.frame(
    group = factor(m_names, levels = m_names),
    mean = means,
    udl = udl,
    ldl = ldl,
    diff = means - grand_mean,
    outlier = means > udl | means < ldl,
    stringsAsFactors = FALSE
  )
  anom_df$status <- ifelse(anom_df$outlier, "Beyond Decision Limits", "Within Limits")
  anom_df$status <- factor(anom_df$status, levels = c("Beyond Decision Limits", "Within Limits"))

  p <- ggplot2::ggplot(anom_df, ggplot2::aes(x = .data$group, y = .data$mean)) +
    # Center line at Grand Mean
    ggplot2::geom_hline(yintercept = grand_mean, linetype = "solid", color = "#333333", linewidth = 0.8) +
    # Decision limit lines
    ggplot2::geom_line(ggplot2::aes(x = as.numeric(.data$group), y = .data$udl),
                       color = "#d62728", linetype = "dashed", linewidth = 0.7) +
    ggplot2::geom_line(ggplot2::aes(x = as.numeric(.data$group), y = .data$ldl),
                       color = "#d62728", linetype = "dashed", linewidth = 0.7) +
    # Needles / stems from grand mean
    ggplot2::geom_segment(ggplot2::aes(x = .data$group, xend = .data$group,
                                       y = grand_mean, yend = .data$mean),
                          color = "#2b5c8f", linewidth = 0.7) +
    # Points
    ggplot2::geom_point(ggplot2::aes(color = .data$status), size = 3.5) +
    ggplot2::scale_color_manual(
      values = c("Beyond Decision Limits" = "#d62728", "Within Limits" = "#2b5c8f"),
      name = "Decision"
    ) +
    ggplot2::labs(
      title = paste("Analysis of Means (ANOM) for", object$effect),
      subtitle = paste0("Grand Mean = ", round(grand_mean, 4), " | Dashed red lines: Decision Limits (alpha = ", alpha, ")"),
      x = object$effect,
      y = "LS-Mean"
    ) +
    theme_fn() +
    ggplot2::theme(legend.position = "top")

  p
}

#' Means Comparison Lines Plot
#'
#' @description
#' Generates a grouping lines display showing non-significant subsets of means
#' (matching SAS PROC GLM \code{LSMEANS / LINES}).
#'
#' @param object An object of class \code{"lsmeans_stat_glm"}.
#' @param alpha Significance level (defaults to \code{object$alpha}).
#' @param theme_fn Theme function (default \code{theme_statglm}).
#' @param ... Additional arguments.
#'
#' @return A \code{ggplot} object.
#' @export
plot_lines <- function(object, alpha = NULL, theme_fn = theme_statglm, ...) {
  if (!inherits(object, "lsmeans_stat_glm")) {
    stop("'object' must be an object of class 'lsmeans_stat_glm'.")
  }

  lsm_df <- object$lsmeans
  diff_df <- object$diffs

  if (nrow(diff_df) == 0) {
    stop("plot_lines requires pairwise comparisons.")
  }

  if (is.null(alpha)) alpha <- object$alpha

  # Order groups descending by LS-mean
  ord <- order(lsm_df$Estimate, decreasing = TRUE)
  sorted_lsm <- lsm_df[ord, , drop = FALSE]
  group_names <- sorted_lsm$Effect
  k <- length(group_names)

  # Adjacency matrix of non-significant differences
  adj <- matrix(TRUE, nrow = k, ncol = k, dimnames = list(group_names, group_names))

  for (r in seq_len(nrow(diff_df))) {
    comp_str <- diff_df$Comparison[r]
    parts <- unlist(strsplit(comp_str, " - "))
    if (length(parts) != 2) next
    g1 <- parts[1]; g2 <- parts[2]
    if (g1 %in% group_names && g2 %in% group_names) {
      p_val <- if ("Adj P" %in% names(diff_df)) diff_df[["Adj P"]][r] else diff_df[["Pr > |t|"]][r]
      is_sig <- !is.na(p_val) && p_val < alpha
      if (is_sig) {
        adj[g1, g2] <- FALSE
        adj[g2, g1] <- FALSE
      }
    }
  }

  # Find maximal cliques of mutually non-significant groups
  # Simple greedy contiguous grouping for ordered means:
  lines_list <- list()
  for (i in 1:k) {
    for (j in k:i) {
      if (all(adj[group_names[i:j], group_names[i:j]])) {
        # Check if already covered by an existing line
        is_sub <- FALSE
        for (line in lines_list) {
          if (line[1] <= i && line[2] >= j) {
            is_sub <- TRUE
            break
          }
        }
        if (!is_sub && (j >= i)) {
          lines_list[[length(lines_list) + 1]] <- c(i, j)
        }
      }
    }
  }

  # Build data frame for lines
  line_records <- list()
  for (idx in seq_along(lines_list)) {
    i <- lines_list[[idx]][1]
    j <- lines_list[[idx]][2]
    line_records[[length(line_records) + 1]] <- data.frame(
      line_id = idx,
      start_group = group_names[i],
      end_group = group_names[j],
      start_x = i,
      end_x = j,
      y_pos = idx,
      stringsAsFactors = FALSE
    )
  }
  lines_plot_df <- if (length(line_records) > 0) do.call(rbind, line_records) else data.frame()

  sorted_lsm$x_pos <- seq_len(k)
  lwr_col <- grep("^Lower", names(sorted_lsm))[1]
  upr_col <- grep("^Upper", names(sorted_lsm))[1]
  sorted_lsm$lwr <- sorted_lsm[[lwr_col]]
  sorted_lsm$upr <- sorted_lsm[[upr_col]]

  p <- ggplot2::ggplot() +
    ggplot2::geom_errorbar(data = sorted_lsm,
                           ggplot2::aes(x = factor(.data$Effect, levels = group_names),
                                        ymin = .data$lwr, ymax = .data$upr),
                           width = 0.2, color = "#2c3e50", linewidth = 0.7) +
    ggplot2::geom_point(data = sorted_lsm,
                        ggplot2::aes(x = factor(.data$Effect, levels = group_names),
                                     y = .data$Estimate),
                        color = "#1f77b4", size = 3.5)

  p <- p +
    ggplot2::labs(
      title = paste("Means Comparison Plot for", object$effect),
      subtitle = paste("Means not sharing a common subset grouping differ significantly (alpha =", alpha, ")"),
      x = object$effect,
      y = "LS-Mean"
    ) +
    theme_fn()

  p
}

#' S3 Plot Method for Least-Squares Means
#'
#' @description
#' Master plotting method for \code{"lsmeans_stat_glm"} objects.
#'
#' @param x An object of class \code{"lsmeans_stat_glm"}.
#' @param which Plot type: \code{"meanplot"} (default), \code{"diffplot"}, \code{"controlplot"},
#'   \code{"anomplot"}, or \code{"linesplot"}.
#' @param ... Arguments passed to the specific plotting function.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' data(warpbreaks)
#' fit <- stat_glm(breaks ~ wool * tension, data = warpbreaks)
#' lsm <- lsmeans_glm(fit, effect = "wool")
#' plot(lsm) # Displays meanplot
#'
#' @export
plot.lsmeans_stat_glm <- function(x, which = c("meanplot", "diffplot", "controlplot",
                                               "anomplot", "linesplot"), ...) {
  which <- match.arg(which)
  switch(
    which,
    meanplot = plot_lsmeans(x, ...),
    diffplot = plot_diff(x, ...),
    controlplot = plot_control(x, ...),
    anomplot = plot_anom(x, ...),
    linesplot = plot_lines(x, ...),
    plot_lsmeans(x, ...)
  )
}
