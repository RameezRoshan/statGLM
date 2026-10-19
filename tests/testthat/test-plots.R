test_that("theme_statglm produces valid ggplot2 theme", {
  th1 <- theme_statglm()
  expect_s3_class(th1, "theme")

  th_none <- theme_statglm(grid = "none")
  expect_s3_class(th_none, "theme")

  th_both <- theme_statglm(grid = "both")
  expect_s3_class(th_both, "theme")
})

test_that("plot_diagnostics produces paneled and unpacked plots matching SAS PLOTS=DIAGNOSTICS", {
  data(warpbreaks)
  fit <- stat_glm(breaks ~ wool * tension, data = warpbreaks)

  # Default paneled grob
  p_panel <- plot_diagnostics(fit)
  expect_true(inherits(p_panel, "gtable") || inherits(p_panel, "grob"))

  # Unpacked individual ggplots
  p_unpacked <- plot_diagnostics(fit, unpack = TRUE, label = TRUE)
  expect_type(p_unpacked, "list")
  expect_equal(length(p_unpacked), 8)
  expect_named(p_unpacked, c(
    "ResidualByPredicted", "RStudentByPredicted", "ObservedByPredicted",
    "RStudentByLeverage", "CooksDPlot", "QQPlot", "ResidualHistogram", "RFPlot"
  ))

  for (nm in names(p_unpacked)) {
    expect_s3_class(p_unpacked[[nm]], "ggplot")
  }
})

test_that("plot_residuals produces residual vs covariate plots matching SAS PLOTS=RESIDUALS", {
  fit <- stat_glm(mpg ~ wt + hp, data = mtcars)

  p_res <- plot_residuals(fit, unpack = TRUE)
  expect_type(p_res, "list")
  expect_equal(length(p_res), 2)
  expect_s3_class(p_res$wt, "ggplot")
  expect_s3_class(p_res$hp, "ggplot")

  # Combined grob
  p_comb <- plot_residuals(fit, unpack = FALSE)
  expect_true(inherits(p_comb, "gtable") || inherits(p_comb, "grob"))
})

test_that("plot_fit produces continuous predictor fit plot with CLM and CLI bands", {
  fit <- stat_glm(dist ~ speed, data = cars)

  # Default with both CLM and CLI
  p_fit <- plot_fit(fit)
  expect_s3_class(p_fit, "ggplot")

  # Only CLM
  p_clm <- plot_fit(fit, clm = TRUE, cli = FALSE)
  expect_s3_class(p_clm, "ggplot")

  # Only CLI
  p_cli <- plot_fit(fit, clm = FALSE, cli = TRUE)
  expect_s3_class(p_cli, "ggplot")
})

test_that("plot_contour produces 2D surface contour plot with observed points", {
  fit <- stat_glm(mpg ~ wt + hp, data = mtcars)

  p_c1 <- plot_contour(fit, obs = "outlinegradient")
  expect_s3_class(p_c1, "ggplot")

  p_c2 <- plot_contour(fit, obs = "gradient")
  expect_s3_class(p_c2, "ggplot")

  p_c3 <- plot_contour(fit, obs = "none")
  expect_s3_class(p_c3, "ggplot")
})

test_that("plot_ancova produces ANCOVA regression lines by group", {
  df <- data.frame(
    y = c(10, 12, 14, 18, 20, 22),
    group = factor(c("A", "A", "A", "B", "B", "B")),
    x = c(1, 2, 3, 1, 2, 3)
  )
  fit <- stat_glm(y ~ group + x, data = df)

  p_anc <- plot_ancova(fit, clm = TRUE, cli = TRUE)
  expect_s3_class(p_anc, "ggplot")
})

test_that("plot_boxplot produces grouped box plots with mean diamonds and jitter", {
  data(warpbreaks)
  fit <- stat_glm(breaks ~ wool, data = warpbreaks)

  p_box <- plot_boxplot(fit, jitter = TRUE, mean_diamond = TRUE)
  expect_s3_class(p_box, "ggplot")
})

test_that("plot_interaction produces two-way interaction plots", {
  data(warpbreaks)
  fit <- stat_glm(breaks ~ wool * tension, data = warpbreaks)

  p_int <- plot_interaction(fit, clm = TRUE)
  expect_s3_class(p_int, "ggplot")
})

test_that("plot_nested produces nested classification box plots", {
  df <- data.frame(
    y = rnorm(16),
    A = factor(rep(c("A1", "A2"), each = 8)),
    B = factor(rep(c("B1", "B2", "B3", "B4"), each = 4))
  )
  fit <- stat_glm(y ~ A + B, data = df)

  p_nest <- plot_nested(fit)
  expect_s3_class(p_nest, "ggplot")
})

test_that("plot.stat_glm intelligent auto-detection and which dispatcher work properly", {
  data(warpbreaks)

  # One-way ANOVA -> auto boxplot
  fit1 <- stat_glm(breaks ~ wool, data = warpbreaks)
  p1 <- plot(fit1)
  expect_s3_class(p1, "ggplot")

  # Two-way ANOVA -> auto intplot
  fit2 <- stat_glm(breaks ~ wool * tension, data = warpbreaks)
  p2 <- plot(fit2)
  expect_s3_class(p2, "ggplot")

  # Single continuous -> auto fitplot
  fit3 <- stat_glm(dist ~ speed, data = cars)
  p3 <- plot(fit3)
  expect_s3_class(p3, "ggplot")

  # Two continuous -> auto contourfit
  fit4 <- stat_glm(mpg ~ wt + hp, data = mtcars)
  p4 <- plot(fit4)
  expect_s3_class(p4, "ggplot")

  # ANCOVA -> auto ancovaplot
  df_anc <- data.frame(y = 1:6, g = factor(rep(c("A", "B"), 3)), x = 1:6)
  fit5 <- stat_glm(y ~ g + x, data = df_anc)
  p5 <- plot(fit5)
  expect_s3_class(p5, "ggplot")

  # Explicit which dispatcher
  p_diag <- plot(fit2, which = "diagnostics", unpack = TRUE)
  expect_type(p_diag, "list")

  # which = "all"
  all_plots <- plot(fit2, which = "all")
  expect_type(all_plots, "list")
  expect_true("diagnostics" %in% names(all_plots))
  expect_true("intplot" %in% names(all_plots))
})

test_that("LS-means plots (meanplot, diffplot, controlplot, anomplot, linesplot) work properly", {
  data(warpbreaks)
  fit <- stat_glm(breaks ~ wool * tension, data = warpbreaks)
  lsm <- lsmeans_glm(fit, effect = "tension", adjust = "tukey")

  # 1. meanplot
  p_mean <- plot_lsmeans(lsm, cl = TRUE, connect = TRUE, order = "ascending")
  expect_s3_class(p_mean, "ggplot")

  # 2. diffplot (Hsu's diffogram)
  p_diff <- plot_diff(lsm, center = TRUE, nolines = FALSE)
  expect_s3_class(p_diff, "ggplot")

  # 3. controlplot (Dunnett comparisons)
  p_ctrl <- plot_control(lsm, control = "L")
  expect_s3_class(p_ctrl, "ggplot")

  # 4. anomplot (Analysis of Means)
  p_anom <- plot_anom(lsm)
  expect_s3_class(p_anom, "ggplot")

  # 5. linesplot (subset comparison lines)
  p_lines <- plot_lines(lsm)
  expect_s3_class(p_lines, "ggplot")

  # 6. S3 plot.lsmeans_stat_glm
  expect_s3_class(plot(lsm, which = "meanplot"), "ggplot")
  expect_s3_class(plot(lsm, which = "diffplot"), "ggplot")
  expect_s3_class(plot(lsm, which = "controlplot"), "ggplot")
  expect_s3_class(plot(lsm, which = "anomplot"), "ggplot")
  expect_s3_class(plot(lsm, which = "linesplot"), "ggplot")
})

test_that("save_glm_plot successfully saves high-definition output", {
  fit <- stat_glm(dist ~ speed, data = cars)
  p <- plot_fit(fit)
  tmp_file <- tempfile(fileext = ".png")
  save_glm_plot(p, tmp_file, width = 6, height = 4, dpi = 150)
  expect_true(file.exists(tmp_file))
  unlink(tmp_file)
})
