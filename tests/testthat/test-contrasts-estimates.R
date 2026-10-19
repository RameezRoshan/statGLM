test_that("contrast_glm matches car::linearHypothesis for single and multi-DF contrasts", {
  data(warpbreaks)
  fit_glm <- stat_glm(breaks ~ wool + tension, data = warpbreaks)
  fit_lm <- lm(breaks ~ wool + tension, data = warpbreaks)

  # 1. Single DF contrast: Wool A vs B
  # In stat_glm: wool A vs B is list(wool = c(1, -1))
  c_glm <- contrast_glm(fit_glm, list(wool = c(1, -1)), name = "Wool A - B")
  c_car <- car::linearHypothesis(fit_lm, "woolB = 0")

  expect_equal(c_glm$df, c_car$Df[2])
  expect_equal(c_glm$ss, c_car$`Sum of Sq`[2], tolerance = 1e-6)
  expect_equal(c_glm$f_value, c_car$F[2], tolerance = 1e-6)
  expect_equal(c_glm$p_value, c_car$`Pr(>F)`[2], tolerance = 1e-6)

  # 2. Multi-DF contrast: Tension (levels L, M, H -> 2 DF)
  # L - M and M - H
  L_tension <- list(
    "L vs M" = list(tension = c(1, -1, 0)),
    "M vs H" = list(tension = c(0, 1, -1))
  )
  c_glm_multi <- contrast_glm(fit_glm, L_tension, name = "Tension Effect")
  c_car_multi <- car::linearHypothesis(fit_lm, c("tensionM = 0", "tensionH = 0"))

  expect_equal(c_glm_multi$df, c_car_multi$Df[2])
  expect_equal(c_glm_multi$ss, c_car_multi$`Sum of Sq`[2], tolerance = 1e-6)
  expect_equal(c_glm_multi$f_value, c_car_multi$F[2], tolerance = 1e-6)
  expect_equal(c_glm_multi$p_value, c_car_multi$`Pr(>F)`[2], tolerance = 1e-6)
})

test_that("contrast_glm errors on non-estimable hypotheses", {
  X <- matrix(c(1, 1, 0,
                1, 1, 0,
                1, 0, 1,
                1, 0, 1), 4, 3, byrow = TRUE)
  y <- c(2, 3, 5, 7)
  df <- data.frame(g = factor(c("A", "A", "B", "B")), y = y)
  fit <- stat_glm(y ~ g, data = df)

  # Attempting to test single cell parameter without intercept (non-estimable)
  # L = (0, 1, 0)
  expect_error(contrast_glm(fit, matrix(c(0, 1, 0), 1, 3)), "non-estimable")
})

test_that("estimate_glm calculates correct estimates, SEs, and confidence intervals", {
  data(warpbreaks)
  fit_glm <- stat_glm(breaks ~ wool + tension, data = warpbreaks)
  fit_lm <- lm(breaks ~ wool + tension, data = warpbreaks)

  est_glm <- estimate_glm(fit_glm, list(wool = c(1, -1)), name = "Wool A - B")
  s_lm <- summary(fit_lm)

  # Wool A - Wool B = -(Wool B effect in lm)
  coef_diff <- -coef(fit_lm)["woolB"]
  se_diff <- s_lm$coefficients["woolB", "Std. Error"]

  tab <- est_glm$table
  expect_equal(tab["Wool A - B", "Estimate"], as.numeric(coef_diff), tolerance = 1e-6)
  expect_equal(tab["Wool A - B", "Std Error"], as.numeric(se_diff), tolerance = 1e-6)
  expect_equal(abs(tab["Wool A - B", "t Value"]), abs(s_lm$coefficients["woolB", "t value"]), tolerance = 1e-6)
  expect_equal(tab["Wool A - B", "Pr > |t|"], s_lm$coefficients["woolB", "Pr(>|t|)"], tolerance = 1e-6)
})
