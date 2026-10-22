test_that("Type I sums of squares match base R anova", {
  data(warpbreaks)
  fit_glm <- stat_glm(breaks ~ wool * tension, data = warpbreaks)
  fit_lm <- lm(breaks ~ wool * tension, data = warpbreaks)

  an_glm <- anova_glm(fit_glm, type = 1)
  an_base <- anova(fit_lm)

  tab1 <- an_glm$tables[["Type I"]]
  terms <- rownames(tab1)

  for (tm in terms) {
    expect_equal(tab1[tm, "DF"], an_base[tm, "Df"])
    expect_equal(tab1[tm, "Sum of Squares"], an_base[tm, "Sum Sq"], tolerance = 1e-6)
    expect_equal(tab1[tm, "Mean Square"], an_base[tm, "Mean Sq"], tolerance = 1e-6)
    expect_equal(tab1[tm, "F Value"], an_base[tm, "F value"], tolerance = 1e-6)
    expect_equal(tab1[tm, "Pr > F"], an_base[tm, "Pr(>F)"], tolerance = 1e-6)
  }
})

test_that("Type II sums of squares match car::Anova(type = 2)", {
  data(warpbreaks)
  fit_glm <- stat_glm(breaks ~ wool * tension, data = warpbreaks)
  fit_lm <- lm(breaks ~ wool * tension, data = warpbreaks)

  an_glm <- anova_glm(fit_glm, type = 2)
  an_car <- car::Anova(fit_lm, type = 2)

  tab2 <- an_glm$tables[["Type II"]]
  terms <- rownames(tab2)

  for (tm in terms) {
    expect_equal(tab2[tm, "DF"], an_car[tm, "Df"])
    expect_equal(tab2[tm, "Sum of Squares"], an_car[tm, "Sum Sq"], tolerance = 1e-6)
    expect_equal(tab2[tm, "F Value"], an_car[tm, "F value"], tolerance = 1e-6)
    expect_equal(tab2[tm, "Pr > F"], an_car[tm, "Pr(>F)"], tolerance = 1e-6)
  }
})

test_that("Type III sums of squares match car::Anova(type = 3)", {
  data(warpbreaks)
  fit_glm <- stat_glm(breaks ~ wool * tension, data = warpbreaks)
  fit_lm <- lm(breaks ~ wool * tension, data = warpbreaks,
               contrasts = list(wool = contr.sum, tension = contr.sum))

  an_glm <- anova_glm(fit_glm, type = 3)
  an_car <- car::Anova(fit_lm, type = 3)

  tab3 <- an_glm$tables[["Type III"]]
  terms <- rownames(tab3)

  for (tm in terms) {
    expect_equal(tab3[tm, "DF"], an_car[tm, "Df"])
    expect_equal(tab3[tm, "Sum of Squares"], an_car[tm, "Sum Sq"], tolerance = 1e-6)
    expect_equal(tab3[tm, "F Value"], an_car[tm, "F value"], tolerance = 1e-6)
    expect_equal(tab3[tm, "Pr > F"], an_car[tm, "Pr(>F)"], tolerance = 1e-6)
  }
})

test_that("SAS GLM documentation benchmark dataset produces exact numerical outputs", {
  # Benchmark dataset from Chapter 53 (Hypothesis Testing in PROC GLM, p. 75-80)
  a <- factor(c(1, 1, 1, 2, 2, 2, 3, 3, 3, 3))
  b <- factor(c(1, 1, 2, 1, 2, 2, 1, 1, 2, 2))
  y <- c(23.5, 23.7, 28.7, 8.9, 5.6, 8.9, 10.3, 12.5, 13.6, 14.6)
  df <- data.frame(a, b, y)

  fit <- stat_glm(y ~ a * b, data = df)

  # Overall Model ANOVA checks
  expect_equal(fit$ssm, 520.4760000, tolerance = 1e-6)
  expect_equal(fit$sse, 8.3850000, tolerance = 1e-6)
  expect_equal(fit$sst, 528.8610000, tolerance = 1e-6)
  expect_equal(fit$df.model, 5L)
  expect_equal(fit$df.residual, 4L)
  expect_equal(fit$r_square, 0.984145, tolerance = 1e-5)
  expect_equal(fit$root_mse, 1.447843, tolerance = 1e-5)
  expect_equal(fit$cv, 9.633022, tolerance = 1e-4)

  # ANOVA Tables Types 1, 2, and 3
  an_all <- anova_glm(fit, type = c(1, 2, 3))

  # Type I Sums of Squares (p. 78)
  t1 <- an_all$tables[["Type I"]]
  expect_equal(t1["a", "DF"], 2L)
  expect_equal(t1["a", "Sum of Squares"], 494.0310000, tolerance = 1e-5)
  expect_equal(t1["b", "DF"], 1L)
  expect_equal(t1["b", "Sum of Squares"], 10.7142857, tolerance = 1e-5)
  expect_equal(t1["a:b", "DF"], 2L)
  expect_equal(t1["a:b", "Sum of Squares"], 15.7307143, tolerance = 1e-5)

  # Type II Sums of Squares (p. 79)
  t2 <- an_all$tables[["Type II"]]
  expect_equal(t2["a", "DF"], 2L)
  expect_equal(t2["a", "Sum of Squares"], 499.1202857, tolerance = 1e-5)
  expect_equal(t2["b", "DF"], 1L)
  expect_equal(t2["b", "Sum of Squares"], 10.7142857, tolerance = 1e-5)
  expect_equal(t2["a:b", "DF"], 2L)
  expect_equal(t2["a:b", "Sum of Squares"], 15.7307143, tolerance = 1e-5)

  # Type III Sums of Squares (p. 80)
  t3 <- an_all$tables[["Type III"]]
  expect_equal(t3["a", "DF"], 2L)
  expect_equal(t3["a", "Sum of Squares"], 479.1078571, tolerance = 1e-5)
  expect_equal(t3["b", "DF"], 1L)
  expect_equal(t3["b", "Sum of Squares"], 9.4556250, tolerance = 1e-5)
  expect_equal(t3["a:b", "DF"], 2L)
  expect_equal(t3["a:b", "Sum of Squares"], 15.7307143, tolerance = 1e-5)
})

test_that("Effect sizes and confidence limits match documentation formulas", {
  a <- factor(c(1, 1, 1, 2, 2, 2, 3, 3, 3, 3))
  b <- factor(c(1, 1, 2, 1, 2, 2, 1, 1, 2, 2))
  y <- c(23.5, 23.7, 28.7, 8.9, 5.6, 8.9, 10.3, 12.5, 13.6, 14.6)
  df <- data.frame(a, b, y)

  fit <- stat_glm(y ~ a * b, data = df)
  an_eff <- anova_glm(fit, type = 3, effectsize = TRUE, alpha = 0.10)
  t3 <- an_eff$tables[["Type III"]]

  # Partial Eta2: SS / (SS + SSE)
  # For term a: 479.1078571 / (479.1078571 + 8.3850000) = 0.9828
  expect_equal(t3["a", "Partial Eta2"], 479.1078571 / (479.1078571 + 8.385), tolerance = 1e-5)

  # Semipartial Eta2: SS / SST
  # For term a: 479.1078571 / 528.8610000 = 0.9059
  expect_equal(t3["a", "Semipartial Eta2"], 479.1078571 / 528.861, tolerance = 1e-5)
})

test_that("Custom error term (TEST statement) correctly updates F and p-values", {
  data(warpbreaks)
  fit <- stat_glm(breaks ~ wool * tension, data = warpbreaks)

  # Test main effect 'wool' with interaction 'wool:tension' as denominator
  an_custom <- anova_glm(fit, type = 3, error_term = "wool:tension")
  t3_cust <- an_custom$tables[["Type III"]]

  # Check denominator MS used is MS of wool:tension
  ms_interaction <- t3_cust["wool:tension", "Mean Square"]
  expected_f <- t3_cust["wool", "Mean Square"] / ms_interaction
  expect_equal(t3_cust["wool", "F Value"], expected_f, tolerance = 1e-6)
})

test_that("Type III sums of squares handles nested unbalanced designs correctly", {
  # Nested design: factor B nested in factor A, unbalanced levels per branch
  # A has 2 levels: A1 has B in {1, 2}, A2 has B in {3, 4, 5} (total 5 levels of B, but 3 DF nested)
  dat_nested <- data.frame(
    A = factor(c(rep("1", 8), rep("2", 12))),
    B = factor(c(rep("1", 4), rep("2", 4), rep("3", 4), rep("4", 4), rep("5", 4))),
    y = c(10, 12, 11, 13,  15, 14, 16, 17,  20, 22, 21, 23,  25, 24, 26, 27,  30, 31, 32, 29)
  )
  fit_nested <- stat_glm(y ~ A + A:B, data = dat_nested)
  an_res <- anova_glm(fit_nested, type = c(1, 3))

  t1 <- an_res$tables[["Type I"]]
  t3 <- an_res$tables[["Type III"]]

  # Type I and Type III DF for nested factor A:B should be 5 - 2 = 3 DF
  expect_equal(t1["A:B", "DF"], 3L)
  expect_equal(t3["A:B", "DF"], 3L)
  expect_equal(t1["A", "DF"], 1L)
  expect_equal(t3["A", "DF"], 1L)

  # For non-contained term A:B, Type I SS and Type III SS should agree identically
  expect_equal(t3["A:B", "Sum of Squares"], t1["A:B", "Sum of Squares"], tolerance = 1e-6)
})
