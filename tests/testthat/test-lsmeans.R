test_that("lsmeans_glm matches analytical marginal cell means on benchmark data", {
  a <- factor(c(1, 1, 1, 2, 2, 2, 3, 3, 3, 3))
  b <- factor(c(1, 1, 2, 1, 2, 2, 1, 1, 2, 2))
  y <- c(23.5, 23.7, 28.7, 8.9, 5.6, 8.9, 10.3, 12.5, 13.6, 14.6)
  df <- data.frame(a, b, y)

  fit <- stat_glm(y ~ a * b, data = df)
  lsm_a <- lsmeans_glm(fit, effect = "a")

  tab <- lsm_a$lsmeans
  rownames(tab) <- tab$Effect

  # Cell marginal means
  expect_equal(tab["a 1", "Estimate"], 26.15, tolerance = 1e-6)
  expect_equal(tab["a 2", "Estimate"], 8.075, tolerance = 1e-6)
  expect_equal(tab["a 3", "Estimate"], 12.75, tolerance = 1e-6)

  # Standard errors
  expect_equal(tab["a 1", "Std Error"], 0.8866193, tolerance = 1e-5)
  expect_equal(tab["a 2", "Std Error"], 0.8866193, tolerance = 1e-5)
  expect_equal(tab["a 3", "Std Error"], 0.7239216, tolerance = 1e-5)

  # Pairwise differences
  diffs <- lsm_a$diffs
  rownames(diffs) <- diffs$Comparison

  expect_equal(diffs["a 1 - a 2", "Estimate"], 26.15 - 8.075, tolerance = 1e-6)
  expect_equal(diffs["a 1 - a 3", "Estimate"], 26.15 - 12.75, tolerance = 1e-6)
  expect_equal(diffs["a 2 - a 3", "Estimate"], 8.075 - 12.75, tolerance = 1e-6)
})

test_that("lsmeans_glm multiple comparison adjustments work properly", {
  data(warpbreaks)
  fit <- stat_glm(breaks ~ tension, data = warpbreaks)

  lsm_tukey <- lsmeans_glm(fit, effect = "tension", adjust = "tukey")
  lsm_bonf <- lsmeans_glm(fit, effect = "tension", adjust = "bonferroni")
  lsm_none <- lsmeans_glm(fit, effect = "tension", adjust = "none")

  # Raw p-values should be identical
  expect_equal(lsm_tukey$diffs$`Pr > |t|`, lsm_none$diffs$`Pr > |t|`)
  expect_equal(lsm_bonf$diffs$`Pr > |t|`, lsm_none$diffs$`Pr > |t|`)

  # Bonferroni should be >= raw p
  expect_true(all(lsm_bonf$diffs$`Adj P` >= lsm_none$diffs$`Adj P`))

  # Tukey should be between raw and Bonferroni (or <= Bonferroni)
  expect_true(all(lsm_tukey$diffs$`Adj P` >= lsm_none$diffs$`Adj P`))
  expect_true(all(lsm_tukey$diffs$`Adj P` <= lsm_bonf$diffs$`Adj P` + 1e-8))
})
