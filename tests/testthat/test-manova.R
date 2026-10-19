test_that("manova_glm matches base R manova statistics on iris dataset", {
  data(iris)
  fit_glm <- stat_glm(cbind(Sepal.Length, Sepal.Width) ~ Species, data = iris)
  fit_base <- manova(cbind(Sepal.Length, Sepal.Width) ~ Species, data = iris)

  man_res <- manova_glm(fit_glm)
  tab <- man_res$results[["Species"]]
  rownames(tab) <- tab$Statistic

  # 1. Wilks
  base_w <- summary(fit_base, test = "Wilks")$stats["Species", ]
  expect_equal(tab["Wilks' Lambda", "Value"], as.numeric(base_w["Wilks"]), tolerance = 1e-4)
  expect_equal(tab["Wilks' Lambda", "F Value"], as.numeric(base_w["approx F"]), tolerance = 1e-4)
  expect_equal(tab["Wilks' Lambda", "Num DF"], as.numeric(base_w["num Df"]))
  expect_equal(tab["Wilks' Lambda", "Den DF"], as.numeric(base_w["den Df"]))
  expect_equal(tab["Wilks' Lambda", "Pr > F"], as.numeric(base_w["Pr(>F)"]), tolerance = 1e-4)

  # 2. Pillai
  base_p <- summary(fit_base, test = "Pillai")$stats["Species", ]
  expect_equal(tab["Pillai's Trace", "Value"], as.numeric(base_p["Pillai"]), tolerance = 1e-4)
  expect_equal(tab["Pillai's Trace", "F Value"], as.numeric(base_p["approx F"]), tolerance = 1e-4)
  expect_equal(tab["Pillai's Trace", "Num DF"], as.numeric(base_p["num Df"]))
  expect_equal(tab["Pillai's Trace", "Den DF"], as.numeric(base_p["den Df"]))
  expect_equal(tab["Pillai's Trace", "Pr > F"], as.numeric(base_p["Pr(>F)"]), tolerance = 1e-4)

  # 3. Hotelling-Lawley
  base_h <- summary(fit_base, test = "Hotelling-Lawley")$stats["Species", ]
  expect_equal(tab["Hotelling-Lawley Trace", "Value"], as.numeric(base_h["Hotelling-Lawley"]), tolerance = 1e-4)
  expect_equal(tab["Hotelling-Lawley Trace", "F Value"], as.numeric(base_h["approx F"]), tolerance = 1e-4)
  expect_equal(tab["Hotelling-Lawley Trace", "Num DF"], as.numeric(base_h["num Df"]))
  expect_equal(tab["Hotelling-Lawley Trace", "Den DF"], as.numeric(base_h["den Df"]))
  expect_equal(tab["Hotelling-Lawley Trace", "Pr > F"], as.numeric(base_h["Pr(>F)"]), tolerance = 1e-4)

  # 4. Roy's Greatest Root
  base_r <- summary(fit_base, test = "Roy")$stats["Species", ]
  expect_equal(tab["Roy's Greatest Root", "Value"], as.numeric(base_r["Roy"]), tolerance = 1e-4)
  expect_equal(tab["Roy's Greatest Root", "F Value"], as.numeric(base_r["approx F"]), tolerance = 1e-4)
  expect_equal(tab["Roy's Greatest Root", "Num DF"], as.numeric(base_r["num Df"]))
  expect_equal(tab["Roy's Greatest Root", "Den DF"], as.numeric(base_r["den Df"]))
  expect_equal(tab["Roy's Greatest Root", "Pr > F"], as.numeric(base_r["Pr(>F)"]), tolerance = 1e-4)
})
