test_that("random_glm replicates SAS PROC GLM Example 53.8 benchmark", {
  # Milliken & Johnson (1984) Unbalanced Mixed Model (Machine-Person Experiment)
  # 3 machines (fixed), 6 persons (random), unbalanced replications (44 observations)
  m <- c(
    1, 1, 52.0, 1, 2, 51.8, 1, 2, 52.8, 1, 3, 60.0, 1, 4, 51.1, 1, 4, 52.3, 1, 5, 50.9,
    1, 5, 51.8, 1, 5, 51.4, 1, 6, 46.4, 1, 6, 44.8, 1, 6, 49.2, 2, 1, 64.0, 2, 2, 59.7,
    2, 2, 60.0, 2, 2, 59.0, 2, 3, 68.6, 2, 3, 65.8, 2, 4, 63.2, 2, 4, 62.8, 2, 4, 62.2,
    2, 5, 64.8, 2, 5, 65.0, 2, 6, 43.7, 2, 6, 44.2, 2, 6, 43.0, 3, 1, 67.5, 3, 1, 67.2,
    3, 1, 66.9, 3, 2, 61.5, 3, 2, 61.7, 3, 2, 62.3, 3, 3, 70.8, 3, 3, 70.6, 3, 3, 71.0,
    3, 4, 64.1, 3, 4, 66.2, 3, 4, 64.0, 3, 5, 72.1, 3, 5, 72.0, 3, 5, 71.1, 3, 6, 62.0,
    3, 6, 61.4, 3, 6, 60.5
  )
  mat <- matrix(m, ncol = 3, byrow = TRUE)
  dat <- data.frame(
    machine = factor(mat[, 1]),
    person = factor(mat[, 2]),
    rating = mat[, 3]
  )

  fit <- stat_glm(rating ~ machine * person, data = dat)
  rnd <- random_glm(fit, effects = ~ person + machine:person, test = TRUE)

  # Check Expected Mean Squares coefficients (Goodnight & Speed 1978)
  # From SAS PROC GLM Output 53.8.3:
  # machine: Var(Error) + 2.137 Var(machine*person) + Q(machine)
  # person: Var(Error) + 2.2408 Var(machine*person) + 6.7224 Var(person)
  # machine*person: Var(Error) + 2.3162 Var(machine*person)
  ems <- rnd$ems_matrix
  expect_equal(ems["machine", "Error"], 1.0, tolerance = 1e-6)
  expect_equal(ems["machine", "machine:person"], 2.1370, tolerance = 1e-4)
  expect_equal(ems["machine", "person"], 0.0, tolerance = 1e-6)

  expect_equal(ems["person", "Error"], 1.0, tolerance = 1e-6)
  expect_equal(ems["person", "person"], 6.7224, tolerance = 1e-4)
  expect_equal(ems["person", "machine:person"], 2.2408, tolerance = 1e-4)

  expect_equal(ems["machine:person", "Error"], 1.0, tolerance = 1e-6)
  expect_equal(ems["machine:person", "person"], 0.0, tolerance = 1e-6)
  expect_equal(ems["machine:person", "machine:person"], 2.3162, tolerance = 1e-4)

  # Check Mixed Model ANOVA Tests (Output 53.8.4)
  tests <- rnd$tests
  rownames(tests) <- tests$Source

  # 1. machine
  expect_equal(unname(tests["machine", "DF"]), 2)
  expect_equal(unname(tests["machine", "Type III SS"]), 1238.197626, tolerance = 1e-5)
  expect_equal(unname(tests["machine", "Mean Square"]), 619.098813, tolerance = 1e-5)
  expect_equal(unname(tests["machine", "F Value"]), 16.57, tolerance = 1e-2)
  expect_equal(unname(tests["machine", "Pr > F"]), 0.0007, tolerance = 1e-3)
  expect_equal(unname(tests["machine", "Denom DF"]), 10.036, tolerance = 1e-3)
  expect_equal(unname(tests["machine", "Denom MS"]), 37.370384, tolerance = 1e-3)

  # 2. person
  expect_equal(unname(tests["person", "DF"]), 5)
  expect_equal(unname(tests["person", "Type III SS"]), 1011.053834, tolerance = 1e-5)
  expect_equal(unname(tests["person", "Mean Square"]), 202.210767, tolerance = 1e-5)
  expect_equal(unname(tests["person", "F Value"]), 5.17, tolerance = 1e-2)
  expect_equal(unname(tests["person", "Pr > F"]), 0.0133, tolerance = 1e-2)
  expect_equal(unname(tests["person", "Denom DF"]), 10.015, tolerance = 1e-3)
  expect_equal(unname(tests["person", "Denom MS"]), 39.143708, tolerance = 1e-3)

  # 3. machine:person
  expect_equal(unname(tests["machine:person", "DF"]), 10)
  expect_equal(unname(tests["machine:person", "Type III SS"]), 404.315028, tolerance = 1e-5)
  expect_equal(unname(tests["machine:person", "Mean Square"]), 40.431503, tolerance = 1e-5)
  expect_equal(unname(tests["machine:person", "F Value"]), 46.34, tolerance = 1e-2)
  expect_lt(unname(tests["machine:person", "Pr > F"]), 0.0001)
  expect_equal(unname(tests["machine:person", "Denom DF"]), 26)
  expect_equal(unname(tests["machine:person", "Denom MS"]), 0.872564, tolerance = 1e-5)
  expect_equal(unname(tests["machine:person", "Error Term"]), "MS(Error)")
})

test_that("Flexible effect specifications yield identical results", {
  data(warpbreaks)
  fit <- stat_glm(breaks ~ wool * tension, data = warpbreaks)

  # Four equivalent ways to specify random effects
  r1 <- random_glm(fit, effects = ~ tension + wool:tension)
  r2 <- random_glm(fit, effects = c("tension", "wool:tension"))
  r3 <- random_glm(fit, effects = c("tension", "wool*tension"))
  r4 <- random_glm(fit, effects = "tension, wool*tension")

  expect_equal(r1$ems_matrix, r2$ems_matrix)
  expect_equal(r1$ems_matrix, r3$ems_matrix)
  expect_equal(r1$ems_matrix, r4$ems_matrix)
  expect_equal(r1$tests, r2$tests)
  expect_equal(r1$tests, r3$tests)
  expect_equal(r1$tests, r4$tests)
})

test_that("Integration with anova_glm and stat_glm works seamlessly", {
  data(warpbreaks)

  # Via anova_glm(..., random = ...)
  fit <- stat_glm(breaks ~ wool * tension, data = warpbreaks)
  an <- anova_glm(fit, type = 3, random = ~ tension + wool:tension)
  expect_s3_class(an$random, "random_stat_glm")
  expect_true(!is.null(an$random$tests))

  # Via stat_glm(..., random = ...)
  fit_rnd <- stat_glm(breaks ~ wool * tension, data = warpbreaks,
                      random = ~ tension + wool:tension)
  expect_equal(fit_rnd$random_effects, c("tension", "wool:tension"))

  # anova_glm automatically picks up fit$random_effects
  an_auto <- anova_glm(fit_rnd, type = 3)
  expect_equal(an_auto$random$tests, an$random$tests)

  # direct random_glm without effects picks up fit$random_effects
  rnd_auto <- random_glm(fit_rnd)
  expect_equal(rnd_auto$tests, an$random$tests)
})

test_that("Balanced two-way mixed model uses exact single error term", {
  set.seed(123)
  df_bal <- expand.grid(
    rep = 1:4,
    B = factor(1:3),
    A = factor(1:2)
  )
  df_bal$y <- 5 + 2 * (df_bal$A == "2") + rnorm(nrow(df_bal))

  fit_bal <- stat_glm(y ~ A * B, data = df_bal)
  # Declare B and A:B as random
  rnd_bal <- random_glm(fit_bal, effects = ~ B + A:B, test = TRUE)

  tests <- rnd_bal$tests
  rownames(tests) <- tests$Source

  # In balanced two-way mixed model:
  # Factor A is tested against MS(A:B) with exact DF = DF(A:B) = 2
  expect_equal(tests["A", "Error Term"], "MS(A:B)")
  expect_equal(tests["A", "Denom DF"], 2)
  expect_equal(tests["A", "Denom MS"], tests["A:B", "Mean Square"], tolerance = 1e-6)
  expect_equal(tests["A", "F Value"], tests["A", "Mean Square"] / tests["A:B", "Mean Square"], tolerance = 1e-6)

  # A:B is tested against MS(Error) with exact DF = 18
  expect_equal(tests["A:B", "Error Term"], "MS(Error)")
  expect_equal(tests["A:B", "Denom DF"], 18)
})

test_that("Q option returns valid numeric quadratic forms", {
  data(warpbreaks)
  fit <- stat_glm(breaks ~ wool * tension, data = warpbreaks)
  rnd_q <- random_glm(fit, effects = ~ tension + wool:tension, q = TRUE)

  expect_true(!is.null(rnd_q$q_matrices))
  expect_true("wool" %in% names(rnd_q$q_matrices))
  expect_true(is.matrix(rnd_q$q_matrices[["wool"]]))
  # Q matrix should be symmetric
  q_wool <- rnd_q$q_matrices[["wool"]]
  expect_equal(q_wool, t(q_wool), tolerance = 1e-7)
})

test_that("Output and S3 methods work as expected", {
  data(warpbreaks)
  fit <- stat_glm(breaks ~ wool * tension, data = warpbreaks)
  rnd <- random_glm(fit, effects = ~ tension + wool:tension)

  # as.data.frame
  df_tests <- as.data.frame(rnd)
  expect_s3_class(df_tests, "data.frame")
  expect_equal(nrow(df_tests), 3)

  # print without test
  rnd_notest <- random_glm(fit, effects = ~ tension + wool:tension, test = FALSE)
  df_ems <- as.data.frame(rnd_notest)
  expect_equal(names(df_ems), c("Source", "Expected Mean Square"))

  # print methods run cleanly
  output <- capture.output(print(rnd))
  expect_true(any(grepl("Expected Values of Type III Mean Squares", output)))
  expect_true(any(grepl("Tests of Hypotheses for Mixed Model Analysis of Variance", output)))
})

test_that("Input validation and error handling", {
  data(warpbreaks)
  fit <- stat_glm(breaks ~ wool * tension, data = warpbreaks)

  # Non-existent effect
  expect_error(random_glm(fit, effects = "invalid_effect"), "not found among model terms")

  # Missing effects when not in object
  expect_error(random_glm(fit), "No random effects specified")

  # Non stat_glm object
  expect_error(random_glm(lm(breaks ~ wool, data = warpbreaks)), "must be an object of class 'stat_glm'")

  # Invalid SS type
  expect_error(random_glm(fit, effects = "wool", type = 2), "'type' must be either 1 or 3")
})

test_that("random_glm handles nested random effects with rank-deficient designs", {
  # Nested design: B nested in A, B is random
  dat_nested <- data.frame(
    A = factor(c(rep("1", 8), rep("2", 12))),
    B = factor(c(rep("1", 4), rep("2", 4), rep("3", 4), rep("4", 4), rep("5", 4))),
    y = c(10, 12, 11, 13,  15, 14, 16, 17,  20, 22, 21, 23,  25, 24, 26, 27,  30, 31, 32, 29)
  )
  fit_nested <- stat_glm(y ~ A + A:B, data = dat_nested)
  rnd <- random_glm(fit_nested, effects = ~ A:B, test = TRUE)

  # Check EMS matrix
  expect_equal(rnd$ems_matrix["A:B", "Error"], 1.0)
  expect_equal(unname(rnd$ems_matrix["A:B", "A:B"]), 4.0, tolerance = 1e-5)

  # A:B should test against MS(Error)
  tests <- rnd$tests
  rownames(tests) <- tests$Source
  expect_equal(tests["A:B", "Error Term"], "MS(Error)")
  expect_equal(unname(tests["A:B", "DF"]), 3)
})
