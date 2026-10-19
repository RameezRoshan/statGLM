test_that("stat_glm matches lm on standard linear regression", {
  data(mtcars)
  fit_glm <- stat_glm(mpg ~ wt + qsec + hp, data = mtcars)
  fit_lm <- lm(mpg ~ wt + qsec + hp, data = mtcars)

  # Fitted values and residuals
  expect_equal(as.vector(fitted(fit_glm)), as.vector(fitted(fit_lm)), tolerance = 1e-7)
  expect_equal(as.vector(residuals(fit_glm)), as.vector(residuals(fit_lm)), tolerance = 1e-7)

  # Error SS and MSE
  expect_equal(fit_glm$sse, deviance(fit_lm), tolerance = 1e-7)
  expect_equal(fit_glm$df.residual, fit_lm$df.residual)
  expect_equal(fit_glm$mse, summary(fit_lm)$sigma^2, tolerance = 1e-7)

  # R-square and Adj R-square
  expect_equal(fit_glm$r_square, summary(fit_lm)$r.squared, tolerance = 1e-7)
  expect_equal(fit_glm$adj_r_square, summary(fit_lm)$adj.r.squared, tolerance = 1e-7)

  # Overall F and p-value
  lm_f <- summary(fit_lm)$fstatistic
  expect_equal(fit_glm$f_value, as.numeric(lm_f[1]), tolerance = 1e-7)
  expect_equal(fit_glm$df.model, as.integer(lm_f[2]))
  expect_equal(fit_glm$df.residual, as.integer(lm_f[3]))
})

test_that("stat_glm matches lm with weights", {
  data(mtcars)
  w <- mtcars$cyl / mean(mtcars$cyl)
  fit_glm <- stat_glm(mpg ~ wt + disp, data = mtcars, weights = w)
  fit_lm <- lm(mpg ~ wt + disp, data = mtcars, weights = w)

  expect_equal(as.vector(fitted(fit_glm)), as.vector(fitted(fit_lm)), tolerance = 1e-7)
  expect_equal(fit_glm$sse, deviance(fit_lm), tolerance = 1e-7)
  expect_equal(fit_glm$mse, summary(fit_lm)$sigma^2, tolerance = 1e-7)
  expect_equal(fit_glm$r_square, summary(fit_lm)$r.squared, tolerance = 1e-7)
})

test_that("stat_glm S3 methods work as expected", {
  data(warpbreaks)
  fit <- stat_glm(breaks ~ wool + tension, data = warpbreaks)

  # summary
  s <- summary(fit)
  expect_s3_class(s, "summary.stat_glm")
  expect_true(is.data.frame(s$coefficients))

  # coef
  cf <- coef(fit)
  expect_true(is.numeric(cf))
  expect_equal(length(cf), ncol(fit$design$X_ind))

  # vcov
  v <- vcov(fit)
  expect_equal(dim(v), c(length(cf), length(cf)))

  # predict
  pred_orig <- predict(fit)
  expect_equal(pred_orig, fitted(fit))

  new_data <- data.frame(wool = factor(c("A", "B"), levels = c("A", "B")),
                         tension = factor(c("L", "M"), levels = c("L", "M", "H")))
  pred_new <- predict(fit, newdata = new_data)
  expect_equal(length(pred_new), 2)
})
