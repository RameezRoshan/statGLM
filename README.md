# statGLM: Modern Statistical General Linear Models in R

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R-CMD-check](https://img.shields.io/badge/R--CMD--check-passing-brightgreen.svg)](#)

`statGLM` is a production-grade, clean-room R package implementing classical general linear models, estimable functions, sequential and marginal sums of squares (Types I, II, and III), least-squares means (population marginal means), custom linear hypotheses, and multivariate analysis of variance (MANOVA).

Developed strictly from foundational mathematical definitions (Searle 1971, Goodnight 1978, 1979), `statGLM` provides rigorous, reliable analytical workflows with complete numerical fidelity, zero proprietary dependencies, and a clean S3 class interface.

---

## Key Features

- **Cell-Indicator & Over-Parameterized Parameterization**: Fits fixed-effects models directly using indicator columns for all factor levels, retaining natural parameter mappings ($\mu, \alpha_1, \alpha_2, \ldots$) while resolving normal equations via reflexive generalized inverses ($g_2$) or Moore-Penrose pseudo-inverses.
- **Estimability Verification**: Evaluates whether arbitrary single or multi-degree-of-freedom hypotheses $L \beta = c$ satisfy the Gauss-Markov estimability condition $L (X'X)^- (X'X) = L$.
- **Type I, Type II, and Type III Sums of Squares**:
  - **Type I (Sequential)**: Incremental reduction in error sums of squares adding terms in user-specified order.
  - **Type II (Hierarchical)**: Reduction in error sums of squares adding each term after all non-containing terms, respecting principle of marginality.
  - **Type III (Orthogonal/Partial)**: Tests each effect adjusting for all other effects in the model with equal weighting across interacting levels.
- **Effect Sizes & Confidence Intervals**: Partial $\eta^2$, semipartial $\eta^2$, partial $\omega^2$, semipartial $\omega^2$, noncentrality parameters (minimum-variance unbiased and low-MSE), and noncentral $F$ confidence limits.
- **Least-Squares Means (LS-Means)**: Automated construction of population marginal means with covariate adjustments (`at`) and pairwise group comparisons with adjustments for multiple testing (Tukey-Kramer, Bonferroni, Scheffé, Šidák).
- **Custom Linear Hypotheses (`contrast_glm` & `estimate_glm`)**: Multi-DF $F$-tests and single-DF $t$-tests with point estimates, standard errors, and confidence intervals. Supports both matrix and intuitive named-list specifications.
- **Multivariate Analysis of Variance (MANOVA)**: Evaluates multivariate response models via Wilks' Lambda, Pillai's Trace, Hotelling-Lawley Trace, and Roy's Greatest Root with exact and Rao approximate $F$-tests.
- **Custom Error Terms**: Evaluates hypotheses using specified model effects as the denominator error term (analogous to classical `TEST H = effect E = error` statements).

---

## Mathematical Foundations

### 1. General Linear Model Specification

For $n$ observations and $p$ model columns:

$$y = X \beta + \epsilon, \quad \epsilon \sim \mathcal{N}(0, \sigma^2 W^{-1})$$

The normal equations are given by:

$$(X' W X) \beta = X' W y$$

A solution vector is computed via generalized inverse $G = (X' W X)^-$:

$$b = G X' W y$$

Under the reflexive $g_2$ inverse (via the Goodnight SWEEP operator), redundant columns corresponding to linearly dependent parameters are explicitly set to zero, matching traditional over-parameterized model displays.

### 2. Estimable Functions

A linear combination $L \beta$ is **estimable** if and only if:

$$\| L G (X' W X) - L \|_\infty \le \epsilon$$

For any estimable function, the estimate $\hat{\theta} = L b$ and its variance $\mathrm{Var}(\hat{\theta}) = \sigma^2 L G L'$ are strictly **invariant** to the choice of generalized inverse $G$.

### 3. General Linear Hypothesis Test

For hypothesis $H_0: L \beta = c$ with $\mathrm{rank}(L) = r$:

$$\mathrm{SS}(H) = (L b - c)' [L G L']^{-1} (L b - c)$$

$$F = \frac{\mathrm{SS}(H) / r}{\mathrm{MSE}} \sim \mathcal{F}(r, \mathrm{DF}_E)$$

### 4. Least-Squares Means (Marginal Means)

For effect $E$, the linear combination vector $L$ is constructed by:
1. Setting the intercept coefficient to 1.
2. Setting the level of effect $E$ to 1 (0 for other levels).
3. Assigning equal weights $1/k$ across levels of containing and interacting factors.
4. Setting continuous covariates to their sample mean $\bar{x}$ (or user-defined reference values).

---

## Installation

You can install the development version of `statGLM` directly from GitHub or build from source:

```r
# Install remotes if needed
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# Install statGLM from GitHub
remotes::install_github("RameezRoshan/statGLM")

# Or install locally
devtools::install(".")
```

---

## Quick-Start Guide

### 1. Model Fitting

```r
library(statGLM)

# Load example data
data(warpbreaks)

# Fit general linear model
fit <- stat_glm(breaks ~ wool * tension, data = warpbreaks)

# Display model fit and overall ANOVA
summary(fit)
```

### 2. Sequential, Hierarchical, and Orthogonal ANOVA

```r
# Type I (Sequential) Sums of Squares
anova_glm(fit, type = 1)

# Type II (Hierarchical) Sums of Squares
anova_glm(fit, type = 2)

# Type III (Partial) Sums of Squares with Effect Sizes and 90% Confidence Intervals
anova_glm(fit, type = 3, effectsize = TRUE, alpha = 0.10)
```

### 3. Linear Contrasts & Parameter Estimates

```r
# Test difference between Wool A and Wool B using named list
contrast_glm(fit, list(wool = c(1, -1)), name = "Wool A vs B")

# Estimate contrast with standard error, t-test, and 95% confidence interval
estimate_glm(fit, list(wool = c(1, -1)), name = "Wool A - B Difference")

# Multi-DF contrast testing tension main effect (L, M, H)
c_tension <- list(
  "L vs M" = list(tension = c(1, -1, 0)),
  "M vs H" = list(tension = c(0, 1, -1))
)
contrast_glm(fit, c_tension, name = "Tension Effect (2 DF)")
```

### 4. Least-Squares Means & Pairwise Comparisons

```r
# Compute LS-means for main effect
lsmeans_glm(fit, effect = "wool")

# Compute LS-means for interaction with Tukey-Kramer adjustment
lsm_inter <- lsmeans_glm(fit, effect = "wool:tension", adjust = "tukey")
print(lsm_inter)
```

### 5. Multivariate Analysis of Variance (MANOVA)

```r
data(iris)

# Fit multivariate response model
fit_m <- stat_glm(cbind(Sepal.Length, Sepal.Width, Petal.Length, Petal.Width) ~ Species, data = iris)

# Evaluate multivariate statistics (Wilks, Pillai, Hotelling-Lawley, Roy)
manova_glm(fit_m)
```

---

## Validation & Benchmark Consistency

`statGLM` includes an extensive unit test suite (`tests/testthat/`) verifying numerical identity against:
- Standard R `lm()` and `aov()` for baseline fits and Type I sums of squares.
- `car::Anova(..., type = 2)` for hierarchical Type II sums of squares.
- `car::Anova(..., type = 3)` for orthogonal Type III sums of squares.
- `car::linearHypothesis()` for general linear contrasts.
- Base R `manova()` for multivariate test statistics and Rao $F$-approximations.
- Published benchmark datasets from classical GLM literature, replicating outputs to $10^{-6}$ precision.

Run the test suite locally with:

```r
devtools::test()
```

---

## License

This package is licensed under the [MIT License](LICENSE).
