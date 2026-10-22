# statGLM: Statistical General Linear Models in R

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
- **Random Effects & Expected Mean Squares (`random_glm`)**: Complete implementation of the SAS GLM `RANDOM` statement. Computes Expected Mean Squares (EMS) coefficients (Goodnight & Speed 1978), quadratic forms in fixed effects (`Q` option), and mixed model hypothesis tests using exact error denominators or Satterthwaite pseudo-$F$ approximations (`TEST` option). Seamlessly integrated into `stat_glm(..., random = ...)` and `anova_glm(..., random = ...)`.
- **Modern Publication-Ready Graphics (SAS PROC GLM ODS Graphics)**:
  - Master S3 `plot()` methods for `stat_glm` and `lsmeans_stat_glm` with intelligent model design auto-detection.
  - Complete suite of 14 SAS ODS Graphics displays: Fit diagnostics panel (`DIAGNOSTICS`), residuals vs covariates (`RESIDUALS`), continuous predictor fit curves with CLM/CLI bands (`FITPLOT`), 2D response surface contours (`CONTOURFIT`), ANCOVA plots with confidence limits (`ANCOVAPLOT`), grouped box plots with jitter and mean diamonds (`BOXPLOT`), two-way interaction plots (`INTPLOT`), nested box plots (`NESTPLOT`), Hsu's pairwise diffograms (`DIFFPLOT`), LS-means confidence intervals/bands (`MEANPLOT`), Dunnett control comparison forest plots (`CONTROLPLOT`), Analysis of Means (`ANOMPLOT`), and subset comparison lines (`LINES`).
  - Native high-definition vector graphics using `ggplot2`, styled with the clean academic `theme_statglm()` and exported via `save_glm_plot()`.

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

### 5. Expected Mean Squares & Satterthwaite Approximation (`RANDOM` Statement)

For a mixed linear model with fixed effects $\beta_0$ and random components $\beta_1, \ldots, \beta_k \sim \mathcal{N}(0, \sigma_j^2 I)$, the expected value of any quadratic form $\mathrm{SS}_L = (L b)' [L (X'X)^- L']^{-1} (L b)$ is derived via the Goodnight & Speed (1978) formulation:

$$\mathbb{E}(\mathrm{SS}_L) = \beta_0' C_0' C_0 \beta_0 + \sum_{j=1}^k \mathrm{tr}(C_j' C_j) \sigma_j^2 + \mathrm{rank}(L) \sigma^2$$

where $C = M L = (C_0 \mid C_1 \mid \cdots \mid C_k)$, and $M$ is the inverse Cholesky factor satisfying $M [L (X'X)^- L'] M' = I$. Dividing by $\mathrm{rank}(L)$ yields the Expected Mean Square $\mathbb{E}(\mathrm{MS}_L)$.

When testing hypotheses, if no single mean square isolates the required expectation under $H_0$, a synthetic linear combination $\mathrm{MS}_{\mathrm{den}} = \sum_i a_i \mathrm{MS}_i$ is constructed, and denominator degrees of freedom are evaluated via Satterthwaite's (1946) approximation:

$$\mathrm{DF}_{\mathrm{den}} = \frac{\left( \sum_i a_i \mathrm{MS}_i \right)^2}{\sum_i \frac{(a_i \mathrm{MS}_i)^2}{\mathrm{DF}_i}}$$

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

# Load simulated fully factorial aquaculture harvest dataset (1,500 observations)
data(harv)

# Fit general linear model (Factorial ANCOVA with initial weight covariate)
fit <- stat_glm(bw ~ bw0 + stock * culturetype + sex, data = harv)

# Display model fit and overall ANOVA
summary(fit)
```

### 2. Sequential, Hierarchical, and Orthogonal ANOVA

```r
# Type I (Sequential) Sums of Squares
anova_glm(fit, type = 1)

# Type II (Hierarchical) Sums of Squares
anova_glm(fit, type = 2)

# Type III (Partial) Sums of Squares with Effect Sizes and 95% Confidence Intervals
anova_glm(fit, type = 3, effectsize = TRUE, alpha = 0.05)
```

### 3. Linear Contrasts & Parameter Estimates

```r
# Fit main-effects model for single-DF linear hypotheses
fit_main <- stat_glm(bw ~ bw0 + stock + culturetype + sex, data = harv)

# Test genetic difference between Stock 1 and Stock 2
contrast_glm(fit_main, list(stock = c(1, -1)), name = "Stock 1 vs 2")

# Estimate contrast with standard error, t-test, and 95% confidence interval
estimate_glm(fit_main, list(stock = c(1, -1)), name = "Stock 1 - 2 Difference")

# Contrast testing sex effect (Male vs Female)
contrast_glm(fit_main, list(sex = c(1, -1)), name = "Male vs Female")
```

### 4. Least-Squares Means & Pairwise Comparisons

```r
# Compute LS-means for main genetic stock effect
lsmeans_glm(fit, effect = "stock")

# Compute LS-means for stock x culturetype interaction with Tukey-Kramer adjustment
lsm_inter <- lsmeans_glm(fit, effect = "stock:culturetype", adjust = "tukey")
print(lsm_inter)
```

### 5. Multivariate Analysis of Variance (MANOVA)

```r
# Fit multivariate response model (initial stocking weight bw0 and final harvest weight bw)
fit_m <- stat_glm(cbind(bw0, bw) ~ stock * culturetype, data = harv)

# Evaluate multivariate statistics (Wilks' Lambda, Pillai's Trace, Hotelling-Lawley, Roy's Root)
manova_glm(fit_m)
```

### 6. Mixed Models & Expected Mean Squares (`random_glm`)

Declare random effects using standard formulas or character vectors (matching SAS GLM `RANDOM effects / TEST;`):

```r
# Aquaculture selective breeding mixed model: Sires as random genetic family effect
fit_rnd <- stat_glm(bw ~ bw0 + stock + culturetype + sex + sire, data = harv)

# Compute Expected Mean Squares (EMS) and Satterthwaite mixed-model ANOVA tests
rnd <- random_glm(fit_rnd, effects = ~ sire, test = TRUE)
print(rnd)

# Or seamlessly via anova_glm()
anova_glm(fit_rnd, type = 3, random = ~ sire)
```

### 7. Modern Publication-Ready Visualizations

`statGLM` provides full graphical fidelity matching SAS `PROC GLM` ODS Graphics specifications, rendered via `ggplot2` with publication-standard typography and vector resolution:

```r
library(ggplot2)

# --- A. Master plot.stat_glm() with Intelligent Design Auto-Detection ---
# Auto-detects two-way factorial design -> produces interaction plot (PLOTS=INTPLOT)
plot(fit)

# Auto-detects one-way layout -> produces box plot (PLOTS=BOXPLOT)
fit_1way <- stat_glm(bw ~ stock, data = harv)
plot(fit_1way)

# Auto-detects single continuous predictor -> produces fit curve with CLM & CLI bands (PLOTS=FITPLOT)
fit_reg <- stat_glm(bw ~ bw0, data = harv)
plot(fit_reg)

# --- B. Comprehensive 8-Panel Fit Diagnostics (PLOTS=DIAGNOSTICS) ---
# Displays: Residual vs Fitted (LOESS), RStudent vs Fitted, Observed vs Predicted (1:1),
# RStudent vs Leverage, Cook's D needles (4/n cutoff), Normal Q-Q, Histogram, & RF-spread plot
plot(fit, which = "diagnostics", label = TRUE)

# Or extract as individual ggplot objects
diag_plots <- plot(fit, which = "diagnostics", unpack = TRUE)

# --- C. Residuals vs Initial Weight Covariate (PLOTS=RESIDUALS) ---
plot(fit, which = "residuals", smooth = TRUE)

# --- D. Analysis of Covariance (ANCOVA) Plots (PLOTS=ANCOVAPLOT) ---
fit_anc <- stat_glm(bw ~ stock + bw0, data = harv)
plot(fit_anc, which = "ancovaplot", clm = TRUE)

# --- E. 2D Response Surface Contour Plots (PLOTS=CONTOURFIT) ---
fit_surf <- stat_glm(bw ~ bw0 + age, data = harv)
plot(fit_surf, which = "contourfit", obs = "outlinegradient")

# --- F. LS-Means Multiple Comparisons & Hsu's Diffogram (PLOTS=DIFFPLOT) ---
lsm <- lsmeans_glm(fit, effect = "stock:culturetype", adjust = "tukey")

# 1. LS-Means with 95% Confidence Intervals & connecting lines (PLOTS=MEANPLOT)
plot(lsm, which = "meanplot", connect = TRUE)

# 2. Hsu's (1996) Pairwise Diffogram (Mean-Mean Scatter Plot)
# Segments crossing the 45-degree diagonal y = x are not statistically significant!
plot(lsm, which = "diffplot", center = TRUE)

# 3. Dunnett-Style Control Comparisons (PLOTS=CONTROLPLOT)
plot(lsm, which = "controlplot", control = "1:1")

# 4. Analysis of Means Display (PLOTS=ANOMPLOT)
# Grand mean central line with Upper and Lower Decision Limits (UDL / LDL)
plot(lsm, which = "anomplot")

# --- G. High-Definition Publication Export ---
# Save crisp 300+ DPI graphics ready for academic journals
p_diff <- plot(lsm, which = "diffplot")
save_glm_plot(p_diff, "figure_diffogram.png", width = 7, height = 6, dpi = 300)
save_glm_plot(p_diff, "figure_diffogram.pdf", width = 7, height = 6)
```

---

## Validation & Benchmark Consistency

`statGLM` includes an extensive unit test suite (`tests/testthat/`) verifying numerical identity against:
- Standard R `lm()` and `aov()` for baseline fits and Type I sums of squares.
- `car::Anova(..., type = 2)` for hierarchical Type II sums of squares.
- `car::Anova(..., type = 3)` for orthogonal Type III sums of squares.
- `car::linearHypothesis()` for general linear contrasts.
- Base R `manova()` for multivariate test statistics and Rao $F$-approximations.
- SAS PROC GLM Example 53.8 (Milliken & Johnson 1984 unbalanced mixed model) reproducing Expected Mean Squares coefficients and Satterthwaite pseudo-$F$ test statistics to exact decimal precision.
- Published benchmark datasets from classical GLM literature, replicating outputs to $10^{-6}$ precision.

Run the test suite locally with:

```r
devtools::test()
```

---

## License

This package is licensed under the [MIT License](LICENSE).
