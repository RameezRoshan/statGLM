# statGLM 0.2.0

* Implemented comprehensive, modern, high-definition, publication-ready plotting capabilities matching all SAS `PROC GLM` ODS Graphics specifications:
  - **`plot.stat_glm()`**: Master S3 plot method with intelligent model design auto-detection and flexible `which` dispatcher.
  - **`plot_diagnostics()` (`PLOTS=DIAGNOSTICS`)**: 8-panel fit diagnostic display (Residuals vs Predicted with LOESS trend, RStudent vs Predicted, Observed vs Predicted with 1:1 identity line, RStudent vs Leverage with cutoffs, Cook's D needle plot with $4/n$ threshold, Normal Q-Q plot, Residual histogram with normal and kernel densities, and Residual-Fit spread plot). Supports `label` and `unpack` options.
  - **`plot_residuals()` (`PLOTS=RESIDUALS`)**: Scatter plots of residuals vs continuous covariates with optional LOESS smoothers and confidence ribbons.
  - **`plot_fit()` (`PLOTS=FITPLOT`)**: Continuous predictor fit plot with observed data, fitted curve, confidence limits for the mean (CLM), and prediction limits (CLI).
  - **`plot_contour()` (`PLOTS=CONTOURFIT`)**: 2D continuous predictor response surface contour plot and heatmap overlaid with observed data points (`outlinegradient`, `gradient`, `outline`, `none`).
  - **`plot_ancova()` (`PLOTS=ANCOVAPLOT`)**: Analysis of covariance plot showing fitted regression lines and confidence bands by factor level.
  - **`plot_boxplot()` (`PLOTS=BOXPLOT`)**: Grouped box plots with jittered observations and mean diamonds with 95% CI error bars.
  - **`plot_interaction()` (`PLOTS=INTPLOT`)**: Two-way ANOVA interaction plot with cell means, connecting lines, and confidence limits.
  - **`plot_nested()` (`PLOTS=NESTPLOT`)**: Nested classification box plots with clustered factor levels.
  - **`plot.lsmeans_stat_glm()`**: Master S3 plot method for LS-means objects.
  - **`plot_lsmeans()` (`PLOTS=MEANPLOT`)**: Least-squares means with error bars, confidence bands (`clband`), connecting lines (`connect`), and sorting (`ascending`, `descending`).
  - **`plot_diff()` (`PLOTS=DIFFPLOT`)**: Hsu's (1996) pairwise diffogram (mean-mean scatter plot) displaying pairwise comparison segments against the diagonal line of equality $y = x$.
  - **`plot_control()` (`PLOTS=CONTROLPLOT`)**: Dunnett-style comparisons against a designated reference control group.
  - **`plot_anom()` (`PLOTS=ANOMPLOT`)**: Analysis of Means display comparing LS-means against the grand mean with decision limits (UDL / LDL).
  - **`plot_lines()` (`LINES` option)**: Means comparison lines display showing non-significant subset groupings.
  - **`theme_statglm()`**: Clean, modern publication-ready ggplot2 theme with crisp borders and subtle gridlines.
  - **`save_glm_plot()`**: High-definition export helper supporting PDF, SVG, PNG, TIFF at 300+ DPI.
* Enhanced `predict.stat_glm()` to support `se.fit = TRUE` and `interval = c("none", "confidence", "prediction")`.

# statGLM 0.1.0

* Initial release of `statGLM`.
* Implemented clean-room general linear models with cell-indicator parameterization (`stat_glm`).
* Implemented reflexive generalized inverses via Goodnight symmetric SWEEP operator (`ginv_g2`, `sweep_matrix`) and Moore-Penrose pseudo-inverses via SVD (`ginv_svd`).
* Implemented estimability verification for linear hypotheses (`is_estimable`).
* Implemented Type I (sequential), Type II (hierarchical), and Type III (orthogonal/partial) sums of squares (`anova_glm`).
* Added effect sizes (partial/semipartial eta-squared and omega-squared) and noncentrality parameter confidence intervals via noncentral F inversion.
* Implemented general linear contrasts (`contrast_glm`) and estimable linear parameter functions with t-tests and confidence intervals (`estimate_glm`).
* Implemented least-squares means (population marginal means) and pairwise differences with multiple comparison adjustments (Tukey-Kramer, Bonferroni, Scheffe, Sidak) (`lsmeans_glm`).
* Implemented multivariate analysis of variance (`manova_glm`) supporting Wilks' Lambda, Pillai's Trace, Hotelling-Lawley Trace, and Roy's Greatest Root with exact and Rao approximate F-tests.
* Implemented SAS `PROC GLM` `RANDOM` statement functionality (`random_glm`) computing Expected Mean Squares (EMS) coefficients via Goodnight & Speed (1978), quadratic forms in fixed effects (Q option), and mixed model hypothesis tests with exact denominators and Satterthwaite pseudo-F approximations (TEST option). Seamlessly integrated into `stat_glm(..., random = ...)` and `anova_glm(..., random = ...)`.
* Comprehensive test suite comparing numerical outputs against base R (`lm`, `aov`, `manova`), `car::Anova`, and SAS PROC GLM published benchmarks.
