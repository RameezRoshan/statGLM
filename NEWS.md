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
* Comprehensive test suite comparing numerical outputs against base R (`lm`, `aov`, `manova`) and `car::Anova`.
