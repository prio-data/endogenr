#' endogenr: A Dynamic Endogenous Simulator of Statistical Models
#'
#' `endogenr` dynamically simulates a system of statistical / mathematical
#' models. The pipeline is three stages: [setup_system()] validates the panel
#' and the system, builds the dependency DAG, computes a topological execution
#' order, and materialises the `(unit x time x inner_sims)` grid;
#' [fit_system()] produces `nsim` coefficient draws (parameter / model
#' uncertainty); [simulate_system()] runs the dynamic loop in topological order,
#' drawing one predictive (innovation) realisation per `(unit, time, sim)` row.
#' The ensemble over `nsim x inner_sims` trajectories is the intended predictive
#' distribution.
#'
#' @section Conceptual contract:
#' Each `(unit, sim)` series **is** one coherent trajectory: `predict.*` reads
#' every lag from the same `(unit, sim)` series, so multi-step error propagation
#' (regressors that are themselves simulated become random) is captured by Monte
#' Carlo over `inner_sims`. The row-expansion architecture and the scheduling
#' are the load-bearing, well-tested core.
#'
#' @section Known issues - unclear semantics:
#' \itemize{
#'   \item **What the ensemble estimates.** `nsim` (parameter) x `inner_sims`
#'     (innovation) reads as a clean nested Monte Carlo *only if* each family
#'     separates the two. It does not (see below), so the ensemble is **not** a
#'     calibrated posterior-predictive of a single coherent generative object.
#'     Treat it as a nested-MC predictive distribution with per-family draws.
#'   \item **`nsim` vs `inner_sims` differ by family.** For `linear`/`glm`/
#'     `heterolm` with `min_window`, parameter uncertainty enters across draws.
#'     For `parametric_distribution`, `univariate_fable`, and `exogen` the model
#'     is fit **once and shared across all `nsim` draws**, so those components
#'     carry **zero parameter uncertainty** regardless of `nsim`. This asymmetry
#'     is real and currently undocumented at the call site.
#' }
#'
#' @section Known issues - arbitrary constants / heuristics:
#' \itemize{
#'   \item `min_window` fallback `10L`, duplicated in the sliding-window
#'     [fit_system()], [forecast_coefficients()], and `.fit_window_grid()`
#'     (`R/systemsim.R`, `R/coef_forecast.R`).
#'   \item `max_years = 50` default in `decay_since_event()`/`time_since_event()`/
#'     `intensity_decay()` (`R/formula_helpers.R`). These transforms depend on
#'     the whole prior series, so `.required_history()` reports `Inf` for them.
#'   \item Winkler `level = 50` default in [get_accuracy()] (`R/systemsim.R`) - a
#'     50% interval is an unusual reporting default.
#'   \item `est_mb > 1000` memory-warning threshold in `prepare_simulation_data()`
#'     (`R/systemsim.R`).
#'   \item Hard-coded ribbon quantiles `c(.05, .5, .95)` in [plotsim()] and
#'     `.cf_forecast_one()` (`R/coef_forecast.R`); `.cf_forecast_one()` also
#'     requires at least 3 window anchors.
#' }
#'
#' @section Known issues - available generalisations:
#' \itemize{
#'   \item **Unified predictive-draw interface.** The per-family draw logic
#'     (`getpi`, `getpi_glm`, `predict.heterolm`'s `rnorm`,
#'     `.sample_from_fitdist`, fable's `generate`) could be a single
#'     `draw_predictive(model, newdata, n_param, n_innov)` contract that
#'     explicitly separates parameter from innovation uncertainty.
#'   \item **Pooled vs panel-heterogeneous estimation.** `linear`/`glm` are
#'     pooled OLS/MLE with no unit effects unless the user writes `factor(unit)`.
#'     First-class fixed/random effects could be offered.
#' }
#'
#' @section Known issues - where brute force should become statistics:
#' \itemize{
#'   \item **Random-window refit (`min_window`)** "confounds recency with sample
#'     size" (`R/coef_forecast.R`). Principled replacements are prototyped in
#'     [forecast_coefficients()]: state-space time-varying-parameter models
#'     (dlm / KFAS) or a fixed-window residual/parametric coefficient bootstrap.
#'     Wiring the TVP path into the simulator is recommended.
#'   \item **Independent per-row innovations.** Each `(unit, time, sim)` row
#'     draws an independent innovation, so the engine ignores **contemporaneous
#'     cross-unit** correlation (common shocks) and **cross-equation** residual
#'     correlation. The statistical object is a SUR / seemingly-unrelated
#'     residual covariance (optionally spatially structured); joint innovation
#'     draws from an estimated covariance would capture it. The gap is
#'     demonstrated in `tests/testthat/test-uncertainty.R`.
#'   \item **Distribution-parameter uncertainty ignored.** `parametric_distribution`
#'     draws from the point-estimate MLE, and the shared independent fits carry
#'     no parameter uncertainty; both could propagate it via the MLE asymptotic
#'     vcov or a parametric bootstrap.
#'   \item **Double/triple-counted parameter uncertainty** for `linear` + `boot`
#'     + `min_window`: the random window perturbs coefficients, the residual
#'     bootstrap perturbs them again, and `predict.lm(se.fit)` adds parameter
#'     variance a third time on top of the residual scale. This is very likely
#'     over-dispersed; `tests/testthat/test-uncertainty.R` measures (rather than
#'     silently re-architects) the effect.
#'   \item **GLM proportion draws.** For `binomial`/`quasibinomial` outcomes with
#'     no trial count in the grid, `getpi_glm()` samples a Beta with mean `mu`
#'     and a precision derived from the dispersion. The Beta variance is capped
#'     at the Bernoulli value `mu(1 - mu)`, so quasibinomial overdispersion
#'     beyond that cannot be represented without a trial count and collapses to
#'     ~Bernoulli draws. Unsupported families fall back to a link-scale
#'     (parameter-only) draw with a one-time warning.
#' }
#'
#' @keywords internal
"_PACKAGE"
NULL
