# endogenr 0.1.0.9000

## New features

- **`build_model("cross_section", outcome ~ I(expr))` — per-period
  cross-sectional aggregates.** `expr` is evaluated across all units at each
  simulated time step, separately within each simulation draw (e.g. a
  simulated productivity frontier
  `front ~ I(quantile(lgdppc_l1, 0.9, na.rm = TRUE))`). A scalar result is
  copied to every unit; a result with one value per unit (e.g. `rank(x)`) is
  kept per unit; any other length errors. `deterministic` evaluates per unit,
  so such aggregates previously failed with "variable lengths differ".
  Time-series functions (`lag`, `diff`, rolling, cumulative) are rejected at
  `build_model()` time; reference an already-lagged column instead.
  `draw_predictive()` errors for `cross_section` like `deterministic`.

- **`draw_predictive()` — unified predictive-draw interface.** One exported
  S3 generic, `draw_predictive(model, newdata, n_param, n_innov, ...)`,
  replaces the per-family draw logic (`getpi`, `getpi_glm`,
  `predict.heterolm`'s inline `rnorm`, point-estimate `.sample_from_fitdist`
  draws) and explicitly separates **parameter** uncertainty (`n_param`;
  `0` = condition on point estimates) from **innovation** uncertainty
  (`n_innov`; `0` = conditional mean, subsuming `what = "expectation"`).
  Returns a `nrow(newdata) x max(n_param,1)*max(n_innov,1)` matrix; column
  `(j-1)*K + k` is parameter draw `j`, innovation draw `k`. Methods for
  `linear`, `glm`, `heterolm`, `gamlss`, `glmmTMB`, `parametric_distribution`,
  and `univariate_fable` (whole coherent `generate()` paths per column;
  takes `ctx` and `horizon`); `deterministic`/`exogen`/`spatial_lag` error
  (no stochastic predictive distribution). The linear/glm/glmmTMB methods
  take `param_scope = c("draw", "row")`: `"draw"` shares one parameter
  deviate per column (coherent replicates for analysis), `"row"` draws an
  independent deviate per row — the engine's row-expansion convention, which
  `predict.*` now uses internally so **ensemble statistics are unchanged**
  (glm/heterolm/gamlss/glmmTMB engine draws are RNG-stream-identical;
  `linear` draws come from the identical predictive t distribution via a
  chi-square/normal mixture, but consume a different RNG stream under a
  fixed seed). `heterolm`/`gamlss`/`univariate_fable` carry no
  parameter-uncertainty path and warn once per session when `n_param > 0`.
  Internally, `getpi`/`get_sepi` are replaced by `.lm_predictive_draws()`
  and `getpi_glm` by `.glm_predictive_draws()` (both keep the documented
  per-row vs per-column deviate rationale).

- **`build_model("parametric_distribution", …, param_uncertainty = TRUE)` —
  opt-in MLE parameter uncertainty.** Draws one parameter vector per
  simulation from the fitted distribution's asymptotic MVN(estimate, vcov)
  (shared across that simulation's rows and times), so `fitdist()` models no
  longer have to carry zero parameter uncertainty. Invalid draws (e.g.
  negative `df` leaving the parameter space) fall back to the point
  estimates with a one-time warning. Default `FALSE` preserves the previous
  point-estimate behavior exactly.

- **Custom distributions work without attaching endogenr.** The
  location-scale Student-t functions `dt_ls`/`pt_ls`/`qt_ls`/`rt_ls` are
  exported again (reversing the 0.1.0.9000 internalization):
  `fitdistrplus::fitdist()` resolves `d<distname>`/`p<distname>` **by name
  from its own namespace chain** (`… -> globalenv -> search path`), so the
  functions must be reachable from the search path — the same visibility
  contract `actuar` uses. In addition, fitting wraps the `fitdist()` call in
  an internal shim that temporarily injects the package's d/p/q/r functions
  into the global environment (and always cleans up, restoring any shadowed
  binding), so `endogenr::build_model("parametric_distribution",
  distribution = "t_ls", …)` now works even when the package is loaded but
  **not attached** — previously this failed with "The dt_ls function must be
  defined".

- **Automatic, window-aware starting values for `fitdist()`.**
  `distribution = "t_ls"` no longer requires `start`: starting values are
  derived from the actual training-window data (median/MAD location-scale, a
  kurtosis-matched `df` clamped to `[2.5, 100]`). `start` may also be a
  `function(x)` returning a named list, evaluated by `fitdist()` on the
  training data at fit time — the correct tool under `run_experiments()` and
  sliding-window refits, where a precomputed list would refer to the wrong
  window. When `fitdist()` still cannot derive starting values, the error now
  says exactly what to pass to `build_model()`.

- **`build_model("glmmTMB", …)` — glmmTMB mixed-effects models.** Adds
  first-class support for `glmmTMB::glmmTMB()` models inside the endogenr
  dynamic simulation loop. Supports lme4-style random-effects bars
  (`(1 + lag(x) | group)`), glmmTMB covariance-structure wrappers
  (`ar1(times + 0 | group)`, `us`, `exp`, `mat`, …), a separate
  `dispformula` (per-row heteroscedasticity / overdispersion), and a
  `ziformula` for zero-inflation. Families with a response-scale draw:
  `gaussian`, `poisson`, `binomial`, `Gamma`, `nbinom1`, `nbinom2`,
  `beta`, `betabinomial`, `t`, `lognormal`, `skewnormal`,
  `truncated_poisson`, `truncated_nbinom1`, `truncated_nbinom2`, and
  `tweedie` (the last requires the `tweedie` package, else it falls back
  to the mean); other unsupported families fall back to the
  conditional mean with a one-time warning. Grouping factors and cov-struct
  coordinates are ordinary predictors: like any predictor they must be
  produced by some model — add an `exogen` (e.g.
  `build_model("exogen", formula = ~region)`) to carry the grouping column
  into the forecast horizon, or group by a panel key. Setup errors if a
  grouping column has no producer. Temporal covariance structures
  (`ar1`/`ou`/…) are forecast correctly multi-step: at each forecast step the
  whole forecast-so-far block is predicted in one call so glmmTMB applies the
  proper `phi^k` correlation decay, and the dispersion prediction passes
  `allow.new.levels` so the predictive-interval scale is correct for
  covariance-structure models. The covariance coordinate must be carried into
  the horizon (e.g. via an `exogen`) and be contiguous and unit-spaced.
  Requires the `glmmTMB` CRAN package (`install.packages("glmmTMB")`).

  ```r
  build_model("glmmTMB",
              formula     = gdppc_grwt ~ lag(dem) + lag(log(gdppc)) +
                              (1 + lag(dem) | region),
              dispformula = ~ lag(dem),
              family      = stats::gaussian())
  ```


- **`build_model("gamlss", …)` — GAMLSS distributional-regression models.**
  Adds support for `gamlss::gamlss()` models (Generalized Additive Models for
  Location, Scale and Shape). Each of the four distribution parameters (`mu`,
  `sigma`, `nu`, `tau`) gets its own formula, enabling modelling of
  heteroscedasticity, skewness, and kurtosis as functions of covariates.
  Supports all `gamlss.dist` families (Normal, BCT, BCCG, Gamma, …) with
  automatic draw via `r<FAMILY>()`. Smoother terms (`pb()`, `cs()`, `lo()`)
  and grouping terms (`random()`, `ra()`, `re()`) are handled correctly: the
  dependency graph is built from a smoother/bar-free representation. Grouping
  factors are ordinary predictors that must be produced by some model — add an
  `exogen` (e.g. `build_model("exogen", formula = ~region)`) to carry the
  grouping column into the forecast horizon, or group by a panel key. Setup
  errors if a grouping column has no producer. Requires the `gamlss` CRAN
  package (`install.packages("gamlss")`).

  ```r
  build_model("gamlss",
              formula       = gdppc_grwt ~ pb(lag(dem)) + random(region),
              sigma.formula = ~ lag(dem),
              family        = gamlss.dist::NO())
  ```

- **`build_model(…, bounds = c(lower, upper))` — per-model output bounds.**
  Any model's simulated outcome can be constrained to `[lower, upper]`: the
  dynamic simulator clamps the predicted value at every forecast step before it
  feeds later steps, stabilizing autoregressive feedback (e.g.
  `bounds = c(-1, 1)` on a growth rate). A draw that is still non-finite after
  clamping is reset to a finite in-range value (the midpoint of two finite
  bounds, or the finite bound of a one-sided `c(0, Inf)` limit), so a bounded
  outcome never propagates `NaN`/`Inf`. This prevents the
  `(subscript) logical subscript too long` crash seen when a divergent
  autoregressive term (e.g. a gamlss `sigma.formula` containing
  `abs(lag(outcome))`) drove a predictor to `Inf`/`NaN`. Default `NULL` leaves
  all behavior unchanged.

- **`fe_ar()` — AR(1) time fixed effects with parameter uncertainty; new
  default.** `setup_param()` (and `simulate_system()` without
  `scenario_params`) now bakes time fixed effects (`factor(timevar)`) with
  `fe_ar()`: per `(outer draw, inner sim)` trajectory it draws the year-effect
  series once from `MVN(tau_hat, V_tau)`, fits an AR(1) by OLS (`rho` clamped
  to `[-1, 1]`), and simulates forward from the last training-year effect by
  residual bootstrap, stepping through any gap between a window's last year
  and the forecast origin. Falls back to `fe_resample()` when fewer than 3
  year effects are estimated. `fe_resample()` and `fe_converge()` (time FE)
  now also resample from each trajectory's MVN-perturbed effects. The effect
  covariance is recorded at fit time (`model$time_fe$vcov`; linear via
  `vcov()`, gamlss from the stored `mu.qr`) and carried as
  `time_fe$vcov_by_draw`; when unavailable, draws use the point estimates
  and warn once. Time-FE draws under a fixed seed differ from earlier
  versions.

## Bug fixes

- **FE scenario slicing no longer fails with `inner_sims = 1` or
  `horizon = 1`.** `.slice_scenario()` indexed the baked arrays without
  `drop = FALSE`, so time- and per-trajectory unit-FE offsets lost a
  dimension and simulation aborted with "incorrect number of dimensions".

- **Clearer error for `T`/`F` shorthand in formulas.** `T` and `F` are
  variables in R, not logical literals, so a formula such as
  `y ~ I(quantile(x, 0.9, na.rm = T))` was rejected by
  `validate_system_closure()` with an unexplained "missing from the input
  data: T". The error now appends a note telling you to write `TRUE`/`FALSE`.

- **Execution order no longer drops variables whose names contain `"lag_"`.**
  The lag-vertex strip in `get_execution_order()` matched the substring
  `"lag_"` anywhere in a vertex name, so an outcome named e.g. `flag_war`
  was misclassified as a lagged input, removed from the execution order, and
  its model silently never ran — the whole forecast degraded to `NA`. The
  match is now anchored to the `lag_` prefix.

- **User columns can no longer shadow internal filters (time column named
  `t`, etc.).** Every `data[data[[idx]] == t]`-style subset in the predict
  methods, the engine, validation, and scoring evaluated its expression
  inside data.table's column scope, where a user column with the same name
  as a function argument silently shadows it. A panel whose time column is
  literally named `t` degenerated the per-step filter to *all rows*, made
  the update-join overwrite observed training rows with `NA`, and cascaded
  `NA` through the forecast. All row subsets now go through a scope-safe
  helper (`.dt_rows()`, using data.table's `env=` substitution), so panel
  column names can never collide with internals. data.table (>= 1.14.2) is
  now required.

- **GLM predictive draws no longer share one parameter-uncertainty draw
  across all rows.** In the row-expansion path (one draw per row),
  `getpi_glm()` drew a *single* t-quantile and applied it to every
  `(unit, sim)` row at a time step, so the link-scale parameter-uncertainty
  component was perfectly correlated across units *and* inner simulations —
  under-dispersing the ensemble and breaking the documented `lm` parity for
  gaussian GLMs. It now draws an independent t per row, matching `getpi()`
  (linear) and the per-row normal draw in the glmmTMB predict method.

- **`diff()` in formulas works.** The registry advertised `diff` as a
  supported time-series function (and `.required_history()` composed its
  depth), but materialization crashed with a length-mismatch error because
  `base::diff()` returns `n - lag*differences` values. Materialization now
  uses a length-stable positional `diff` that pads the leading rows with
  `NA`, aligned row-for-row with the input series.

- **`shift()` and `lead()` in formulas work without attaching data.table or
  dplyr.** Both were in the time-series registry but had no binding in the
  materialization environment, failing with *"could not find function"*
  unless the user happened to have data.table/dplyr attached (and silently
  using pooled semantics was never possible — but resolution depended on the
  session). The materialization environment (and
  `inject_positional_lag()`, used by deterministic models) now binds
  positional `lag`/`lead`/`shift`/`diff` uniformly; `shift` delegates to
  `data.table::shift` (type-preserving, supports `n`/`fill`/`type`).
  `lag(x, -1)` now errors with a clear message pointing to `lead()`.

- **`shift(x)` is classified as a prior-period dependency.** The
  dependency-graph edge classifier only recognised `lag()` as lag context,
  so a model using `shift(x)` got a same-period `x -> y` edge — producing
  spurious "same-period dependency cycle" errors for legal mutually-lagged
  systems. Backward `shift()` (no `type=`, or `type = "lag"`) now sets lag
  context exactly like `lag()`; `shift(type = "lead")` stays same-period.

- **Intercept-only `gamlss`/`glmmTMB` formulas no longer fail closure
  validation.** A parameter formula with no predictors (e.g.
  `gdppc_grwt ~ 1`) built its dependency-graph RHS as the symbol `` `1` ``
  instead of the numeric intercept, so `all.vars()` reported a phantom
  predictor named `1`. `validate_system_closure()` then aborted with
  *"variables ... missing from the input data: 1"*. The empty-predictor
  graph RHS is now the numeric intercept `1`, so intercept-only mu/sigma/nu/tau
  (gamlss) and main/disp/zi (glmmTMB) formulas validate and simulate.

- **glmmTMB simulation is faster.** `predict.glmmTMB_endogenr()` no longer
  issues a redundant `predict(type = "disp")` rebuild at every forecast step
  when `dispformula` is trivial (`~1`, the default): the constant dispersion
  `sigma(fitted)` is cached at fit time (alongside the family name and
  link-inverse). This drops one of the two `predict.glmmTMB` rebuilds per step
  — about 40% of raw per-step prediction cost, ~20% less `simulate_system()`
  wall time in a typical gaussian configuration, and proportionally more for
  `ar1`/cov-struct models whose disp call re-predicts a growing block. Output
  is numerically identical; models with a non-trivial `dispformula` are
  unaffected.

- **`lag()` now preserves factor predictors.** A factor wrapped in `lag()`
  (e.g. `lag(conflict)`) was coerced to its integer codes during panel
  materialization, so it entered models as a single numeric term instead of
  expanding into dummy columns. The positional-lag shift is now type-preserving,
  so lagged factors (and `Date`/character columns) keep their class and levels
  across all model types (`linear`, `glm`, `heterolm`, `glmmTMB`, `gamlss`,
  `deterministic`).

## Performance

- **Fitted models no longer carry their constructor frame.** The regression
  families (`linear`, `glm`, `glmmTMB`, `gamlss`) stored the stage-2 fit
  step as a closure on the model object; its environment dragged the input
  training data, the materialized copy, and a circular reference to the
  model itself into every serialized draw shipped to parallel workers
  (measured: a stripped linear model serialized at ~3.5 MB, of which
  ~3.3 MB was the dead closure). The fit step is now a file-level function
  per family, cutting the per-draw payload by an order of magnitude under
  `future::multisession` plans.

- **Deterministic models evaluate only the history they need.**
  `predict.deterministic()` re-derived its evaluation formula and evaluated
  the model frame over the *entire* simulation grid at every forecast step
  (O(T x grid) per step). It now caches the prepared formula at build time
  and subsets to the composed `.required_history()` window like every other
  predict method.

- **`simulate_system()` assigns `.sim` arithmetically.** The unique draw id
  was built with `unique()` + a full-table `merge` (which also re-sorted the
  result by `(.id, sim)`); it is now computed in place as
  `(.id - 1) * inner_sims + sim`, preserving the `(unit, sim, time)` row
  order and skipping the merge copy. The numbering is unchanged.

- **`validate_panel()` contiguity check is a single grouped pass.** The
  per-unit scan re-filtered the whole table once per unit
  (O(units x rows)); it now runs one data.table by-group pass.

- **`predict.spatial_lag()` splits the multi-sim frame once** instead of
  re-scanning it per sim id (O(sims^2 x units) -> O(sims x units)), and
  `.score_draws()` computes both Winkler quantiles in one `quantile()`
  call.

## Maintainability

- Removed dead code: `.max_lag_depth()` (superseded by
  `.required_history()`) and `select_col_per_row()` (no callers), plus
  their tests. The duplicated positional-lag implementation and the
  ten copies of the `.pt#` alias-rename block now share single helpers
  (`.pt_bind_ts_fns()`, `.pt_apply_aliases()`).
  `get_coefficients()`/`plot_coefficients()` docs now state that `glmmTMB`
  and `gamlss` coefficients are included.

## Migration guide

The three-stage pipeline replaces the old two-stage API. The minimum
working translation:

```r
# Old (main)
sim <- setup_simulator(models, data, ...) |> simulate_endogenr(nsim = 50, ...)

# New (HEAD)
setup  <- setup_system(models, data, ...)     # validate + DAG + grid
fitted <- fit_system(setup, nsim = 50)        # store 50 coefficient draws
sim    <- simulate_system(fitted)             # predict-only, no fitting
```

`parallel` and `ncores` arguments are removed. Set a `future` plan before
fitting and wrap calls in `progressr::with_progress()` for a progress bar:

```r
future::plan(future::multisession, workers = 4)
progressr::with_progress({
  fitted <- fit_system(setup, nsim = 50)
  sim    <- simulate_system(fitted)
})
future::plan(future::sequential)
```

Model constructors are replaced by `build_model(type, ...)`:

```r
# Old
linearmodel(formula = y ~ lag(x), boot = "resid")

# New
build_model("linear", formula = y ~ lag(x), boot = "resid")
```

---

## Breaking changes

* **Three-stage pipeline.** The two-stage `setup_simulator()` +
  `simulate_endogenr()` API is replaced by `setup_system()` +
  `fit_system()` + `simulate_system()`. Fitting is now a distinct stage
  that stores `nsim` coefficient draws; `simulate_system()` is
  predict-only and derives `nsim` from the stored draws.
  `inner_simulation()` (the per-draw fitting loop) is removed.

* **Parallelism moved to the caller.** `simulate_endogenr(parallel,
  ncores)` is removed. The caller sets a `future` plan; both
  `fit_system()` and `simulate_system()` use `future.apply` internally.
  Wrap the pipeline in `progressr::with_progress()` for a progress bar.

* **`build_model()` replaces model constructors.** The eight
  constructor functions (`linearmodel()`, `glmmodel()`,
  `heterolmmodel()`, `exogenmodel()`, `parametric_distribution_model()`,
  `univariate_fable_model()`, `spatial_lag_model()`,
  `deterministicmodel()`) are now internal. Use
  `build_model(type, ...)` instead, where `type` is one of `"linear"`,
  `"glm"`, `"heterolm"`, `"exogen"`, `"parametric_distribution"`,
  `"univariate_fable"`, `"spatial_lag"`, `"deterministic"`. S3 fitting
  is dispatched via the new `fit_model()` generic.

* **De-exported helpers.** The following low-level functions are now
  internal (marked `@keywords internal`): `bootstraplm`, `bootstrapglm`,
  `getpi`, `getpi_glm`, `get_sepi`, `parse_formula`, `func_in_term`,
  `prepare_simulation_data`, `process_independent_models`,
  `process_dependent_models`, `pt_ls`, `qt_ls`, `rt_ls`,
  `select_col_per_row`, `update_dependency_graph`, `new_endogenmodel`,
  `create_distribution_object`, `create_panel_frame`,
  `fit_parametric_distribution_model`, `get_independent_models`.

* **`scoringRules` replaces `fabletools` for scoring.** CRPS,
  MAE, and Winkler scores in `get_accuracy()` / `get_lh_accuracy()` now
  use `scoringRules::crps_sample()` instead of
  `fabletools::accuracy()`. `fabletools` is no longer an Import;
  it remains in Suggests for the `univariate_fable` model family.

---

## New features

### Experiment grid

* **`run_experiments()`** runs the full `setup_system()` →
  `fit_system()` → `simulate_system()` pipeline across a Cartesian grid
  of model systems, training starts, forecast starts, and (optionally)
  window configurations, stacking results into a single `data.table`.
  Each experiment is stamped with `.experiment`, `model`, `train_start`,
  `test_start`, and `window_cfg` columns. Supports `parallel =
  "experiments"` (grid mapped over workers) or `parallel = "draws"`
  (single-experiment parallelism).

* **`vary_model(base, variants)`** builds a named list of full model
  systems by swapping specs into a base system by outcome variable,
  producing the exact shape accepted by `run_experiments()`.

* **`get_experiment_accuracy(experiment_results, outcome, truth)`**
  scores every experiment in a `run_experiments()` result with
  per-experiment `test_start` awareness, returning CRPS, MAE, and
  Winkler scores keyed by `.experiment`.

* **`run_experiments()` gains a `windows` axis via `window_config()`.**
  Pass a single `window_config()` to apply it uniformly, or a named list
  to cross window approach as a fourth grid axis alongside `models`,
  `train_start`, and `test_start`. Results carry a `window_cfg` column;
  the `experiments` attribute includes `window_cfg`, `window`, and
  `window_policy` columns. The default (`windows = NULL`) reproduces
  historical behaviour (random windows, latest policy).

### Coefficient-trajectory forecasting

* **`forecast_coefficients(object, ...)`** estimates a deterministic
  coefficient path for each `linear`/`glm` spec by refitting on a grid
  of expanding or rolling training windows, then projects every
  coefficient forward over the forecast horizon as a Gaussian random walk
  (optionally with drift). The forward distribution fans out from the
  forecast origin's estimation uncertainty (vcov) and the empirical
  evolution covariance of the path.

* **`plot_coefficient_forecast(object, ...)`** plots the observed
  coefficient path plus forecast band for each term, faceted by outcome.

### Long-horizon reduced-form benchmark (redesigned)

* **`lead_horizon(x, h)`** leads a vector `h` steps ahead within a
  group (the long-horizon analogue of `lag()`). Used on the formula LHS
  to declare the forecasting target.

* **`setup_long_horizon(data, formulas, horizons, ...)`** fits pooled
  OLS/bootstrap models for every `(formula, horizon)` combination using a
  two-stage pooled panel design (Stage 1: within-unit time-series
  materialisation; Stage 2: pooled `lm()` fit, so `poly()`, `bs()`,
  `factor()`, interactions, and `-1` are all handled correctly). Returns
  a `lh_setup` object.

* **`forecast_long_horizon(lh_setup, data, ...)`** generates
  probabilistic forecasts (`nsim` bootstrap refits × `inner_sims`
  predictive draws per unit and horizon).

* **`get_lh_accuracy(lh_forecasts, truth, ...)`** scores long-horizon
  draws against observed truth with CRPS, MAE, and Winkler scores,
  per-variant outcome-aware. Supports `"native"` and `"model"` scoring
  scales with optional inverse transform.

* **`cv_long_horizon(data, formulas, horizons, test_starts, ...)`**
  runs the full pipeline across multiple training folds and returns
  accuracy scores averaged across folds.

* **`compare_approaches(lh_accuracy, sim_accuracy)`** stacks long-horizon
  and dynamic-simulation accuracy results for side-by-side comparison.

* Helper: **`create_lh_data()`** and **`fit_lh_model()`** are exported
  for users who want finer control over the two-stage fitting.

### Panel and system validation

* **`validate_panel(data, ctx, test_start, model_outcomes = NULL)`**
  checks integer time steps, contiguous series per unit, origin coverage
  (at least one unit must have data at `test_start - 1`), complete
  initial state for units present at the origin, and sort order. Throws
  informative errors.

* **Unbalanced panels are supported.** Units may enter late or exit
  early; each unit's series must still be contiguous. Units without
  data at `test_start - 1` are used for training only and are excluded
  from the simulation grid (with a message). `setup_system()` errors
  when a unit present at the forecast origin enters too late to supply
  the history depth the model formulas require.

* **`validate_system_closure(models, data_columns)`** checks that every
  RHS predictor is produced by some model, no two models share an
  outcome, and every formula-referenced variable exists as a data column.

* Both validations run automatically inside `setup_system()`.

### Coefficient introspection

* **`get_coefficients(fitted_system)`** returns a tidy `data.table` of
  coefficient draws (mean, SD, per-draw columns) for all `linear`/`glm`
  specs in a fitted system.

* **`plot_coefficients(fitted_system, ...)`** plots per-term coefficient
  distributions as violin plots, faceted by outcome.

* **`plot_estimates(fitted_system, ...)`** plots coefficient estimates
  with uncertainty ribbons.

### Sliding-window fitting

* **`fit_system()` gains `window`, `width`, and `step` arguments** for
  rolling (`"rolling"`) and expanding (`"expanding"`) window fitting in
  addition to the existing random-window bootstrap (`"random"`).
  Internal helpers `.fit_window_grid()` and `.window_weight_matrix()`
  manage the grid of window-end anchors and per-horizon weighting.

* **`simulate_system()` gains `window_policy`, `decay`, and `weights`
  arguments** controlling how multiple window fits are blended at each
  forecast step: `"latest"` (default), `"equal"`, or `"decay"`.

### `panel_context` object

* **`panel_context(unit, time, sim = NULL)`** is a lightweight metadata
  carrier that stores panel dimension names separately from the data,
  preventing silent metadata loss on data operations.

* Accessors: `ctx_unit(ctx)`, `ctx_time(ctx)`, `ctx_sim(ctx)`,
  `ctx_keys(ctx)`.

---

## Dependencies

* **Imports gained:** `broom`, `future`, `future.apply`, `igraph`,
  `janitor`, `rlang`, `scoringRules`. `broom` moved from Suggests to
  Imports.
* **Imports dropped:** `R6` (removed entirely).
* **Suggests gained:** `fable`, `fitdistrplus`, `tidyr`, `tsibble`.
* **Suggests dropped:** `pbapply`, `conflicted`.

---

## Bug fixes

* **Bootstrap fits now estimate on the same sample as a plain
  regression.** `bootstraplm()`/`bootstrapglm()` previously ran
  `na.omit()` on the entire materialized training table, so an `NA` in
  any unrelated panel column silently dropped that row from every
  `boot = "resid"`/`"wild"` fit — but not from a plain `lm()`/`glm()`
  on the same window. Regression fits are now restricted to the
  model's own columns before fitting, so the bootstrap estimation
  sample equals `lm()`'s complete cases on the model terms. Previously
  this could materially attenuate coefficients (e.g. conflict effects)
  whenever unrelated columns had missing values.

* **Bootstrap refits now assign the resampled response to a dedicated
  `.boot_y` column.** Previously the resampled response was written to
  `as.character()` of the formula LHS; a call LHS (e.g. `log(y)`) made
  the refit error out or silently ignore the resampled response,
  collapsing every "draw" onto the base fit (fake zero parameter
  uncertainty). Poisson-family refits on the intentionally continuous
  link-scale reconstruction no longer emit `non-integer` warnings
  (muffled locally inside `bootstrapglm()`).

* **`get_train_window()` now always returns `end <= test_start - 1`**
  (the documented contract; previously a zero random decrement labelled
  the window as ending at `test_start`) **and handles `min_window`
  equal to the full training range** (previously `sample.int(0)`
  errored). Random-window RNG streams shift accordingly.

* **`poly(dem, 2)` and other data-dependent design bases now work inside
  panel model formulas.** Previously, `linearmodel()`, `glmmodel()`, and
  `heterolmmodel()` evaluated the full formula per group via
  `model.frame()`, so `poly(dem, 2)` (and any basis that needs ≥
  degree+1 unique values) would error with `'degree' must be less than
  number of unique points` whenever a unit had fewer unique values than
  the degree — a guaranteed failure when `.build_mat_cache()` probed
  with only 2 rows per group. The root cause: data-dependent bases must
  be evaluated **pooled** (matching what a base `lm()` on the same data
  produces), not per unit. The fix is a two-stage approach already used
  by `longhorizon`: Stage 1 materialises only the true within-unit
  time-series sub-expressions (`lag`, `diff`, `rollmean`, …) per group;
  Stage 2 fits a pooled `lm`/`glm` on the rewritten formula, so
  `poly()`, `bs()`, `factor()`, interactions and `-1` are all handled by
  base R's formula machinery and reproduced coherently at predict time
  via `predict.lm()`'s stored `predvars`/`xlevels`.

* **Coherent predict-time basis for heterolm.** `predict.heterolm()` now
  reconstructs `model.matrix()` columns using the training-time `terms`
  object (with `predvars` and `xlevels`), ensuring polynomial/spline
  bases and factor contrasts match those from fitting exactly.

* **`spatial_lag_model()` now correctly handles island units.** The
  constructor normalises both the manual `integer(0)` representation and
  the `sfdep`/`spdep` `0L` sentinel to `0L`, and detects islands via a
  predicate robust to both encodings. Previously, `integer(0)` islands
  caused `sfdep::st_lag()` to error with `"zero length neighbour
  vector"`, and `0L` islands silently returned `0` instead of
  `island_default`.

* **`predict.linear()`, `predict.glm_endogenr()`, and
  `predict.heterolm()` now size the per-unit history window by
  `.required_history()`**, an AST walk that **composes** nested
  time-series depths instead of taking the maximum of individual depths
  (`.max_lag_depth()`). Composed transforms such as `lag(lag(x))`,
  `lag(rollmeanr(x, k), n)`, and cumulative transforms (`cumsum()`,
  `decay_since_event()`) are no longer silently truncated, which
  previously produced `NA` or wrong values in forecast cells. Cumulative
  and since-event transforms now correctly use the full per-unit history.
  The same sizing replaces the hard-coded 20-row history in
  `.predict_lh()`.

* **`predict.univariate_fable()` now draws coherent sample paths** with
  `fabletools::generate()` on the fitted mable (innovations simulated
  forward with the model's temporal structure), rather than stitching
  together independent per-horizon draws from the marginal `forecast()`
  distributions. Each `(unit, sim)` trajectory is a single
  internally-consistent path.

* **`getpi_glm()` adds response-scale dispersion** in addition to
  link-scale parameter uncertainty. A Gaussian GLM's prediction interval
  again matches the equivalent `lm` interval; Poisson/quasipoisson
  draws are non-negative integers with the right mean/variance; Gamma
  draws are positive with variance `dispersion * mu^2`;
  binomial/quasibinomial proportion outcomes are drawn from a Beta with
  mean `mu`. Unsupported families fall back to the previous
  link-scale-only draw with a one-time warning.

---

## Architecture

* `R/panel_transform.R` is now the **canonical panel design layer**. The
  time-series function registries (`lag`, rolling, cumulative, decay
  functions) are defined once as `.pt_ts_fns`, `.pt_roll_fns`,
  `.pt_cum_fns` and shared by both materialisation (`panel_materialize()`)
  and history-depth/edge extraction (`systemgraph.R`). The previous
  `systemgraph.R`-local `.rh_cum_fns` / `.rh_roll_fns` duplicates are
  removed.

* **Full `data.table` migration.** Simulation data structures throughout
  the pipeline use `data.table` operations (by-reference updates,
  `CJ()`-based grid construction, `setattr()` for metadata). The
  `panel_context` object separates panel metadata from the data carrier.

* **Naming contract.** Time-series synthetic columns (`.pt#` internal
  keys) are renamed to human-readable aliases — `janitor`-cleaned
  deparsed source expressions — before the pooled fit, so coefficient
  names from `build_model("linear")` / `build_model("glm")` match what
  base `lm()` would produce: `lag(x)` → `lag_x`,
  `lag(log(gdppc))` → `lag_log_gdppc`. The internal `.pt#` symbols
  never appear in `get_coefficients()` output or `$coefs`.

* **Formulas preserve interactions and full term structure.** Previously,
  `derive_naive_formula()` rebuilt the formula additively from
  materialised column names, silently dropping interaction terms
  (`g:x`, `g*x`). The fitted `linear` and `glm` models now estimate
  interactions identically to base `lm()`/`glm()`. For `heterolm`,
  interaction terms are pre-expanded via `model.matrix()` into flat
  predictor columns before fitting.

* **`factor(<panel-unit-id>)` in model formulas** (e.g. unit fixed
  effects via `y ~ factor(gwcode) + x`) now materialises correctly.
  Previously, the panel key column was absent from the per-group `.SD`
  inside `model.frame()`, causing `"object 'gwcode' not found"`. The new
  internal `.model_frame_by_group()` helper recycles the scalar
  group-key value into each per-group frame when the formula references
  it.

* **`build_model()` + `fit_model()` S3 dispatch.** Every model family
  now has a spec class (`endogenr_spec`) and a `fit_model.<type>_spec()`
  method. This decouples specification from fitting, enables
  `forecast_coefficients()` to refit specs on arbitrary window subsets,
  and makes it straightforward to add new families.

---

## Performance

* **`simulate_system()` now returns forecast-only rows** (time ≥
  `test_start`). Previously the result included the full training period
  replicated across every draw; all in-repo consumers (`get_accuracy()`,
  `sim_to_dist()`, `plotsim()`) already filtered to the forecast window.
  The per-worker trim shrinks peak `rbindlist` allocation proportionally
  to `horizon / n_times_total` (≈ 2–4× for typical horizons).

* **Simulation grid is pruned to referenced columns.** `setup_system()`
  collects every variable named in any model formula or variance formula
  and keeps only those columns (plus group/time/sim keys) in the
  simulation data. Columns present in the input panel but not referenced
  by any model are dropped.

* **Stored lm/glm/heterolm fits shed fitting-only payloads.**
  `.strip_fit_data()` removes `$model`, `$effects`, and
  `$fitted.values` from the inner `lm` object; `$linear.predictors` and
  `$y` from `glm`; and `$frame$X`, `$frame$Z`, `$frame$y` from
  `hetero_fit`. Fields needed by `predict.*(..., se.fit = TRUE)` and
  `predict.hetero_fit(newdata = ...)` are preserved.

* **`.apply_ts_map()` gains a `copy = FALSE` path** for the predict
  route. `predict.linear()`, `predict.glm_endogenr()`, and
  `predict.heterolm()` all call `.history_subset()` before
  `.apply_ts_map()`, which already returns a fresh subset; the second
  deep-copy is now skipped on those paths. The `panel_materialize()`
  fit-time path continues to copy by default.

---

## Documentation

* `?endogenr` carries a "Known issues" / statistical-roadmap section
  documenting unclear ensemble semantics, arbitrary constants, available
  generalisations (unified predictive-draw interface, panel-heterogeneous
  estimation), and brute-forced pieces that warrant proper statistical
  treatment (TVP/state-space coefficient paths, joint SUR/spatial
  innovation draws, distribution-parameter uncertainty, and the likely
  over-dispersion of `linear` + `boot` + `min_window`).

---

## Tests

* Added a layered test harness: regression tests for the three fixes
  above (`test-history-window.R`, `test-univariate-fable.R`,
  `test-glm-model.R`) and four behavioural pillars (`test-dataflow.R`,
  `test-estimation.R`, `test-uncertainty.R`, `test-stepwise.R`) covering
  scheduling, estimation, uncertainty calibration, and step-wise dynamic
  forecasting. Deterministic unit tests are always on; replication-heavy
  statistical tests are gated behind `ENDOGENR_SLOW_TESTS=true` (and
  `skip_on_cran()`) via `helper-skip.R`, with known-truth DGPs in
  `helper-dgp.R` and scoring helpers in `helper-stats.R`. The branch
  ships 21 test files vs 5 on `main`.
