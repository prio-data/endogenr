# endogenr 0.1.0.9000

## New features

* **`run_experiments()` gains a `windows` axis via the new `window_config()`
  constructor.** The fit-side window approach (`"random"`, `"rolling"`,
  `"expanding"` plus `width` and `step`) and the sim-side window policy
  (`"latest"`, `"equal"`, `"decay"`, or custom `weights`) are now first-class
  experiment variations. Pass a single `window_config()` to apply it uniformly,
  or a named list to cross window approach as a fourth grid axis alongside
  `models`, `train_start`, and `test_start`. Results carry a `window_cfg`
  column and the `experiments` attribute now includes `window_cfg`, `window`,
  and `window_policy` columns. The default (`windows = NULL`) reproduces
  historical behaviour exactly (random windows, latest policy).

## Bug fixes

* **`poly(dem, 2)` and other data-dependent design bases now work inside panel
  model formulas.** Previously, `linearmodel()`, `glmmodel()`, and
  `heterolmmodel()` evaluated the full formula per group via `model.frame()`,
  so `poly(dem, 2)` (and any basis that needs ≥ degree+1 unique values) would
  error with `'degree' must be less than number of unique points` whenever a
  unit had fewer unique values than the degree — a guaranteed failure when
  `.build_mat_cache()` probed with only 2 rows per group. The root cause:
  data-dependent bases must be evaluated **pooled** (matching what a base
  `lm()` on the same data produces), not per unit. The fix is a two-stage
  approach already used by `longhorizon`: Stage 1 materialises only the true
  within-unit time-series sub-expressions (`lag`, `diff`, `rollmean`, …) per
  group; Stage 2 fits a pooled `lm`/`glm` on the rewritten formula, so
  `poly()`, `bs()`, `factor()`, interactions and `-1` are all handled by base
  R's formula machinery and reproduced coherently at predict time via
  `predict.lm()`'s stored `predvars`/`xlevels`.

* **Coherent predict-time basis for heterolm.** `predict.heterolm()` now
  reconstructs `model.matrix()` columns using the training-time `terms` object
  (with `predvars` and `xlevels`), ensuring polynomial/spline bases and factor
  contrasts match those from fitting exactly.

## Architecture

* `R/panel_transform.R` is now the **canonical panel design layer**. The
  time-series function registries (`lag`, rolling, cumulative, decay functions)
  are defined once as `.pt_ts_fns`, `.pt_roll_fns`, `.pt_cum_fns` and shared
  by both materialisation (`panel_materialize()`) and history-depth/edge
  extraction (`systemgraph.R`). The previous `systemgraph.R`-local
  `.rh_cum_fns` / `.rh_roll_fns` duplicates are removed.

* **Naming contract.** Time-series synthetic columns (`.pt#` internal keys) are
  now renamed to human-readable aliases — `janitor`-cleaned deparsed source
  expressions — before the pooled fit, so coefficient names from
  `linearmodel()` / `glmmodel()` match what base `lm()` would produce:
  `lag(x)` → `lag_x`, `lag(log(gdppc))` → `lag_log_gdppc`. The internal
  `.pt#` symbols never appear in `get_coefficients()` output or `$coefs`.

* Formulas now preserve **interactions and full term structure** in model fits.
  Previously, `derive_naive_formula()` rebuilt the formula additively from
  materialized column names, silently dropping interaction terms (`g:x`, `g*x`).
  The fitted `linear` and `glm` models now estimate interactions identically to
  base `lm()`/`glm()`. For `heterolm`, interaction terms are pre-expanded via
  `model.matrix()` into flat predictor columns before fitting (required because
  `heterolm::hetero()` uses column-name lookups rather than R formula algebra).

* `factor(<panel-unit-id>)` in model formulas (e.g. unit fixed effects via
  `y ~ factor(gwcode) + x`) now materializes correctly. Previously, the panel
  key column was absent from the per-group `.SD` inside `model.frame()`, causing
  `"object 'gwcode' not found"`. The new internal `.model_frame_by_group()`
  helper recycles the scalar group-key value into each per-group frame when the
  formula references it. Stacking single-level per-group factors unions to all
  levels across the full dataset, so factor coding and prediction are correct.

* `spatial_lag_model()` now correctly handles island units (geographic units with
  no neighbours). The constructor normalises both the manual `integer(0)`
  representation and the `sfdep`/`spdep` `0L` sentinel to `0L`, and detects
  islands via a predicate robust to both encodings. Previously, `integer(0)`
  islands caused `sfdep::st_lag()` to error with `"zero length neighbour vector"`,
  and `0L` islands silently returned `0` instead of `island_default` because
  `lengths(nb) == 0` evaluated to `FALSE` for the length-1 sentinel. Both paths
  now correctly return `island_default` (or `NA_real_` when unset).

* `predict.linear()`, `predict.glm_endogenr()`, and `predict.heterolm()` now
  size the per-unit history window by a new internal `.required_history()` AST
  walk that **composes** nested time-series depths, instead of taking the
  maximum of individual depths (`.max_lag_depth()`). Composed transforms such as
  `lag(lag(x))`, `lag(rollmeanr(x, k), n)`, and any cumulative transform
  (`cumsum()`, `decay_since_event()`, ...) are no longer silently truncated to
  too short a window, which previously produced `NA` or wrong values in forecast
  cells. Cumulative/since-event transforms now correctly use the full per-unit
  history. The same sizing replaces the hard-coded 20-row history in the
  long-horizon `.predict_lh()`.

* `predict.univariate_fable()` now draws **coherent sample paths** with
  `fabletools::generate()` on the fitted mable (innovations simulated forward
  with the model's temporal structure), rather than stitching together
  independent per-horizon draws from the marginal `forecast()` distributions.
  Each `(unit, sim)` trajectory is now a single internally-consistent path, so
  the autocorrelation an ETS/ARIMA implies survives into the endogenous
  equations that read it via `lag()`.

* `getpi_glm()` now adds **response-scale dispersion** in addition to the
  link-scale parameter uncertainty. A Gaussian GLM's prediction interval again
  matches the equivalent `lm` interval; Poisson/quasipoisson draws are
  non-negative integers with the right mean/variance; Gamma draws are positive
  with variance `dispersion * mu^2`; binomial/quasibinomial proportion outcomes
  are drawn from a Beta with mean `mu` (see `?endogenr` for the proportion
  assumption). Unsupported families fall back to the previous link-scale-only
  draw with a one-time warning.

## Performance

* **`simulate_system()` now returns forecast-only rows** (time ≥ `test_start`).
  Previously the result included the full training period replicated across
  every draw; all in-repo consumers (`get_accuracy()`, `sim_to_dist()`,
  `plotsim()`) already filtered to the forecast window, so their behaviour is
  unchanged. The per-worker trim shrinks peak `rbindlist` allocation and
  `object.size(res)` proportionally to `horizon / n_times_total`
  (≈ 2–4× for typical horizons).

* **Simulation grid is pruned to referenced columns.** `setup_system()` now
  collects every variable named in any model formula or variance formula and
  keeps only those columns (plus group/time/sim keys) in the simulation data.
  Columns present in the input panel but not referenced by any model are
  dropped, shrinking the grid, each per-draw copy, and the final result.

* **Stored lm/glm/heterolm fits shed fitting-only payloads.** `.strip_fit_data()`
  now removes `$model`, `$effects`, and `$fitted.values` from the inner `lm`
  object; `$linear.predictors` and `$y` from `glm`; and `$frame$X`, `$frame$Z`,
  `$frame$y` from `hetero_fit` (the training design matrices). The fields needed
  by `predict.*(..., se.fit = TRUE)` and `predict.hetero_fit(newdata = ...)` are
  preserved; parity is verified by new tests.

* **`.apply_ts_map()` gains a `copy = FALSE` path** for the predict route.
  `predict.linear()`, `predict.glm_endogenr()`, and `predict.heterolm()` all
  call `.history_subset()` before `.apply_ts_map()`, which already returns a
  fresh subset; the second deep-copy is now skipped on those paths. The
  `panel_materialize()` fit-time path continues to copy by default.

## Documentation

* `?endogenr` now carries a "Known issues" / statistical-roadmap section
  documenting the unclear ensemble semantics, the arbitrary constants, the
  available generalisations (unified predictive-draw interface, panel-
  heterogeneous estimation), and the brute-forced pieces that warrant a proper
  statistical treatment (TVP/state-space coefficient paths, joint SUR/spatial
  innovation draws, distribution-parameter uncertainty, and the likely
  over-dispersion of `linear` + `boot` + `min_window`).

## Tests

* Added a layered test harness: regression tests for the three fixes above
  (`test-history-window.R`, `test-univariate-fable.R`, `test-glm-model.R`) and
  four behavioural pillars (`test-dataflow.R`, `test-estimation.R`,
  `test-uncertainty.R`, `test-stepwise.R`) covering scheduling, estimation,
  uncertainty calibration, and step-wise dynamic forecasting. Deterministic
  unit tests are always on; replication-heavy statistical tests are gated behind
  `ENDOGENR_SLOW_TESTS=true` (and `skip_on_cran()`) via `helper-skip.R`, with
  known-truth DGPs in `helper-dgp.R` and scoring helpers in `helper-stats.R`.
