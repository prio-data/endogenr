# endogenr 0.1.0.9000

## Bug fixes

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
