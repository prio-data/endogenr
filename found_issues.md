# Issues found but left out of scope — resolved items noted inline

Discovered while building the test harness and fixing the three approved defects
(history-window under-read, fable path-incoherence, GLM dispersion). These were
**not** part of that scope, so they are recorded here rather than fixed. The
broad statistical roadmap (arbitrary constants, SUR/TVP, parameter-uncertainty
propagation, unified draw interface, pooled-vs-heterogeneous estimation) already
lives in `?endogenr` ("Known issues"); this file lists the concrete code-level
problems found during the work that are *not* in that roadmap.

---

## ~~1. `predict.spatial_lag` island handling is broken~~ — **FIXED**

**Where:** `R/spatiallagmodel.R` (`spatial_lag_model()` constructor).

**Root cause:** `lengths(nb) == 0` failed to detect real `sfdep` islands (encoded
as `0L`, length 1) and `integer(0)` islands caused `sfdep::st_lag()` to error
before the `island_default` override was reached.

**Fix applied:** `spatial_lag_model()` now normalises both representations (`integer(0)`
and `0L`) to the spdep `0L` sentinel at construction, and detects islands via
`length(x) == 0L || (length(x) == 1L && isTRUE(x[[1]] == 0))`. The stored
`nb`/`wt` are always in the normalised form; `predict.spatial_lag()` is unchanged.

**Tests added:** `test-stepwise.R` — three new `test_that` blocks covering the
`integer(0)` error case, the silent-`0` `0L` case, and mask-correctness at
construction for both encodings and the all-connected no-op regression.

---

## ~~2. `get_train_window()` can return `end = test_start` (doc/impl mismatch)~~ — **FIXED**

**Where:** `R/utilities.R` lines 25–40 (docstring lines 5–6, 12 vs body lines 37–39).

**What:** The docstring promises the window "always ends at `test_start - 1` or
earlier". The implementation is
`end = test_start - stop_decrement` with
`stop_decrement = sample.int(...) - 1`, so a `0` decrement yields
`end = test_start` — i.e. the forecast origin period is nominally inside the
returned window.

**Why it is harmless (today):** fitting always runs on `sys$train_data`, which is
pre-filtered to `time < test_start`, so even `end = test_start` cannot pull a
forecast-period row into a fit. The bug is a **misleading label**, not actual
leakage.

**Fix applied:** the body now caps the random window at the forecast origin
(`end = min(test_start - 1L, test_start - stop_decrement)`) and returns the
full window `[earliest_train_start, test_start - 1]` when
`min_window == full_range` (previously `sample.int(0)` errored). The
`?endogenr` known-issues bullet and the bound expectations in
`tests/testthat/test-estimation.R` were updated; random-window RNG streams
shift accordingly (the branch already declares results non-bit-identical).

---

## ~~3. `bootstrapglm()` on count families: spurious warnings~~ — **warnings FIXED**; SE inflation documented

**Where:** `R/glmmodel.R` (`bootstrapglm()`), exercised by `test-estimation.R`.

**What:** For a Poisson/quasipoisson GLM, the residual bootstrap reconstructs `y`
on the link scale and refits with `glm(..., family = poisson())` on the
**continuous** reconstructed response. This emits a stream of
`non-integer x = ...` warnings, and the bootstrap coefficient SE materially
**over-disperses** relative to the analytic Fisher SE.

**Evidence (this session):**
```
GLM (poisson) x: 0.475 | boot mean: 0.466 | boot sd: 0.1134 | analytic SE: 0.0437
```
The point estimate is recovered (boot mean ≈ MLE), but the spread is ~2.6× the
analytic SE — the working-residual bootstrap on the link scale is not a
calibrated SE estimator for non-Gaussian families.

**Why out of scope:** this is the same family of mis-calibration the roadmap flags
for `linear + boot + min_window` (over-counting parameter uncertainty); it is
documented, not fixed. The non-integer warning is a cosmetic code smell on top.

**Fix applied (warnings):** `bootstrapglm()` now muffles only the
`non-integer` warning around its refit via `withCallingHandlers()` — the
continuous link-scale reconstruction is intentional; other warnings still
surface. The SE-calibration concern is statistical, remains out of scope, and
stays documented in `?endogenr` ("Known issues"): prefer a parametric/cases
bootstrap for GLMs, or treat `boot` SEs for non-Gaussian families as
indicative only.

---

## ~~4. Latent time-column / argument name collision in the predict subset~~ — **FIXED (comprehensively)**

**Where:** originally `R/systemgraph.R` `.history_subset()`; in fact every
`data[data[[idx]] == t]`-style subset across the predict methods, the engine,
validation, and scoring.

**What:** subsets evaluated `data[[idx]] == t` inside the data.table `i`
scope. If the time column is literally named `t` (or shares a name with any
local/argument used in the expression), data.table column scoping shadows the
argument and the filter silently degenerates (e.g. `t == t` → all rows). The
degenerate per-step filter then made the update-join overwrite observed
training rows with `NA`, cascading `NA` through the whole forecast.

**Status:** the earlier claim that fixing `.history_subset()` closed this was
wrong — the per-step time filter in every `predict.*` method (and a dozen
other engine/validation subsets) carried the same collision. All row subsets
now route through `.dt_rows()` (`R/utilities.R`), which applies a mask
computed outside any data.table scope via the `env=` substitution mechanism.
Regression tests: `test-systemsim.R` ("full pipeline works when the time
column is literally named 't'") and the `flag_war` end-to-end test in
`test-systemgraph.R`.

---

## 5. `R CMD check` baseline (updated after panel-context refactor)

`devtools::check()` run against `HEAD` (`refactor/panel-context`) reports:
**0 errors, 0 warnings, 3 notes.** The following are pre-existing or deferred:

- ~~**WARNING — invalid license file pointer.**~~ **FIXED** — `LICENSE` stub
  (two-line MIT year/holder) created alongside `LICENSE.md`; `DESCRIPTION`
  keeps `MIT + file LICENSE`.
- ~~**WARNING — undeclared `actuar`.**~~ **FIXED** — `actuar` added to
  `Suggests` in `DESCRIPTION`.
- ~~**NOTE — non-standard top-level files** (`data_raw/`, `examples/`,
  `lm_prediction_interval_test.R`, `package_test.R`, `HANDOFF.md`,
  `found_issues.md`, `.understand-anything/`).~~ **FIXED** — all added to
  `.Rbuildignore`.
- **NOTE — S3 generic/method consistency.** Every `predict.*` method uses
  `model` as its first argument instead of `object` (the `predict` generic's
  first formal). Fix: rename the first arg to `object` (and update call sites).
  Deferred — internal-only dispatch, tests pass, no user-visible impact.
- **NOTE — "no visible binding for global variable".** The usual `data.table` /
  `rlang` / `stats` NSE symbols (`:=`, `.SD`, `.data`, `predict`, `na.omit`,
  `residuals`, `terms`, the quantile column names, etc.). Fix: add
  `utils::globalVariables(...)` and the `importFrom("stats", ...)` the check
  suggests. Deferred — cosmetic; does not affect runtime.
- **NOTE — `R (>= 4.1.0)` dependency.** `R CMD build` auto-added because source
  files use `|>` and `\(...)` shorthand. Deferred — simply bump `Depends:
  R (>= 4.1.0)` in `DESCRIPTION` when targeting older R versions becomes a
  requirement.

---

## 6. `examples/growth_forecasting.R` requires Remotes packages not on CRAN

**Where:** `examples/growth_forecasting.R` — uses `vdemdata` and `poldat`,
both in `DESCRIPTION Remotes`, not installed in a fresh environment.

**What:** With those packages missing, the data filtering step produces 0 units,
`setup_system()` prints an empty panel, and the subsequent `ggplot()` call fails
with `object 'p.value' not found` (empty dataset → empty coefficients table with
correct schema but 0 rows — ggplot2 treats missing columns as absent). This is
not a code regression: the README pipeline, which uses the bundled `example_data`,
runs end-to-end without error.

**Suggested fix:** Either gate the vdemdata sections with
`if (requireNamespace("vdemdata", quietly = TRUE))` or add a note at the top of
the file explaining the data requirements. Not a merge blocker.

---

# Issues from the 2026-07 code audit — noted, not changed

Found during the full-codebase audit that fixed the execution-order `lag_`
substring match, the data.table scope collisions (item 4 above), the GLM
shared-parameter-draw, the `diff()`/`shift()`/`lead()` materialization gaps,
and the fit-closure serialization bloat. The items below were deliberately
left alone; each says why. The audit also re-confirmed the
**parameter-uncertainty asymmetry** across families (`linear`/`glm` add a
per-row se.fit draw; `heterolm`/`gamlss` sample at point-estimate
parameters) — that stays tracked in `?endogenr` ("Known issues"), not here.

---

## 7. `simulate_system()` copies the full training history once per draw

**Where:** `R/systemsim.R`, both `future_lapply` branches of
`simulate_system()` (`sim <- data.table::copy(simulation_data)` per outer
draw).

**What:** every outer draw deep-copies the entire simulation grid, including
all training rows back to `train_start`, but the dynamic loop only ever reads
the composed history depth (`.required_history()`) behind each forecast step.
With a 40-year training window and lag-1 models, ~85% of the copied rows are
never read; copy time and per-worker memory scale with `nsim x grid`.

**Suggested fix:** at setup, compute the max **finite** history need across
specs (the `.check_entry_depth()` walk already does this). When every spec's
need is finite, trim the per-draw copy to
`[test_start - max_need, test_start + horizon - 1]`; any spec with a
cumulative/since-event term (`Inf` need) keeps the full history — the same
fallback `.history_subset()` uses.

**Why not changed:** it alters the simulation-grid contract (models see a
shorter frame), so it needs its own regression test proving trimmed vs full
output is bit-identical for finite-need systems before landing. Pure perf,
no correctness impact today.

---

## 8. `intensity_decay()` is O(n²) per group — and re-evaluated per step

**Where:** `R/formula_helpers.R` (`intensity_decay()`).

**What:** the nested loop re-accumulates every past event for each row i,
O(n²) per (unit, sim) series. Because it is registered as a cumulative ts
function (`Inf` history), the predict path re-evaluates it over the FULL
per-unit history at **every** forecast step: total cost
O(horizon x sims x units x n²) for formulas that use it.

**Why not changed:** the per-event floor
`max(exp(-lambda*dt), exp(-lambda*max_years))` breaks the plain geometric
recurrence `total_i = total_{i-1} * exp(-lambda) + I_i`. A correct O(n)
version needs a two-part decomposition — a decaying sum over events younger
than `max_years` plus `floor x cumsum(intensity)` over older ones. That is
correctness-sensitive in a statistical helper and deserves its own targeted
tests (equality against the brute-force version over random event/intensity
series, boundary at `dt == max_years`).

---

## 9. User columns named `lag_*` share the synthetic dependency-graph namespace

**Where:** `R/systemgraph.R` (`.edges_from_formula()` / `parse_formula()`:
`paste0("lag_", var)` vertices).

**What:** the substring bug (a variable named `flag_war` being stripped as a
lag vertex) is fixed — the strip is now anchored to the `lag_` **prefix**.
But a data column literally named `lag_x`, used as an outcome or predictor,
still lands on the same vertex the graph synthesises for `lag(x)`: the two
are conflated silently (e.g. `y ~ lag(x)` plus an `exogen ~lag_x` share one
vertex, distorting ordering/closure checks).

**Suggested fix:** either move the synthetic namespace out of the legal
column-name space (e.g. a `.lag:` prefix — column names can contain it, but
never by accident), or make `validate_system_closure()` error when a data
column starts with `lag_` **and** its suffix also appears inside a `lag()`
call in some formula.

**Why not changed:** low collision probability, and renaming the vertex
namespace touches ordering, closure validation, and their tests; not worth
coupling to the audit fixes. Silent when hit, though — hence recorded.

---

## 10. Slow statistical tests never run in a default suite

**Where:** `tests/testthat/helper-skip.R` (`skip_if_not_slow()`); 20 tests
across `test-estimation.R`, `test-uncertainty.R`, `test-stepwise.R`,
`test-univariate-fable.R`, `test-glm-model.R`, `test-glmmTMB.R`,
`test-gamlss.R`.

**What:** the calibration-coverage, PI-widening, coefficient-recovery,
ar1-decay, and moment-matching claims are exercised only with
`ENDOGENR_SLOW_TESTS=1`. A default `devtools::test()` (and any CI that does
not set the variable) never checks them, so a statistical regression can land
silently while the fast suite stays green.

**Suggested fix:** add a scheduled (cron) CI job running the suite with
`ENDOGENR_SLOW_TESTS=1`; keep PR runs fast. The full slow suite was verified
green during this audit (0 failed / 2285 passed), so enabling it costs only
~30s of runtime today.