# endogenr Refactoring Hand-off

## Current State (2026-05-21)

**Branch**: `refactor/panel-context`  
**Tests**: 135 passing, 0 failures  
**`devtools::check()`**: 0 errors, warnings (pre-existing missing Rd descriptions), notes (cosmetic)

### What has been done

**Phases 0–3** (all committed) of a 7-phase refactoring plan. The full plan is in `~/.claude/plans/glistening-hopping-llama.md`.

#### Phase 0: Bug Fixes (committed as `7a45fad`)

- **Parse formula fix** (`R/systemgraph.R`): Rewrote the dependent-model branch of `parse_formula()` to use `.edges_from_formula()` which iterates per formula term. Previously, `ifelse()` recycling misclassified variables that appeared both lagged and unlagged (e.g., `y ~ lag(x) + x` would only produce `lag_x -> y`, missing `x -> y`).

- **Max history fix** (`R/systemgraph.R`): Added `.max_lag_depth(formula)` which walks the formula AST to extract lag depths from `lag(expr, n)` calls and `k` from rolling functions. Replaces the old regex heuristic (`str_extract_all(formula, "[0-9]+")`) which broke on `I(100)`, variable names with digits, etc. Model constructors now compute `model$max_history` at construction time.

#### Phase 1: panel_context + data.table Migration (committed as `ecaa130`)

**Core change**: Replaced tsibble metadata (key_vars/index_var) and dplyr operations with `panel_context` + plain data.table throughout the package.

Files changed (all uncommitted):

| File | What changed |
|------|-------------|
| `R/panel_context.R` | **NEW**. `panel_context()`, `ctx_unit()`, `ctx_time()`, `ctx_sim()`, `ctx_keys()`. Also contains `.datatable.aware = TRUE`. |
| `R/utilities.R` | `inject_positional_lag()` (new). `create_panel_frame()` rewritten to accept `ctx`, uses `split()`/`lapply()` with `model.frame()` per group. |
| `R/systemsim.R` | Full rewrite. `prepare_simulation_data()` uses `CJ()`. `process_independent/dependent_models()` use data.table update-joins. `simulate_endogenr()` returns `data.table`. `sim_to_dist()` uses list-columns (no distributional). `get_accuracy()` uses `scoringRules::crps_sample()` (no fabletools). `plotsim()` uses direct ggplot2. |
| `R/linearmodel.R` | Accepts `ctx`, uses `create_panel_frame(formula, data, ctx)`, data.table subsetting. |
| `R/glmmodel.R` | Same pattern as linear. |
| `R/heterolmmodel.R` | Same pattern, dual formula (mean + variance) handling preserved. |
| `R/deterministicmodel.R` | `predict.deterministic()` uses `inject_positional_lag()` + data.table. |
| `R/exogenmodel.R` | Uses `CJ` grid expansion, returns data.table. |
| `R/univariate_fable_model.R` | Accepts `ctx`, converts to tsibble internally for fable fitting, returns data.table. fable/fabletools/tsibble moved to Suggests. |
| `R/parametric_distribution_model.R` | Replaced `distributional::generate()` with `.sample_from_fitdist()` that calls base R `r*` functions directly. `create_distribution_object()` removed. |
| `R/spatiallagmodel.R` | Accepts `ctx`, returns data.table instead of tsibble. |
| `R/longhorizon.R` | Migrated to data.table. `get_lh_accuracy()` uses `scoringRules::crps_sample()` instead of `fabletools::accuracy()`. |
| `R/basemodel.R` | Updated example (removed tsibble assumption). `build_model()` unchanged structurally — still uses `purrr::partial()`. |
| `DESCRIPTION` | Added `scoringRules`, `rlang`, `broom`, `ggplot2`, `igraph`, `janitor`, `purrr`, `future` to Imports. Moved `fable`, `fabletools`, `tsibble`, `distributional` to Suggests. Removed `R6`, `pbapply`. |
| `NAMESPACE` | Added exports for `panel_context`, `ctx_*`, `inject_positional_lag`. Removed `create_distribution_object`. |
| `tests/testthat/test-systemsim.R` | Rewritten: `make_test_data()` returns data.table, assertions check for `data.table` class, added `get_accuracy` test, tsibble input acceptance test. |

#### Phase 2: Validations (uncommitted)

| File | What changed |
|------|-------------|
| `R/validate.R` | **NEW**. `validate_panel()` checks integer time steps, contiguous series, balanced panel, complete initial state, sorted data. `validate_system_closure()` checks no duplicate outcomes, every RHS variable available from model outcomes or data columns. |
| `R/systemsim.R` | `setup_simulator()` calls `validate_panel()` and `validate_system_closure()` before fitting. `prepare_simulation_data()` warns when estimated grid > 1 GB. `process_independent_models()` and `process_dependent_models()` wrap predict calls in `tryCatch` with model name + time step in error messages. |
| `R/basemodel.R` | `build_model()` stores `formula` and `variance` as attributes on the partial function, so validations can access them before fitting. |
| `NAMESPACE` | Added exports for `validate_panel`, `validate_system_closure`. |
| `tests/testthat/test-validate.R` | **NEW**. 16 tests covering all validation paths + integration tests via `setup_simulator()`. |

#### Phase 3: Model Spec/Fit Separation (uncommitted)

| File | What changed |
|------|-------------|
| `R/basemodel.R` | `build_model()` returns `endogenr_spec` list (not `purrr::partial`). Added `fit_model()` S3 generic. |
| `R/linearmodel.R` | Added `fit_model.linear_spec()`. `linearmodel()` kept as constructor. |
| `R/glmmodel.R` | Added `fit_model.glm_spec()`. `glmmodel()` kept. |
| `R/heterolmmodel.R` | Added `fit_model.heterolm_spec()`. `heterolmmodel()` kept. |
| `R/deterministicmodel.R` | Added `fit_model.deterministic_spec()`. |
| `R/exogenmodel.R` | Added `fit_model.exogen_spec()`. |
| `R/parametric_distribution_model.R` | Added `fit_model.parametric_distribution_spec()`. |
| `R/univariate_fable_model.R` | Added `fit_model.univariate_fable_spec()`. |
| `R/spatiallagmodel.R` | Added `fit_model.spatial_lag_spec()`. |
| `R/systemsim.R` | `setup_simulator()` calls `fit_model()` directly (no partials). `inner_simulation()` accepts `specs` + `pre_fitted`, refits linear/glm/heterolm per outer sim. `simulate_endogenr()` simplified — no pre-fit loop. `_subset` type variants eliminated. Setup returns `specs` instead of `models`. |
| `R/validate.R` | `validate_system_closure()` reads `$formula` and `$type` from specs directly. |
| `R/plot_estimates.R` | Replaced `purrr::map_dfr`/`map2_chr` with base `lapply`/`do.call(rbind, ...)`/`mapply`. |
| `DESCRIPTION` | Removed `purrr` from Imports. |
| `NAMESPACE` | Added `fit_model` export + S3method registrations for all spec types. |
| `tests/testthat/test-build_model.R` | Rewritten: tests spec structure instead of function/partial. Added `fit_model` tests. |
| `tests/testthat/test-systemsim.R` | Updated: checks for `specs` instead of `models` in setup result. |

### Gotchas discovered during Phase 3

1. **Roxygen `@export` on S3 methods generates bogus exports.** Using `@export` on `fit_model.linear_spec` caused roxygen to export words from the title as functions (e.g., `export(Fit)`, `export(a)`, `export(linear)`). **Fix**: Use `@exportS3Method` instead of `@export` on S3 method definitions.

2. **Separate fit context needed.** `setup_simulator()` stores two contexts: `ctx` (with `sim = "sim"`, for predict) and `fit_ctx` (without sim, for fitting). `inner_simulation()` uses `fit_ctx` when re-fitting models so the training data (which has no sim column) works correctly.

3. **`_subset` type variants eliminated.** Previously `setup_simulator()` created synthetic types like `"linear_subset"`, `"glm_subset"` etc. Now `subset` is just a parameter passed to `fit_model()`. The `spec$type` stays clean (`"linear"`, `"glm"`, `"heterolm"`).

### Gotchas discovered during Phase 2

1. **`build_model()` returns closures, not lists.** (Now resolved by Phase 3 — specs are plain lists.)

2. **`model$variance_formula` on closures.** (Now resolved — specs have `$args$variance` directly.)

### Gotchas discovered during Phase 1

1. **`.datatable.aware = TRUE` is mandatory.** Without it, data.table's `[` method is not dispatched correctly from package namespace functions. `.SD` becomes empty, `data[data[[col]] >= val]` falls back to `[.data.frame` and fails with "undefined columns selected". This single line in `R/panel_context.R` fixes everything.

2. **Shared ctx mutation bug.** `setup_simulator()` originally did `ctx$sim <- "sim"` to add the sim dimension after fitting. But `ctx` is a list shared by reference with the model partials (via `purrr::partial(x, ctx = ctx)`). When `inner_simulation()` re-fits models, they receive a ctx that now has `sim = "sim"` but the training data has no "sim" column. **Fix**: create a separate `sim_ctx <- panel_context(unit, time, sim = "sim")` instead of mutating the original.

3. **`model.frame()` + data.table `[, by]` scoping.** `model.frame()` uses `match.call()` and `eval()` internally, which can have scoping issues inside data.table's `[, expr, by]` when called from a package namespace. Adding `.datatable.aware = TRUE` fixes most cases. The `create_panel_frame` function now uses `split()`/`lapply()` as a belt-and-suspenders approach, which is simpler and avoids the issue entirely.

4. **NA handling in `get_accuracy()`.** `scoringRules::crps_sample()` does not accept NA values. Simulation draws can contain NAs (e.g., when an exogenous variable is unavailable in the forecast period). The scoring functions now filter NAs before calling crps_sample and use `na.rm = TRUE` in aggregation.

---

## What remains (Phases 4-6)

### Phase 4: Performance (Pre-materialize)

- Split `create_panel_frame()` into `materialize_formula()` (pure data transform) and `derive_naive_formula()` (explicit outcome arg).
- Pre-materialize formulas before the step loop in `process_dependent_models()`. In the step loop, only update lag columns and slice rows — don't re-call `create_panel_frame()`. This eliminates `create_panel_frame()` from the O(H * M) hot loop.
- Cache naive formulas at model construction time.

### Phase 5: Parallelism Fix

- Remove internal `future::plan()` management — let users set their own plan.
- Replace the manual `future::future()` loop with `future.apply::future_lapply()`.
- Add `progressr` support.
- Deprecate `parallel`/`ncores` args with message pointing to `future::plan()`.

### Phase 6: NAMESPACE Cleanup + Contracts

- Mark internal functions with `@keywords internal`, remove `@export`.
- Keep exported: `build_model`, `setup_simulator`, `simulate_endogenr`, `fit_model`, `sim_to_dist`, `get_accuracy`, `plotsim`, `plot_estimates`, `get_execution_order`, `get_train_window`, long-horizon functions, helper functions (`decay_since_event`, etc.).
- Update `CONTRACTS.md` to reflect all changes.
- Run `devtools::document()` to regenerate NAMESPACE from roxygen.

### Phase dependency chain

```
Phase 0 (done) -> Phase 1 (done) -> Phase 2 (done) -> Phase 3 (done) -> Phase 4 -> Phase 5 -> Phase 6
```

Each phase should keep tests green. After Phase 3 (spec/fit separation), the package API is stable and Phases 4-6 are internal optimizations/cleanup.

---

## How to continue

1. **Pick up Phase 4** (pre-materialize performance). Split `create_panel_frame()` into `materialize_formula()` and `derive_naive_formula()`. Pre-materialize before the step loop in `process_dependent_models()`.

2. **Phases 5-6** (parallelism fix + NAMESPACE cleanup) are incremental improvements. Phase 5 replaces the manual future loop with `future.apply::future_lapply()`. Phase 6 marks internal functions with `@keywords internal`.

---

## Key architectural decisions made

| Decision | Rationale |
|----------|-----------|
| data.table as internal carrier | Faster than tibble for large panels, `setkeyv()` for O(log n) joins, no metadata loss on operations |
| `panel_context` as metadata carrier | Decouples metadata from data object. Survives copy, subset, merge without silent loss |
| Output type is `data.table` | Users convert to tibble/tsibble if needed. Simpler than maintaining dual output |
| Keep fable/fabletools in Suggests | `univariate_fable_model` genuinely needs fable. Guarded by `requireNamespace()` |
| scoringRules for accuracy | Direct CRPS/MAE computation without fabletools. Already in CRAN |
| `.sample_from_fitdist()` replaces distributional | Direct `r*()` calls from fitdistrplus results. No distributional dependency for parametric models |
| `split()/lapply()` in `create_panel_frame()` | Avoids `model.frame()` + data.table `[, by]` scoping issues from package namespace. Simpler, equally fast for typical panel sizes |
| Spec/fit separation over `purrr::partial()` | Specs are transparent lists (`$type`, `$formula`, `$args`). `fit_model()` generic enables clean re-fitting per outer sim. No closure opacity, no `attr()` hacks. `purrr` dependency removed. |
