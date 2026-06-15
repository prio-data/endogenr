# HANDOFF: Scope-A full panel design layer unification

This document captures the deferred scope-A unification plan for the panel
design layer introduced in the commit that fixed `poly(dem, 2)` and aligned the
naming contract (scope B). Read in conjunction with `R/panel_transform.R` and
the companion NEWS.md entry.

---

## What scope B delivered

- `R/panel_transform.R` is the canonical panel design layer with:
  - `.pt_ts_fns` / `.pt_roll_fns` / `.pt_cum_fns` — single registry
  - `panel_materialize()` — Stage-1 per-unit ts extraction
  - `.pt_make_aliases()` — `.pt#` → human-readable alias map
  - `.pt_alias_formula()` — formula rewriting for alias substitution
  - `.rewrite_from_map()` — rewrite a sub-formula against an existing map
- `linearmodel` / `glmmodel` / `heterolmmodel` — two-stage fit + predict
- `systemgraph.R` imports `.pt_cum_fns` / `.pt_roll_fns` (registry unified);
  edge-extraction and `.required_history()` logic unchanged
- `longhorizon.R` already used `panel_materialize()`; no changes needed

## What scope A still needs to do

### 1. Derive dependency edges and `.required_history()` from the layer's AST walk

Currently `systemgraph.R` has its own `.classify_term_vars()` / `.edges_from_formula()`
and a separate `.req_hist_expr()` / `.required_history()` depth walker. These
are correct but duplicate intent with the `panel_materialize()` AST walk.

**Full unification would:**
- Extend `panel_materialize()` (or a companion function) to output, in ONE
  pass over the formula AST:
  1. The ts-extraction map (`ts_map`) — already done
  2. Per-variable lag classification (`lagged` vs `plain`) feeding dependency edges
  3. Required history depth for each ts sub-expression
- Retire `.classify_term_vars()` / `.edges_from_formula()` in `systemgraph.R`
  (replace with a thin wrapper over the layer's output)
- Retire the standalone `.req_hist_expr()` / `.required_history()` walker
  (depth now comes from the same map used for materialization)

**Risk:** `.classify_term_vars()` operates on raw term labels from
`stats::terms()` (one label per term, including interactions like `lag(x):z`),
while the layer's rewriter operates on the full expression AST. The two views
are not identical; care is needed around terms that contain BOTH lagged and
plain variables. The existing tests in `test-systemgraph.R` and
`test-history-window.R` must stay green throughout.

### 2. Collapse the three naming surfaces to one documented model

The current design has three naming surfaces:
| Surface | Example | Where |
|---|---|---|
| Synthetic ts columns (internal) | `.pt1` | `panel_materialize()` map |
| Fit column / coef name | `lag_x` | `lm()` coefficient |
| Dependency-graph vertex | `lag_x` (same-period `x`; lagged `lag_x`) | `igraph::V()$name` |

The readable alias (`lag_x`) currently bridges the gap between `.pt1` and the
dependency vertex. Full unification would:
- Document this three-layer model explicitly in `R/panel_transform.R`
- Guarantee that the alias produced by `.pt_make_aliases()` matches exactly
  what the systemgraph uses as vertex names (today they happen to agree for
  simple expressions like `lag(x)` but this is not enforced)
- Potentially add a `panel_design()` S3 object (the planned but deferred object
  from the original plan) that bundles `ts_map`, `alias_map`, `edges`,
  `required_history`, and `fit_formula` as a single inspectable artefact

### 3. Deprecate / remove remaining per-group public utilities

The following exported helpers date from before the two-stage approach:

| Function | Issue |
|---|---|
| `create_panel_frame()` | Internally calls `.build_mat_cache()` which uses `head(.SD, 2L)` per group — same per-group bug for any data-dependent formula term. No longer called internally. |
| `materialize_formula()` | Evaluates the full formula per group via `model.frame()`. Same issue. |
| `derive_naive_formula()` | Constructs a lm formula from janitor-cleaned column names (old naming scheme). |
| `inject_positional_lag()` | Low-level helper; still useful standalone but the layer's `.apply_ts_map()` injects the positional lag internally |

**Recommended migration path:**
1. Add `Deprecated:` section to each Rd page pointing to `panel_materialize()`
2. In a subsequent breaking release, remove them

**Caution:** These are exported; removing them is a user-visible API change that
needs a deprecation cycle.

### 4. Test migration notes

Tests to update when scope A lands:

- `test-systemgraph.R` — if `.edges_from_formula()` is retired, the tests that
  call it directly need to call the new wrapper (or the layer's output)
- `test-history-window.R` — if `.required_history()` is retired, tests calling
  it need updating; the contract (same numeric results) must hold
- `test-panel_transform.R` — add tests for the new unified object (edges, depths
  output in a single pass)
- `test-formula-interactions.R` — currently tests both `derive_naive_formula()`
  (public) and `linearmodel`/`glmmodel` coef names; after deprecation, the
  `derive_naive_formula()` tests can be moved to a deprecation-testing module

### 5. Ordering and risks

Recommended order:
1. Add `panel_design()` S3 object (inspectable, round-trippable) as internal
2. Port `.required_history()` to derive from `panel_design$ts_map` using
   `.req_hist_expr()` on the map values — lowest blast radius
3. Port `.edges_from_formula()` to derive from the panel_design walk — higher
   blast radius (systemgraph integration tests)
4. Deprecate public helpers with a one-version warning
5. Remove in next major version

**Biggest risk:** The `lag_`-vertex naming in the dependency graph and the
`lag_x` alias in the fit formula happen to coincide for simple expressions but
diverge for complex ones. For example, `lag(log(gdppc))` becomes the alias
`lag_log_gdppc` (fit column) but the dependency vertex is `lag_gdppc` (only the
outermost lagged variable is registered). Scope A must decide: do we unify the
vertex name to `lag_log_gdppc` (alias) or keep `lag_gdppc` (just the variable)?
The former is more consistent but changes the graph vertex namespace.
