# Panel-aware transform module -------------------------------------------------
#
# Separates two concerns that `model.matrix()` conflates:
#   (1) panel-aware data handling  -- time-series transforms evaluated PER UNIT
#                                     on time-ordered data (lag, decay, etc).
#   (2) model specification        -- design construction evaluated POOLED
#                                     (factor contrasts, poly/spline bases,
#                                     interactions, intercept handling).
#
# Stage 1 (here) materialises the maximal time-series sub-expressions of a
# formula into synthetic `.pt#` columns, per group and in time order, and
# rewrites the formula to reference those columns. Stage 2 (the caller) fits a
# pooled `lm`/`bootstraplm` with the rewritten formula, so `factor`, `poly`,
# splines, interactions and `-1` are all resolved correctly on the pooled data
# and reproduced by `predict.lm` via its `predvars`/`xlevels`.

# Cumulative time-series functions: require the entire prior series (depth Inf).
# Exported as `.pt_cum_fns` so systemgraph.R can import a single canonical list.
.pt_cum_fns <- c(
  "cumsum", "cumprod", "cummax", "cummin",
  "decay_since_event", "time_since_event", "intensity_decay"
)

# Rolling/window time-series functions: require a finite window of history.
# Includes all zoo/data.table roll* variants. Exported as `.pt_roll_fns` so
# systemgraph.R can import the canonical list (replaces .rh_roll_fns there).
.pt_roll_fns <- c(
  "rollmean", "rollmeanr", "rollmeanl",
  "rollsum",  "rollsumr",  "rollsuml",
  "rollmax",  "rollmaxr",  "rollmaxl",
  "rollmedian", "rollmedianr", "rollmedianl",
  "rollapply",  "rollapplyr",  "rollapplyl",
  "frollmean", "frollsum", "frollapply", "frollmax", "frollmin"
)

# Time-series functions: extracted to a per-group `.pt#` column. Matched by
# bare name or `pkg::fn`. The whole call (e.g. the entire `lag(...)`) is
# materialised, so anything wrapped in one is always evaluated within group.
# This is the single canonical registry — consumed by both panel_materialize()
# and systemgraph.R (.required_history / edges).  Callers may extend this set
# per fit via the `ts_fns` argument of panel_materialize().
.pt_ts_fns <- c(
  "lag", "lead", "lead_horizon", "shift",
  .pt_roll_fns,
  .pt_cum_fns,
  "diff"
)

# Design operators: NOT extracted. The call head is kept and its arguments are
# recursed into, so time-series terms nested inside them are still extracted
# (e.g. `poly(lag(x), 2)` -> `poly(.pt1, 2)`). Listed for documentation; the
# rewrite treats them identically to any other non-time-series call.
.pt_design_ops <- c(
  "factor", "ordered", "poly", "ns", "bs", "cut", "C", "interaction"
)

# Positional within-group lag, matching inject_positional_lag(): `lag(x, n)`
# shifts `x` forward by `n` rows, padding the head with NA. Type-preserving
# (index-based) so factor/Date/character columns keep their class and levels.
.pt_positional_lag <- function(x, n = 1L) {
  n   <- as.integer(n)
  len <- length(x)
  if (n >= len) return(x[rep(NA_integer_, len)])
  x[c(rep(NA_integer_, n), seq_len(len - n))]
}

# Name of the function a call invokes, stripping a `pkg::`/`pkg:::` qualifier.
# Returns NA for a call whose head is itself a complex expression.
.pt_call_name <- function(call) {
  head <- call[[1L]]
  if (is.call(head) &&
        (identical(head[[1L]], as.name("::")) ||
           identical(head[[1L]], as.name(":::")))) {
    return(as.character(head[[3L]]))
  }
  if (is.symbol(head)) return(as.character(head))
  NA_character_
}

# Recursively rewrite an expression, extracting maximal time-series
# sub-expressions into `.pt#` symbols. `state` is an environment accumulator
# carrying `$ts_fns` (the active name set), `$map` (named list `{.pt#: ts_expr}`),
# `$keys` (deparse -> symbol, for dedup) and `$counter`. Returns the rewritten
# expression; the extracted map is read from `state$map` by the caller.
.rewrite_panel_formula <- function(expr, state) {
  # Symbols and literals pass through untouched.
  if (!is.call(expr)) return(expr)

  fname <- .pt_call_name(expr)

  # Time-series call: replace the whole expression with a `.pt#` symbol,
  # deduplicating identical sub-expressions by their deparse.
  if (!is.na(fname) && fname %in% state$ts_fns) {
    key <- deparse(expr, width.cutoff = 500L)
    existing <- state$keys[[key]]
    if (!is.null(existing)) return(as.symbol(existing))
    state$counter <- state$counter + 1L
    sym <- paste0(".pt", state$counter)
    state$map[[sym]] <- expr
    state$keys[[key]] <- sym
    return(as.symbol(sym))
  }

  # Any other call (design operator, log/I/arithmetic, `:`/`*`/`^`): keep the
  # head, recurse into the arguments, and rebuild. Argument names are preserved.
  parts <- as.list(expr)
  if (length(parts) > 1L) {
    for (i in 2:length(parts)) {
      parts[[i]] <- .rewrite_panel_formula(parts[[i]], state)
    }
  }
  as.call(parts)
}

# Materialise a `.pt#` map into per-group, time-ordered columns. Returns `data`
# (keyed by group then time) with one new column per map entry. When `copy =
# TRUE` (default) a deep copy is made first so the caller's data is not
# modified. Pass `copy = FALSE` when the caller already owns a private copy
# (e.g. `.history_subset()` output on the predict path) to avoid the redundant
# allocation.
.apply_ts_map <- function(map, data, groupvar, timevar, env = parent.frame(), copy = TRUE) {
  if (!data.table::is.data.table(data)) {
    data <- data.table::as.data.table(as.data.frame(data))
    copy <- FALSE  # as.data.table already produced a new object
  }
  if (copy) data <- data.table::copy(data)
  data.table::setkeyv(data, c(groupvar, timevar))
  if (length(map) == 0L) return(data)

  ts_env <- new.env(parent = env)
  ts_env$lag <- .pt_positional_lag

  for (sym in names(map)) {
    expr <- map[[sym]]
    data[, (sym) := eval(expr, .SD, ts_env), by = groupvar]
  }
  data
}

# Rewrite a formula into pooled form and materialise its time-series terms.
# Both sides are rewritten through the shared map, so a time-series expression
# appearing on the LHS (rare, e.g. `log(diff(.target))`) is handled too.
#
# `ts_fns` is an optional character vector of extra function names to treat as
# time-series (unioned with the built-in `.pt_ts_fns`), letting an advanced
# caller register a custom within-group function used unwrapped in the formula.
#
# Returns a list with:
#   * `data`    -- a copy of `data` with the `.pt#` columns added,
#   * `formula` -- the rewritten formula (same environment as `formula`),
#   * `map`     -- the named `{.pt#: ts_expr}` list.
#' @keywords internal
panel_materialize <- function(formula, data, groupvar, timevar,
                              ts_fns = NULL) {
  env <- rlang::f_env(formula)

  state <- new.env(parent = emptyenv())
  state$ts_fns  <- union(.pt_ts_fns, ts_fns)
  state$map     <- list()
  state$keys    <- list()
  state$counter <- 0L

  lhs <- rlang::f_lhs(formula)
  rhs <- rlang::f_rhs(formula)
  new_lhs <- if (is.null(lhs)) NULL else .rewrite_panel_formula(lhs, state)
  new_rhs <- .rewrite_panel_formula(rhs, state)

  rewritten <- rlang::new_formula(new_lhs, new_rhs, env = env)
  mat <- .apply_ts_map(state$map, data, groupvar, timevar, env)

  list(data = mat, formula = rewritten, map = state$map)
}

# Build readable alias names for a ts_map: converts the deparsed source
# expression for each `.pt#` entry to a janitor-cleaned column name.
# Returns a named character vector `{".pt#" -> "lag_x"}`.
.pt_make_aliases <- function(ts_map) {
  if (length(ts_map) == 0L) return(stats::setNames(character(0L), character(0L)))
  exprs <- vapply(ts_map, deparse, character(1L), width.cutoff = 200L)
  dummy_df  <- stats::setNames(
    as.data.frame(rep(list(1L), length(exprs)), check.names = FALSE),
    exprs
  )
  raw_aliases <- names(janitor::clean_names(dummy_df))
  aliases     <- make.unique(raw_aliases, sep = "_")
  stats::setNames(aliases, names(ts_map))
}

# Recursively substitute `.pt#` symbols in a formula with their readable aliases.
.pt_alias_formula <- function(formula, alias_map) {
  if (length(alias_map) == 0L) return(formula)
  sub_pt <- function(expr) {
    if (is.symbol(expr)) {
      nm <- as.character(expr)
      if (nm %in% names(alias_map)) return(as.symbol(alias_map[[nm]]))
      return(expr)
    }
    if (!is.call(expr)) return(expr)
    as.call(lapply(as.list(expr), sub_pt))
  }
  lhs <- rlang::f_lhs(formula)
  rhs <- rlang::f_rhs(formula)
  rlang::new_formula(
    if (is.null(lhs)) NULL else sub_pt(lhs),
    sub_pt(rhs),
    env = rlang::f_env(formula)
  )
}

# Rewrite `formula` using the ts_map from a prior panel_materialize() call,
# without materialising any data. New ts sub-expressions not already in the map
# are extracted into additional .pt# slots (counter continues from the map
# length). This lets sub-formulas (mean/variance sides of heterolm) share the
# same materialization without a second `.apply_ts_map()` pass.
.rewrite_from_map <- function(formula, ts_map, ts_fns = NULL) {
  state <- new.env(parent = emptyenv())
  state$ts_fns  <- union(.pt_ts_fns, ts_fns)
  state$map     <- ts_map
  # Invert: deparse(expr) -> .pt# sym for dedup lookup
  state$keys    <- if (length(ts_map) == 0L) {
    list()
  } else {
    stats::setNames(
      as.list(names(ts_map)),
      vapply(ts_map, deparse, character(1L), width.cutoff = 500L)
    )
  }
  state$counter <- length(ts_map)

  lhs <- rlang::f_lhs(formula)
  rhs <- rlang::f_rhs(formula)
  new_lhs <- if (is.null(lhs)) NULL else .rewrite_panel_formula(lhs, state)
  new_rhs <- .rewrite_panel_formula(rhs, state)
  rlang::new_formula(new_lhs, new_rhs, env = rlang::f_env(formula))
}
