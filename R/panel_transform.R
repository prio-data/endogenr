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

# Time-series functions: extracted to a per-group `.pt#` column. Matched by
# bare name or `pkg::fn`. The whole call (e.g. the entire `lag(...)`) is
# materialised, so anything wrapped in one is always evaluated within group.
# Callers may extend this set per fit via the `ts_fns` argument of
# panel_materialize() to register custom within-group functions.
.pt_ts_fns <- c(
  "lag", "lead", "lead_horizon", "shift",
  "rollmean", "rollsum", "rollapply", "frollmean", "frollsum", "frollapply",
  "cumsum", "cumprod", "cummax", "cummin", "diff",
  "decay_since_event", "time_since_event", "intensity_decay"
)

# Design operators: NOT extracted. The call head is kept and its arguments are
# recursed into, so time-series terms nested inside them are still extracted
# (e.g. `poly(lag(x), 2)` -> `poly(.pt1, 2)`). Listed for documentation; the
# rewrite treats them identically to any other non-time-series call.
.pt_design_ops <- c(
  "factor", "ordered", "poly", "ns", "bs", "cut", "C", "interaction"
)

# Positional within-group lag, matching inject_positional_lag(): `lag(x, n)`
# shifts `x` forward by `n` rows, padding the head with NA.
.pt_positional_lag <- function(x, n = 1L) {
  c(rep(NA_real_, n), utils::head(x, -n))
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

# Materialise a `.pt#` map into per-group, time-ordered columns. Returns a copy
# of `data` (keyed by group then time) with one new column per map entry. `env`
# is the formula environment used to resolve time-series functions other than
# `lag` (which is overridden with the positional shift).
.apply_ts_map <- function(map, data, groupvar, timevar, env = parent.frame()) {
  if (!data.table::is.data.table(data)) {
    data <- data.table::as.data.table(as.data.frame(data))
  }
  data <- data.table::copy(data)
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
