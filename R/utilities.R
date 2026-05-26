#' Draw a random training window
#'
#' Used by `setup_simulator()` (and the inner refit loop in
#' `simulate_endogenr()`) to draw a random training window for `linear`,
#' `glm`, and `heterolm` models when `min_window` is set. The window is
#' anchored inside `[earliest_train_start, test_start - 1]`, with random
#' (possibly zero) padding on either side, and is guaranteed to be at
#' least `min_window` time steps long.
#'
#' @param earliest_train_start Integer. The earliest time step that may
#'   appear in the training window.
#' @param test_start Integer. The first time step in the forecast period;
#'   the training window always ends at `test_start - 1` or earlier.
#' @param min_window Integer or `NULL`. Minimum window length. When `NULL`
#'   the full range `[earliest_train_start, test_start - 1]` is returned
#'   unchanged.
#'
#' @return A list with `start`, `end`, and `window` (the latter `NULL`
#'   when `min_window` is `NULL`).
#' @family build
#' @export
#'
#' @examples
#' set.seed(1)
#' get_train_window(1970, 2010, min_window = 10)
get_train_window <- function(earliest_train_start, test_start, min_window = NULL){
  if(is.null(min_window)){
    return(list("start" = earliest_train_start, "end" = test_start -1, "window" = NULL))
  }
  full_range <- test_start - earliest_train_start
  if(min_window > full_range){
    stop("Train min_window must be smaller or equal to largest possible train-set.")
  }
  if(min_window < 1){
    stop("Train min_window must be 1 or larger")
  }

  start_increment <- sample.int(full_range-min_window, 1) - 1
  stop_decrement <- sample.int(full_range-min_window-start_increment, 1) - 1
  return(list("start" = earliest_train_start+start_increment, "end" = test_start - stop_decrement))
}

#' Inject a positional lag function into a formula's environment
#'
#' Replaces `lag()` in the formula's evaluation environment with a positional
#' shift function: `c(rep(NA, n), head(x, -n))`. This ensures `lag()` in model
#' formulas performs a within-group positional shift rather than `stats::lag()`.
#'
#' @param formula An R formula.
#' @return The formula with modified environment.
#' @family formula_helpers
#' @export
inject_positional_lag <- function(formula) {
  .positional_lag <- function(x, n = 1L) {
    c(rep(NA_real_, n), utils::head(x, -n))
  }
  formula_env <- new.env(parent = environment(formula))
  formula_env$lag <- .positional_lag
  environment(formula) <- formula_env
  formula
}

#' Materialize formula terms into a data.table
#'
#' Evaluates formula terms per group using `model.frame()` with positional lag
#' injection. Returns the materialized data.table with cleaned column names.
#'
#' For repeated calls with the same formula (e.g. inside a simulation loop),
#' pass `.mat_formula` and `.col_mapping` (cached from a prior call) to skip
#' the `update()` / `inject_positional_lag()` / `clean_names()` overhead.
#'
#' @param formula An R formula.
#' @param data A data.table (or data.frame/tsibble — will be coerced).
#' @param ctx A `panel_context` object.
#' @param .mat_formula Optional. A pre-prepared formula (already has index var
#'   appended and positional lag injected). When provided, `formula` is ignored
#'   for evaluation and used only as documentation.
#' @param .col_mapping Optional named character vector. Maps raw `model.frame()`
#'   column names to clean names. When provided, `janitor::clean_names()` is
#'   skipped and `setnames()` is used instead.
#'
#' @return A data.table with cleaned column names.
#' @family formula_helpers
#' @export
materialize_formula <- function(formula, data, ctx,
                                .mat_formula = NULL, .col_mapping = NULL) {
  all_keys <- ctx_keys(ctx)

  if (!data.table::is.data.table(data)) {
    data <- data.table::as.data.table(as.data.frame(data))
  }

  if (is.null(.mat_formula)) {
    idx <- ctx_time(ctx)
    .mat_formula <- stats::update(formula, paste(c(". ~ .", idx), collapse = "+"))
    .mat_formula <- inject_positional_lag(.mat_formula)
  }

  eval_formula <- .mat_formula
  result <- data[, {
    mf <- stats::model.frame(eval_formula, data = .SD, na.action = stats::na.pass)
    data.table::as.data.table(mf)
  }, by = all_keys]

  if (!is.null(.col_mapping)) {
    old <- intersect(names(.col_mapping), names(result))
    if (length(old) > 0L) {
      data.table::setnames(result, old, .col_mapping[old])
    }
  } else {
    result <- janitor::clean_names(result)
  }

  result
}

#' Derive a naive formula from materialized column names
#'
#' Given a materialized data.table (or its column names), determines the outcome
#' and predictor columns by excluding key and index columns, and returns a
#' formula suitable for `lm()` / `glm()`.
#'
#' @param data A data.table, or a character vector of column names.
#' @param outcome Character. Explicit outcome column name. If `NULL`, the first
#'   non-key column is used (positional convention from `model.frame()`).
#' @param ctx A `panel_context` object.
#'
#' @return An R formula.
#' @family formula_helpers
#' @export
derive_naive_formula <- function(data, outcome = NULL, ctx) {
  idx <- ctx_time(ctx)
  all_keys <- ctx_keys(ctx)

  col_names <- if (is.character(data)) data else names(data)
  varnames <- setdiff(col_names, c(all_keys, idx))

  if (is.null(outcome)) {
    outcome <- varnames[1]
    predictors <- varnames[-1]
  } else {
    predictors <- setdiff(varnames, outcome)
  }

  stats::reformulate(predictors, response = outcome)
}

#' Creates a panel data frame based on a formula and data.table
#'
#' Convenience wrapper that calls [materialize_formula()] followed by
#' [derive_naive_formula()]. Returns both the materialized data and the
#' naive formula.
#'
#' @param formula An R formula.
#' @param data A data.table (or data.frame/tsibble — will be coerced).
#' @param ctx A `panel_context` object.
#'
#' @return A list with `data` (data.table) and `naive_formula`.
#' @keywords internal
create_panel_frame <- function(formula, data, ctx) {
  materialized <- materialize_formula(formula, data, ctx)
  naive_formula <- derive_naive_formula(materialized, ctx = ctx)
  list(data = materialized, naive_formula = naive_formula)
}

#' Build a prepared formula and column-name mapping for fast materialization
#'
#' Pre-computes the formula (with index appended + positional lag injected) and
#' the raw → clean column name mapping. These can be passed to
#' [materialize_formula()] via `.mat_formula` and `.col_mapping` to skip
#' per-call overhead.
#'
#' @param formula An R formula.
#' @param data A data.table with at least a few rows (used to discover column
#'   names from `model.frame()`).
#' @param ctx A `panel_context` object.
#'
#' @return A list with `mat_formula` and `col_mapping`.
#' @keywords internal
.build_mat_cache <- function(formula, data, ctx) {
  idx <- ctx_time(ctx)
  all_keys <- ctx_keys(ctx)

  mat_formula <- stats::update(formula, paste(c(". ~ .", idx), collapse = "+"))
  mat_formula <- inject_positional_lag(mat_formula)

  # Run model.frame on a small per-group sample to discover raw column names
  sample_dt <- data[, utils::head(.SD, 2L), by = all_keys]

  raw <- sample_dt[, {
    mf <- stats::model.frame(mat_formula, data = .SD, na.action = stats::na.pass)
    data.table::as.data.table(mf)
  }, by = all_keys]

  raw_names <- names(raw)
  cleaned <- janitor::clean_names(raw)
  clean_names <- names(cleaned)

  # Only map names that actually changed
  changed <- raw_names != clean_names
  col_mapping <- stats::setNames(clean_names[changed], raw_names[changed])

  list(mat_formula = mat_formula, col_mapping = col_mapping)
}
