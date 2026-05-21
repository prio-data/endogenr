#' Draw a random time-window
#'
#' @param earliest_train_start
#' @param test_start
#' @param min_window
#'
#' @return
#' @export
#'
#' @examples
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

#' Creates a panel data frame based on a formula and data.table
#'
#' Evaluates formula terms per group using `model.frame()` with positional lag
#' injection. Returns materialized data and a naive formula for model fitting.
#'
#' @param formula An R formula.
#' @param data A data.table (or data.frame/tsibble — will be coerced).
#' @param ctx A `panel_context` object.
#'
#' @return A list with `data` (data.table) and `naive_formula`.
#' @export
create_panel_frame <- function(formula, data, ctx) {
  grp <- ctx_unit(ctx)
  idx <- ctx_time(ctx)
  sim <- ctx_sim(ctx)
  all_keys <- ctx_keys(ctx)

  # Coerce to data.table if needed
  if (!data.table::is.data.table(data)) {
    data <- data.table::as.data.table(as.data.frame(data))
  }

  # Update formula to include index variable
  formula <- stats::update(formula, paste(c(". ~ .", idx), collapse = "+"))

  # Inject positional lag
  formula <- inject_positional_lag(formula)

  # Per-group model.frame
  result <- data[, {
    mf <- stats::model.frame(formula, data = .SD, na.action = stats::na.pass)
    data.table::as.data.table(mf)
  }, by = all_keys]

  # Clean names
  result <- janitor::clean_names(result)

  # Derive naive formula: first non-key column is outcome, rest are RHS
  varnames <- setdiff(names(result), c(all_keys, idx))
  naive_formula <- stats::reformulate(varnames[-1], response = varnames[1])

  return(list(
    "data" = result,
    "naive_formula" = naive_formula
  ))
}
