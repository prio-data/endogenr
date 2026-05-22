#' Exogenous model
#'
#' Creates a model that copies pre-existing data (e.g. population projections)
#' into the simulation grid. The exogenous variable must already be present in
#' the data passed to [setup_simulator()].
#'
#' @param formula A one-sided formula with the exogenous variable (e.g. `~population`).
#' @param impute_from Integer. The time from which to take values (typically `test_start`).
#' @param newdata A data.frame or data.table containing the exogenous variable for
#'   the full forecast period.
#' @param inner_sims Integer. Number of inner simulations.
#' @param ctx A `panel_context` object.
#'
#' @return An endogenmodel of class `exogen`.
#' @export
#' @exportS3Method
fit_model.exogen_spec <- function(spec, newdata = NULL, ctx = NULL,
                                  test_start = NULL, inner_sims = NULL, ...) {
  exogenmodel(formula = spec$formula, impute_from = test_start,
              newdata = newdata, inner_sims = inner_sims, ctx = ctx)
}

#' @keywords internal
exogenmodel <- function(formula = NULL, impute_from = NULL, newdata = NULL,
                        inner_sims = NULL, ctx = NULL) {
  model <- new_endogenmodel(formula)
  model$impute_from <- impute_from

  class(model) <- c("exogen", class(model))
  model$independent <- TRUE
  model$outcome <- parse_formula(model)$outcome

  # Store the source data for predict-time expansion

  grp <- ctx_unit(ctx)
  idx <- ctx_time(ctx)

  if (!data.table::is.data.table(newdata)) {
    newdata <- data.table::as.data.table(as.data.frame(newdata))
  }

  # Keep only the columns we need: unit, time, outcome variable(s)
  outcome_vars <- all.vars(formula)
  keep_cols <- c(grp, idx, outcome_vars)
  model$source_data <- newdata[newdata[[idx]] >= impute_from, ..keep_cols]

  return(model)
}


#' Prediction function for exogenous models
#'
#' Expands the stored source data across all simulation draws and returns a
#' data.table keyed by (unit, sim, time).
#'
#' @param model An `exogen` endogenmodel.
#' @param data A data.table (the simulation grid, used to extract sim IDs).
#' @param ctx A `panel_context` object.
#' @param test_start Integer. Start of the forecast period.
#' @param horizon Integer. Number of forecast steps.
#' @param inner_sims Integer. Number of inner simulations.
#' @param ... Ignored.
#'
#' @return A data.table with columns: unit, sim, time, and outcome variable(s).
#' @export
predict.exogen <- function(model, data, ctx, test_start, horizon, inner_sims, ...) {
  grp <- ctx_unit(ctx)
  idx <- ctx_time(ctx)

  source <- model$source_data

  # Filter source to the relevant forecast window
  source <- source[source[[idx]] >= test_start & source[[idx]] <= (test_start + horizon - 1)]

  # Expand across sim dimension via CJ
  all_units <- unique(source[[grp]])
  all_times <- unique(source[[idx]])

  grid <- data.table::CJ(
    unit_placeholder = all_units,
    time_placeholder = all_times,
    sim = seq_len(inner_sims),
    sorted = FALSE
  )
  data.table::setnames(grid, c("unit_placeholder", "time_placeholder"), c(grp, idx))

  # Merge source values onto grid (same value for every sim)
  result <- merge(grid, source, by = c(grp, idx), all.x = TRUE)

  return(result)
}
