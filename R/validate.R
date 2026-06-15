#' Validate panel data integrity for simulation
#'
#' Checks that the input data meets the structural requirements for the endogenr
#' simulation system: integer time steps, contiguous series per unit, complete
#' initial state, and correct sort order.
#'
#' Panels may be unbalanced: each unit's series must be contiguous with integer
#' time steps, but units may enter late or exit early. Units with no row at
#' `test_start - 1` are used for fitting only and are excluded from the
#' simulation; all units present at `test_start - 1` must have a complete
#' initial state in every modelled outcome.
#'
#' @param data A data.table.
#' @param ctx A panel_context object.
#' @param test_start Integer. The first forecast time step.
#' @param model_outcomes Character vector or NULL. Outcome variable names from
#'   model formulas. If provided, checks that these columns have no NAs at the
#'   last pre-forecast time step.
#'
#' @return Invisible TRUE if all checks pass. Throws informative errors otherwise.
#' @family validate
#' @export
validate_panel <- function(data, ctx, test_start, model_outcomes = NULL) {
  unit_col <- ctx_unit(ctx)
  time_col <- ctx_time(ctx)

  # 1. Required columns exist
  missing_cols <- setdiff(c(unit_col, time_col), names(data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "),
         call. = FALSE)
  }

  times <- data[[time_col]]

  # 2. Integer (or integer-like) time steps
  if (!is.numeric(times)) {
    stop("Time column '", time_col, "' must be numeric, got ", class(times)[1],
         call. = FALSE)
  }
  if (!isTRUE(all.equal(times, as.integer(times)))) {
    stop("Time column '", time_col, "' must contain integer values",
         call. = FALSE)
  }

  # 3. Contiguous time series per unit
  units <- unique(data[[unit_col]])
  for (u in units) {
    unit_times <- sort(data[[time_col]][data[[unit_col]] == u])
    diffs <- diff(unit_times)
    if (any(diffs != 1L)) {
      gaps <- which(diffs != 1L)
      gap_at <- unit_times[gaps[1]]
      stop("Non-contiguous time series for unit '", u,
           "': gap after time step ", gap_at,
           call. = FALSE)
    }
  }

  # 4. At least one unit must have data at the forecast origin
  origin_units <- unique(data[[unit_col]][data[[time_col]] == (test_start - 1L)])
  if (length(origin_units) == 0L) {
    stop("No units have data at the forecast origin (test_start - 1 = ",
         test_start - 1L, "). Cannot simulate.", call. = FALSE)
  }

  # 5. Complete initial state (no NAs at test_start - 1 for model outcomes)
  if (!is.null(model_outcomes)) {
    existing_outcomes <- intersect(model_outcomes, names(data))
    if (length(existing_outcomes) > 0) {
      initial <- data[data[[time_col]] == (test_start - 1L)]
      if (nrow(initial) > 0) {
        for (col in existing_outcomes) {
          na_units <- initial[[unit_col]][is.na(initial[[col]])]
          if (length(na_units) > 0) {
            stop("Missing initial state at time ", test_start - 1L,
                 " for outcome '", col, "' in units: ",
                 paste(na_units, collapse = ", "),
                 call. = FALSE)
          }
        }
      }
    }
  }

  # 6. Data must be sorted by (unit, time)
  expected_order <- order(data[[unit_col]], data[[time_col]])
  if (!identical(seq_len(nrow(data)), expected_order)) {
    stop("Data must be sorted by (", unit_col, ", ", time_col, "). ",
         "Use data.table::setkeyv(data, c('", unit_col, "', '", time_col, "'))",
         call. = FALSE)
  }

  invisible(TRUE)
}


#' Validate that the model system is closed
#'
#' Checks that every predictor referenced by a model formula is produced by some
#' model (so it is available throughout the forecast window), that every
#' formula-referenced variable exists as a data column, and that no two models
#' share an outcome.
#'
#' @param models A list of models (unfitted specs or fitted model objects with
#'   a `$formula` field).
#' @param data_columns Character vector. Column names available in the data.
#'
#' @return Invisible TRUE if valid. Throws informative errors otherwise.
#' @family validate
#' @export
validate_system_closure <- function(models, data_columns) {
  # Collect all model outcomes
  independent_types <- get_independent_models()
  outcomes <- character(0)
  all_rhs_vars <- character(0)
  all_formula_vars <- character(0)

  for (i in seq_along(models)) {
    model <- models[[i]]

    # Support specs ($formula directly), fitted models ($formula), and legacy partials (attr)
    formula <- model$formula
    if (is.null(formula)) formula <- attr(model, "formula")

    # Determine if independent: check spec$type, then class
    is_independent <- FALSE
    if (!is.null(model$type)) {
      is_independent <- model$type %in% independent_types
    } else {
      is_independent <- any(independent_types %in% class(model))
    }

    if (is_independent) {
      outcome <- all.vars(formula)
      all_formula_vars <- c(all_formula_vars, outcome)
    } else {
      outcome <- all.vars(rlang::f_lhs(formula))
      rhs_vars <- all.vars(rlang::f_rhs(formula))
      all_rhs_vars <- c(all_rhs_vars, rhs_vars)
      all_formula_vars <- c(all_formula_vars, outcome, rhs_vars)

      # Also check variance formula for heterolm
      var_formula <- if (!is.null(model$args)) model$args$variance else model$variance_formula
      is_heterolm <- (!is.null(model$type) && model$type == "heterolm") ||
                     ("heterolm" %in% class(model))
      if (is_heterolm && !is.null(var_formula)) {
        var_rhs <- all.vars(rlang::f_rhs(var_formula))
        all_rhs_vars <- c(all_rhs_vars, var_rhs)
        all_formula_vars <- c(all_formula_vars, var_rhs)
      }
    }

    outcomes <- c(outcomes, outcome)
  }

  # Check for duplicate outcomes
  dup_outcomes <- outcomes[duplicated(outcomes)]
  if (length(dup_outcomes) > 0) {
    stop("Duplicate model outcomes: ", paste(unique(dup_outcomes), collapse = ", "),
         ". Each outcome variable must be modelled by exactly one model.",
         call. = FALSE)
  }

  # Every formula-referenced variable (LHS, RHS, variance) must exist as a data column.
  # Outcomes are produced at simulation time, but the column must already be present
  # so predict-time model.frame() and update-joins have somewhere to write.
  all_formula_vars <- unique(all_formula_vars)
  missing_in_data <- setdiff(all_formula_vars, data_columns)
  if (length(missing_in_data) > 0) {
    stop("The following variables are referenced by model formulas but missing ",
         "from the input data: ", paste(missing_in_data, collapse = ", "),
         ". Add them as columns (use NA for purely derived outputs).",
         call. = FALSE)
  }

  # System closure: every RHS predictor must be produced by a model. The
  # simulation grid carries only model outputs into the forecast window, so a
  # variable referenced as a predictor but produced by no model is dropped to NA
  # after the training period and silently propagates NA through the forecast —
  # immediately when used at the current period, after `n` steps when used at
  # lag `n`. (Lagged self-references such as lag(y) in y's own model are fine: y
  # is its own outcome.)
  all_rhs_vars <- unique(all_rhs_vars)
  unmodeled <- setdiff(all_rhs_vars, outcomes)
  if (length(unmodeled) > 0) {
    stop("The following variables are referenced as predictors but are not ",
         "produced by any model: ", paste(unmodeled, collapse = ", "),
         ". Every predictor must be produced by a model so it is available ",
         "throughout the forecast window. Add a model that produces each, e.g. ",
         "build_model('exogen', formula = ~", unmodeled[1], ").",
         call. = FALSE)
  }

  invisible(TRUE)
}
