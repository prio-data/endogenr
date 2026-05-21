#' Validate panel data integrity for simulation
#'
#' Checks that the input data meets the structural requirements for the endogenr
#' simulation system: integer time steps, contiguous series, balanced panel,
#' complete initial state, and correct sort order.
#'
#' @param data A data.table.
#' @param ctx A panel_context object.
#' @param test_start Integer. The first forecast time step.
#' @param model_outcomes Character vector or NULL. Outcome variable names from
#'   model formulas. If provided, checks that these columns have no NAs at the
#'   last pre-forecast time step.
#'
#' @return Invisible TRUE if all checks pass. Throws informative errors otherwise.
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

  # 4. Balanced panel (all units have the same time range)
  unit_ranges <- tapply(data[[time_col]], data[[unit_col]], range, simplify = FALSE)
  ref_range <- unit_ranges[[1]]
  for (u in names(unit_ranges)) {
    if (!identical(unit_ranges[[u]], ref_range)) {
      stop("Unbalanced panel: unit '", u, "' has time range [",
           unit_ranges[[u]][1], ", ", unit_ranges[[u]][2],
           "] but expected [", ref_range[1], ", ", ref_range[2], "]",
           call. = FALSE)
    }
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
#' Checks that every RHS variable referenced by model formulas is available
#' either as another model's outcome or as a column in the data. Also checks
#' for duplicate outcomes across models.
#'
#' @param models A list of models (unfitted specs or fitted model objects with
#'   a `$formula` field).
#' @param data_columns Character vector. Column names available in the data.
#'
#' @return Invisible TRUE if valid. Throws informative errors otherwise.
#' @export
validate_system_closure <- function(models, data_columns) {
  # Collect all model outcomes
  independent_types <- get_independent_models()
  outcomes <- character(0)
  all_rhs_vars <- character(0)

  for (i in seq_along(models)) {
    model <- models[[i]]
    # Support both build_model partials (attr) and fitted models ($formula)
    formula <- attr(model, "formula")
    if (is.null(formula)) formula <- model$formula

    if (any(independent_types %in% class(model))) {
      outcome <- all.vars(formula)
    } else {
      outcome <- all.vars(rlang::f_lhs(formula))
      rhs_vars <- all.vars(rlang::f_rhs(formula))
      all_rhs_vars <- c(all_rhs_vars, rhs_vars)

      # Also check variance formula for heterolm
      var_formula <- attr(model, "variance_formula")
      if (is.null(var_formula) && is.list(model)) var_formula <- model$variance_formula
      if ("heterolm" %in% class(model) && !is.null(var_formula)) {
        all_rhs_vars <- c(all_rhs_vars, all.vars(rlang::f_rhs(var_formula)))
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

  # Check system closure: every RHS var must be in {outcomes} U {data_columns}
  all_rhs_vars <- unique(all_rhs_vars)
  available <- unique(c(outcomes, data_columns))
  missing_vars <- setdiff(all_rhs_vars, available)

  if (length(missing_vars) > 0) {
    stop("System is not closed. The following RHS variables are not provided ",
         "by any model or data column: ", paste(missing_vars, collapse = ", "),
         call. = FALSE)
  }

  invisible(TRUE)
}
