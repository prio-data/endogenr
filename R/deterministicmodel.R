
#' @exportS3Method
fit_model.deterministic_spec <- function(spec, ctx = NULL, ...) {
  deterministicmodel(formula = spec$formula, ctx = ctx)
}

#' Deterministic model
#'
#' This calculates a deterministic outcome based on an R formula.
#'
#' The formula's first RHS term must be wrapped in `I()`. Everything the
#' predict method needs repeatedly — the `I()` term label, the prepared
#' evaluation formula (index appended, positional time-series functions
#' injected), and the composed history depth of the RHS — is computed once
#' here so the hot simulation loop does no formula surgery.
#'
#' @param formula A two-sided formula of the form `outcome ~ I(expr)`.
#' @param ctx A panel_context object supplying the time column name.
#'
#' @return An endogenmodel of class `deterministic`.
#' @keywords internal
deterministicmodel <- function(formula = NULL, ctx = NULL){
  model <- new_endogenmodel(formula)
  class(model) <- c("deterministic", class(model))
  model$independent <- FALSE

  model$outcome <- parse_formula(model)$outcome

  # The transformed term label, validated once at build time.
  y_star <- attr(stats::terms(formula), "term.labels")[1L]
  if (is.na(y_star) || !grepl("^I\\(", y_star)) {
    stop("The formula to apply must be the first term and wrapped in I()",
         call. = FALSE)
  }
  model$y_star <- y_star

  # Prepared evaluation formula: index appended so the time column survives
  # model.frame(), positional lag/lead/shift/diff injected.
  frm <- stats::update(formula, paste(c(". ~ .", ctx_time(ctx)), collapse = "+"))
  model$mat_formula <- inject_positional_lag(frm)

  # Composed per-unit history depth of the RHS (Inf for cumulative terms),
  # used by predict to evaluate only the rows it needs instead of the full
  # grid at every step.
  model$required_history <- .required_history(formula)

  return(model)
}

#' Predict function for a deterministic model
#'
#' @param model A `deterministic` endogenmodel.
#' @param t Time step to predict.
#' @param data A data.table.
#' @param ctx A panel_context object.
#' @param ... Ignored, accepted for S3 generic consistency.
#'
#' @return A data.table with key + index + outcome columns.
#' @family simulation
#' @export
predict.deterministic <- function(model, t, data, ctx, ...) {
  idx <- ctx_time(ctx)
  all_keys <- ctx_keys(ctx)
  y_star <- model$y_star
  y <- model$outcome

  # Coerce to data.table if needed
  if (!data.table::is.data.table(data)) {
    data <- data.table::as.data.table(as.data.frame(data))
  }

  # Evaluate only the per-unit history window the RHS needs before `t`
  # (matching the other predict methods) instead of the entire grid — the
  # full-grid evaluation made every step O(grid), i.e. the loop O(T x grid).
  data <- .history_subset(data, idx, t, model$required_history)

  # Per-group model.frame — expose keys when the formula references them
  needs_keys <- any(all.vars(model$formula) %in% all_keys)
  result <- .model_frame_by_group(model$mat_formula, data, all_keys, needs_keys)

  # Remove the original outcome column, filter to time t, rename
  if (y %in% names(result)) {
    result[, (y) := NULL]
  }
  result <- .dt_rows(result, result[[idx]] == t)
  data.table::setnames(result, y_star, y)

  # Return only key + index + outcome
  result_cols <- c(all_keys, idx, y)
  result <- result[, ..result_cols]

  return(result)
}
