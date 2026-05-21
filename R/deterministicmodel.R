
#' Deterministic model
#'
#' This calculates a deterministic outcome based on an R formula.
#'
#' @param formula
#' @param ctx A panel_context object (unused but accepted for consistency).
#'
#' @return
#' @export
#'
#' @examples
deterministicmodel <- function(formula = NULL, ctx = NULL){
  model <- new_endogenmodel(formula)
  class(model) <- c("deterministic", class(model))
  model$independent <- FALSE

  model$outcome <- parse_formula(model)$outcome
  return(model)
}

#' Predict function for a deterministic model
#'
#' @param model
#' @param t Time step to predict.
#' @param data A data.table.
#' @param ctx A panel_context object.
#' @param ...
#'
#' @return A data.table with key + index + outcome columns.
#' @export
predict.deterministic <- function(model, t, data, ctx, ...) {
  idx <- ctx_time(ctx)
  all_keys <- ctx_keys(ctx)

  # Get the transformed variable name (y_star) and original outcome variable (y)
  y_star <- attr(stats::terms(model$formula)[1], "term.labels")
  y <- model$outcome

  if (!grepl("^I\\(", y_star)) {
    stop("The formula to apply must be the first term and wrapped in I()")
  }

  # Update formula to include index variable
  frm <- stats::update(model$formula, paste(c(". ~ .", idx), collapse = "+"))

  # Inject positional lag
  frm <- inject_positional_lag(frm)

  # Coerce to data.table if needed
  if (!data.table::is.data.table(data)) {
    data <- data.table::as.data.table(as.data.frame(data))
  }

  # Per-group model.frame
  result <- data[, {
    mf <- stats::model.frame(frm, data = .SD, na.action = stats::na.pass)
    data.table::as.data.table(mf)
  }, by = all_keys]

  # Remove the original outcome column, filter to time t, rename
  if (y %in% names(result)) {
    result[, (y) := NULL]
  }
  result <- result[result[[idx]] == t]
  data.table::setnames(result, y_star, y)

  # Return only key + index + outcome
  result_cols <- c(all_keys, idx, y)
  result <- result[, ..result_cols]

  return(result)
}
