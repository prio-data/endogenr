
#' Deterministic model
#'
#' This calculates a deterministic outcome based on an R formula.
#'
#' @param formula
#'
#' @return
#' @export
#'
#' @examples
deterministicmodel <- function(formula = NULL){
  model <- new_endogenmodel(formula)
  class(model) <- c("deterministic", class(model))
  model$independent <- FALSE

  model$outcome <- parse_formula(model)$outcome
  return(model)
}

#' Predict function for a deterministic model
#'
#' @param model
#' @param t
#' @param data
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
predict.deterministic <- function(model, t, data, ...) {
  # Get key and index variables from tsibble
  grp <- tsibble::key_vars(data)
  idx <- tsibble::index_var(data)

  # Get the transformed variable name (y_star) and original outcome variable (y)
  y_star <- attr(stats::terms(model$formula)[1], "term.labels")
  y <- model$outcome

  # Check if the formula is properly wrapped in I()
  if (!stringr::str_detect(y_star, "I\\(")) {
    stop("The formula to apply must be the first term and wrapped in I()")
  }

  # Update formula to include index variable
  frm <- stats::update(model$formula, paste(c(". ~ .", idx), collapse = "+"))

  # Create model frame by group
  result <- data %>%
    dplyr::group_by(!!!rlang::syms(grp)) %>%
    dplyr::group_modify(~ {
      mf <- stats::model.frame(frm, data = ., na.action = NULL)
      tibble::as_tibble(mf)
    }) %>%
    dplyr::ungroup() %>%
    # Remove the original outcome column if it exists
    dplyr::select(-!!rlang::sym(y)) %>%
    # Filter for specific time point
    dplyr::filter(!!rlang::sym(idx) == t) %>%
    # Rename the transformed variable to the original outcome name
    dplyr::rename(!!rlang::sym(y) := !!rlang::sym(y_star)) %>%
    # Convert to tsibble
    tsibble::as_tsibble(
      key = grp,
      index = idx
    )
  return(result)
}
