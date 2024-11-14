#' Exogenous model
#'
#' This just pastes in an already established quantitative trajectory into the prediction model.
#'
#' @param formula
#' @param impute_from
#' @param newdata
#' @param inner_sims
#'
#' @return
#' @export
#'
#' @examples
exogenmodel <- function(formula = NULL, impute_from = NULL, newdata = NULL, inner_sims = NULL){
  fit_exogenous_model <- function(formula, impute_from, data, inner_sims) {
    # Get key and index variables from tsibble
    grp <- tsibble::key_vars(data)
    idx <- tsibble::index_var(data)

    # Filter data based on impute_from
    filtered_data <- data |>
      dplyr::filter(!!rlang::sym(idx) >= impute_from)

    # Update formula to include group and index variables, removing intercept
    model_formula <- stats::update(
      formula,
      paste("~", grp, idx, ". - 1", sep = "+")
    )

    # Create model matrix
    result <- stats::model.matrix(model_formula, filtered_data) |>
      # Convert to tibble and add back key and index columns
      tibble::as_tibble() |>
      tsibble::as_tsibble(
        key = grp,
        index = idx
      )

    # Create sequence of all times
    all_times <- base::seq(min(result[[idx]]), max(result[[idx]]), by = 1)

    # Get unique groups
    all_groups <- tsibble::key_data(result)[[grp]]

    # Create expanded grid using tidyr
    group_sym <- rlang::sym(grp)
    time_sym <- rlang::sym(idx)
    expanded <- tidyr::expand_grid(
      !!time_sym := all_times,
      !!group_sym := all_groups,
      sim = 1:inner_sims
    )

    # Convert to tsibble
    result <- expanded |>
      tsibble::as_tsibble(
        key = c(grp, "sim"),
        index = idx
      ) |>
      dplyr::left_join(
        result,
        by = c(grp, idx)
      )

    # Return a fresh copy
    return(result)
  }

  grp <- tsibble::key_vars(newdata)
  idx <- tsibble::index_var(newdata)

  model <- new_endogenmodel(formula)
  model$impute_from <- impute_from

  model$fitted <- fit_exogenous_model(formula, impute_from, newdata, inner_sims)

  class(model) <- c("exogen", class(model))
  model$independent <- TRUE

  outcome <- parse_formula(model)$outcome
  return(model)
}


#' Prediction function for exogenous models
#'
#' @param model
#' @param data
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
predict.exogen <- function(model, data, ...){
  return(model$fitted)
}
