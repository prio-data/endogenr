#' Univariate model from fable
#'
#' @param formula
#' @param data
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
univariate_fable_model <- function(formula = NULL, data = NULL, ...){
  model <- new_endogenmodel(formula)

  model$fitted <- df |> fabletools::model(fable::ARIMA(formula))


  class(model) <- c("univariate_fable", class(model))
  model$independent <- TRUE

  outcome <- parse_formula(model)$outcome
  return(model)
}

#' Predict univariate fable model
#'
#' @param model
#' @param horizon
#' @param test_start
#' @param data
#'
#' @return
#' @export
#'
#' @examples
predict.univariate_fable <- function(model, horizon, test_start, data){
  # Get index and key variables from tsibble
  idx <- tsibble::index_var(data)
  grp <- tsibble::key_vars(data)

  forecast <- model$fitted |> fabletools::forecast(h = paste(horizon, "years")) |>
    dplyr::mutate(sample = distributional::generate(!!rlang::sym(model$outcome), 1) |> unlist()) |>
    dplyr::select(dplyr::all_of(c(grp, idx, "sample"))) |>
    tsibble::tsibble(key = grp, index = grp)

  # Create prediction data frame with only necessary columns
  pred_data <- data |>
    dplyr::filter(!!rlang::sym(idx) >= test_start) |>
    dplyr::select(
      !!!rlang::syms(grp),
      !!rlang::sym(idx),
      !!rlang::sym(model$outcome)
    ) |>
    tsibble::tsibble(key = grp, index = grp)

  pred_data <- dplyr::left_join(pred_data, forecast)

  return(pred_data)
}
