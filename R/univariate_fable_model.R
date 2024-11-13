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

  model$outcome <- parse_formula(model)$outcome
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
predict.univariate_fable <- function(model, horizon, test_start, inner_sims, data){
  # Get index and key variables from tsibble
  idx <- tsibble::index_var(data)
  grp <- tsibble::key_vars(data)

  forecast <- model$fitted |> fabletools::forecast(h = paste(horizon, "years")) |>
    dplyr::mutate(samples = distributional::generate(!!rlang::sym(model$outcome), inner_sims)) |>
    dplyr::as_tibble() |>
    dplyr::select(dplyr::all_of(c(grp, idx, "samples"))) |>
    tidyr::unnest(samples) |>
    dplyr::mutate(.sim = rep(1:inner_sims, n()/inner_sims)) |>
    dplyr::rename(!!rlang::sym(model$outcome) := "samples") |>
    tsibble::tsibble(key = c(grp, ".sim"), index = idx)

  return(forecast)
}
