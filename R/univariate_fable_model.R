#' Univariate model from fable
#'
#' Univariate statistical models ("univariate_fable") currently supports exponential smoothing (ETS) and ARIMA models. See <https://otexts.com/fpp3/> for
#' details on these models. These models are estimated independently for each groupvar in [setup_simulator()], and the forecasts are completely
#' independent of the rest of the system (the forecasts are populated in the simulation dataset before the dynamic simulation is calculated). See [fable::ETS()]
#' and [fable::ARIMA()] for how to write the function calls for these models. A simple exponential smoothing model can be set up using
#' build_model("univariate_fable", formula = dem ~ error("A") + trend("N") + season("N"), method = "ets").
#'
#' @param formula An R-formula using the model specification in [fable::ARIMA()] or [fable::ETS()]
#' @param data The training data used to fit the model. A tsibble.
#' @param ... Other arguments from [fable::ARIMA()] or [fable::ETS()]
#'
#' @return A univariate_fable endogenmodel class.
#' @export
univariate_fable_model <- function(formula = NULL, data = NULL, method = "arima", ...){
  model <- new_endogenmodel(formula)

  if(method == "arima"){
    model$fitted <- data |> fabletools::model(fable::ARIMA(formula, ...))
  } else if(method == "ets"){
    model$fitted <- data |> fabletools::model(fable::ETS(formula, ...))
  } else{
    stop("Method not implemented.")
  }


  class(model) <- c("univariate_fable", class(model))
  model$independent <- TRUE

  model$outcome <- parse_formula(model)$outcome
  return(model)
}

#' Predict univariate fable model
#'
#' Generates predictions from a univariate_fable endogenmodel. Only use in dynamic simulation. If you want to use
#' separately, you should use fable directly!
#'
#' @param model a univariate_fable endogenmodel
#' @param horizon the number of future steps to predict
#' @param inner_sims the number of inner simulations. See [setup_simulator()].
#' @param data the data used in forecasting (used to gain information about group and index variables)
#'
#' @return
#' @export
#'
#' @examples
predict.univariate_fable <- function(model, horizon, inner_sims, data, ...){
  # Get index and key variables from tsibble
  idx <- tsibble::index_var(data)
  grp <- tsibble::key_vars(data)
  grp <- grp[!grp %in% "sim"]

  forecast <- model$fitted |> fabletools::forecast(h = paste(horizon, "years")) |>
    dplyr::mutate(samples = distributional::generate(!!rlang::sym(model$outcome), inner_sims)) |>
    dplyr::as_tibble() |>
    dplyr::select(dplyr::all_of(c(grp, idx, "samples"))) |>
    tidyr::unnest(samples) |>
    dplyr::mutate(sim = rep(1:inner_sims, n()/inner_sims)) |>
    dplyr::rename(!!rlang::sym(model$outcome) := "samples") |>
    tsibble::tsibble(key = c(grp, "sim"), index = idx)

  return(forecast)
}
