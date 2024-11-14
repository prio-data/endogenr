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
#' @param model
#' @param horizon
#' @param test_start
#' @param data
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
