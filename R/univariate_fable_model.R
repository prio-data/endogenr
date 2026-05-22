#' Univariate model from fable
#'
#' Univariate statistical models ("univariate_fable") currently supports exponential smoothing (ETS) and ARIMA models. See <https://otexts.com/fpp3/> for
#' details on these models. These models are estimated independently for each groupvar in [setup_simulator()], and the forecasts are completely
#' independent of the rest of the system (the forecasts are populated in the simulation dataset before the dynamic simulation is calculated). See [fable::ETS()]
#' and [fable::ARIMA()] for how to write the function calls for these models. A simple exponential smoothing model can be set up using
#' build_model("univariate_fable", formula = dem ~ error("A") + trend("N") + season("N"), method = "ets").
#'
#' Requires packages `fable`, `fabletools`, and `tsibble` (in Suggests).
#'
#' @param formula An R-formula using the model specification in [fable::ARIMA()] or [fable::ETS()]
#' @param data A data.table, data.frame, or tsibble. Converted to tsibble internally for fable.
#' @param method Character. Either "arima" or "ets".
#' @param ctx A panel_context object.
#' @param ... Other arguments from [fable::ARIMA()] or [fable::ETS()]
#'
#' @return A univariate_fable endogenmodel class.
#' @export
#' @exportS3Method
fit_model.univariate_fable_spec <- function(spec, data = NULL, ctx = NULL, ...) {
  method <- if (!is.null(spec$args$method)) spec$args$method else "arima"
  extra_args <- spec$args[!names(spec$args) %in% "method"]
  do.call(univariate_fable_model, c(
    list(formula = spec$formula, data = data, method = method, ctx = ctx),
    extra_args
  ))
}

#' @export
univariate_fable_model <- function(formula = NULL, data = NULL, method = "arima", ctx = NULL, ...) {
  if (!requireNamespace("fable", quietly = TRUE) ||
      !requireNamespace("fabletools", quietly = TRUE) ||
      !requireNamespace("tsibble", quietly = TRUE)) {
    stop("Packages 'fable', 'fabletools', and 'tsibble' are required for univariate_fable models.")
  }

  model <- new_endogenmodel(formula)

  # Convert to tsibble for fable fitting
  grp <- ctx_unit(ctx)
  idx <- ctx_time(ctx)
  if (!inherits(data, "tbl_ts")) {
    ts_data <- tsibble::as_tsibble(as.data.frame(data), key = grp, index = idx)
  } else {
    ts_data <- data
  }

  if (method == "arima") {
    model$fitted <- ts_data |> fabletools::model(fable::ARIMA(formula, ...))
  } else if (method == "ets") {
    model$fitted <- ts_data |> fabletools::model(fable::ETS(formula, ...))
  } else {
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
#' @param data A data.table (the simulation grid).
#' @param ctx A panel_context object.
#' @param test_start Integer. Not used directly (horizon determines forecast length).
#' @param horizon the number of future steps to predict
#' @param inner_sims the number of inner simulations.
#' @param ...
#'
#' @return A data.table with columns: unit, sim, time, outcome.
#' @export
predict.univariate_fable <- function(model, data, ctx, test_start, horizon, inner_sims, ...) {
  grp <- ctx_unit(ctx)
  idx <- ctx_time(ctx)

  forecast <- model$fitted |>
    fabletools::forecast(h = paste(horizon, "years")) |>
    dplyr::mutate(samples = distributional::generate(!!rlang::sym(model$outcome), inner_sims)) |>
    dplyr::as_tibble() |>
    dplyr::select(dplyr::all_of(c(grp, idx, "samples"))) |>
    tidyr::unnest(samples) |>
    dplyr::mutate(sim = rep(1:inner_sims, dplyr::n() / inner_sims)) |>
    dplyr::rename(!!rlang::sym(model$outcome) := "samples")

  # Return as data.table
  data.table::as.data.table(forecast)
}
