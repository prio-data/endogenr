#' Univariate model from fable
#'
#' Univariate statistical models ("univariate_fable") currently support
#' exponential smoothing (ETS) and ARIMA models. See <https://otexts.com/fpp3/>
#' for details on these models. These models are estimated independently for
#' each panel unit, and their forecasts are populated in the simulation
#' dataset before the dynamic simulation runs. See [fable::ETS()] and
#' [fable::ARIMA()] for how to write the function calls. A simple ETS model
#' is set up via
#' `build_model("univariate_fable", formula = y ~ error("A") + trend("N") + season("N"), method = "ets")`.
#'
#' Requires packages `fable`, `fabletools`, and `tsibble` (in Suggests).
#'
#' @param spec A `univariate_fable_spec` object from [build_model()]. The
#'   `formula` slot follows the model-specification syntax in [fable::ARIMA()]
#'   or [fable::ETS()].
#' @param data A data.table, data.frame, or tsibble. Converted to tsibble
#'   internally for fable.
#' @param ctx A panel_context object.
#' @param ... Additional arguments forwarded to [fable::ARIMA()] or
#'   [fable::ETS()].
#'
#' @return A univariate_fable endogenmodel.
#' @family simulation
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

#' @keywords internal
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
    ts_data <- tsibble::as_tsibble(as.data.frame(data),
                                   key = !!rlang::sym(grp),
                                   index = !!rlang::sym(idx))
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
#' Generates predictions from a univariate_fable endogenmodel. Only use in
#' dynamic simulation. If you want to use separately, you should use fable
#' directly!
#'
#' Sample paths are drawn with [fabletools::generate()], which simulates the
#' fitted ETS/ARIMA forward step by step with correlated innovations. Each
#' returned `(unit, sim)` series is therefore one internally-coherent
#' trajectory whose temporal structure (trend persistence, autocorrelation) is
#' preserved — unlike stitching independent per-horizon draws from the marginal
#' [fabletools::forecast()] distributions, which destroys that structure and
#' then feeds an incoherent path into the endogenous equations via `lag()`.
#'
#' @param model a univariate_fable endogenmodel
#' @param data A data.table (the simulation grid).
#' @param ctx A panel_context object.
#' @param test_start Integer. Not used directly (horizon determines forecast length).
#' @param horizon the number of future steps to predict
#' @param inner_sims the number of inner simulations.
#' @param ... Ignored, accepted for S3 generic consistency.
#'
#' @return A data.table with columns: unit, sim, time, outcome.
#' @family simulation
#' @export
predict.univariate_fable <- function(model, data, ctx, test_start, horizon, inner_sims, ...) {
  grp <- ctx_unit(ctx)
  idx <- ctx_time(ctx)
  sim_col <- if (is.null(ctx$sim)) "sim" else ctx$sim

  # Coherent sample paths: simulate `inner_sims` forward trajectories from the
  # fitted mable. `.rep` (the path id, "1".."inner_sims") becomes `sim`; `.sim`
  # (the simulated value) becomes the outcome column.
  sims <- fabletools::generate(model$fitted, h = horizon, times = inner_sims)

  forecast <- sims |>
    dplyr::as_tibble() |>
    dplyr::transmute(
      !!rlang::sym(grp)           := .data[[grp]],
      !!rlang::sym(idx)           := .data[[idx]],
      !!rlang::sym(sim_col)       := as.integer(.data[[".rep"]]),
      !!rlang::sym(model$outcome) := .data[[".sim"]]
    )

  data.table::as.data.table(forecast)
}
