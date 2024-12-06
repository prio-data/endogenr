#' A baseline structure for a model
#'
#' Do not use directly.
#'
#' @param formula
#'
#' @return A list with class endogenmodel
#' @export
new_endogenmodel <- function(formula){
  structure(
    list(
      formula = formula
    ),
    class = "endogenmodel"
  )
}

#' Build a model
#'
#' This function is used to set up a model that will become part of a simulation-system. endogenr currently supports
#' five different model types: deterministic models (mathematical functions), parametric distributions (using fitdistrplus),
#' linear models (using lm), exogenous inputs (pre-calculated data), and ETS/AARIMA univariate statistical models (using fable).
#'
#' Deterministic models (type = "deterministic") are built using the R formula language, with the resulting variable on the left hand side, and a function
#' potentially using several functions, variables or constants. For these models to work, you will have to wrap the full right hand side
#' using [base::I()]. E.g., build_model("deterministic", formula = gdp ~ I(abs(gdppc*population))).
#'
#' Parametric distributions (type = "parametric_distribution") are fitted using fitdistrplus. The formula must only contain a single right hand side variable.
#' The model uses the full pool of historical observations to fit a global distribution that will be used to forecast all subsequent
#' observations. endogenr uses [fitdistrplus::fitdist()], but the supported distributions is currently limited to "norm" ([stats::pnorm()]),
#' "cauchy" ([stats::pcauchy()]), "gumbel" ([actuar::pgumbel()]), "gamma" ([stats::pgamma()]), "t_ls" ([pt_ls()]), "nbinom" ([stats::pnbinom()]),
#' and "poisson" ([stats::ppoisson()]). Support for other distributions must be coded into [create_distribution_object()], as endogenr is
#' using the distributional package to store distributions and generate forecasts from them. Certain distributions will require start values,
#' or different fitting methods. start, method, discrete, and other arguments to [fitdistrplus::fitdist()] can be added directly to
#' [build_model()]. E.g., build_model("parametric_distribution", formula = ~gdppc_grwt, distribution = "t_ls", start = list(df = 1, mu = mean(train$gdppc_grwt), sigma = sd(train$gdppc_grwt)))
#'
#' Linear models (type = "linear) are using [stats::lm()] to fit the models. However, if you set boot to "resid" or "wild", then a bootstrap method ([bootstraplm()]) is
#' used instead. This function supports residual (residuals are resampled) or wild bootstrap (residuals are multiplied with a Standard normal distribution).
#' Linear models are normally using the full training data to train the model. However, if "min_window" is set in [setup_simulator()],
#' then a random training interval with window at least "min_window" is drawn for each outer simulation ("nsim" in [simulate_endogenr()]) and used
#' as the training set for the linear model. Only linear models are affected by "min_window". Predictions from the linear model are drawn
#' from the full predictive distribution. See [getpi()] and [get_sepi()] for the implementations. The formula in the linear model
#' must consider whether you have the data dynamically available for simulation at any time-step. Dynamically estimated outcomes that
#' are used as input-variables must be lagged using [dplyr::lag()]. Transformations and time-series functions (like [zoo::rollapplyr()])
#' works. Remember to set align = "right" when using these functions, and bear in mind that you will loose some training data when
#' using these functions. E.g., build_model("linear", formula = yjbest ~ lag(zoo::rollsumr(yjbest, k = 5, fill = NA)) + lag(log(gdppc)), boot = "resid").
#' Hierarchical models are currently not supported.
#'
#' Exogenous input (type = "exogen") requires a [tsibble::tsibble()] data.frame where the key and index is the
#' same as "groupvar" and "timevar" in [setup_simulator()], respectively. The formula must only contain a single right hand side variable.
#' Exogenous input must be merged into the data in [setup_simulator()] for the full forecast horizon. E.g., build_model("exogen", formula = ~population)
#'
#' Univariate statistical models ("univariate_fable") currently supports exponential smoothing (ETS) and ARIMA models. See <https://otexts.com/fpp3/> for
#' details on these models. These models are estimated independently for each groupvar in [setup_simulator()], and the forecasts are completely
#' independent of the rest of the system (the forecasts are populated in the simulation dataset before the dynamic simulation is calculated). See [fable::ETS()]
#' and [fable::ARIMA()] for how to write the function calls for these models. A simple exponential smoothing model can be set up using
#' build_model("univariate_fable", formula = dem ~ error("A") + trend("N") + season("N"), method = "ets").
#'
#' @param type One of deterministic, parametric_distribution, linear, exogen, or univariate_fable. See details.
#' @param formula An R formula. See details.
#' @param ...
#'
#' @return A partialised function used by the simulation system. See [purrr::partial()] for details.
#' @export
#'
#' @examples
#' df <- endogenr::example_data |> tsibble::as_tsibble(key = "gwcode", index = "year")
#' train <- df |> dplyr::filter(year>= 1970, year < 2010) # used for starting values in parametric distribution
#' c1 <- yjbest ~ lag(zoo::rollsumr(yjbest, k = 5, fill = NA)) + lag(log(gdppc))
#' model_system <- list(
#'   build_model("deterministic",formula = gdppc ~ I(abs(lag(gdppc)*(1+gdppc_grwt)))),
#'   build_model("deterministic", formula = gdp ~ I(abs(gdppc*population))),
#'   build_model("parametric_distribution", formula = ~gdppc_grwt, distribution = "t_ls", start = list(df = 1, mu = mean(train$gdppc_grwt), sigma = sd(train$gdppc_grwt))),
#'   build_model("linear", formula = c1, boot = "resid"),
#'   build_model("univariate_fable", formula = dem ~ error("A") + trend("N") + season("N"), method = "ets"),
#'   build_model("exogen", formula = ~psecprop),
#'   build_model("exogen", formula = ~population)
#' )
build_model <- function(type, formula, ...) {
  dots <- list(...)

  f <- switch(type,
              "deterministic" = purrr::partial(
                deterministicmodel,
                formula = formula
              ),
              "parametric_distribution" = purrr::partial(
                parametric_distribution_model,
                formula = formula,
                distribution = dots$distribution,
                start = dots$start,
                method = ifelse(is.null(dots$method), "mle", dots$method),
                discrete = ifelse(is.null(dots$discrete), FALSE, dots$discrete),
                fix.arg = dots$fix.arg
              ),
              "linear" = purrr::partial(
                linearmodel,
                formula = formula,
                outcome = dots$outcome,
                boot = dots$boot
              ),
              "exogen" = purrr::partial(
                exogenmodel,
                formula = formula
              ),
              "univariate_fable" = purrr::partial(
                univariate_fable_model,
                formula = formula,
                method = dots$method
              ),
              stop("Unknown model type: ", type)
  )

  class(f) <- c(class(f), type)
  return(f)
}
