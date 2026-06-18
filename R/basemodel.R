#' A baseline structure for a model
#'
#' Do not use directly.
#'
#' @param formula A two-sided R formula stored on the model spec.
#'
#' @return A list with class endogenmodel
#' @keywords internal
new_endogenmodel <- function(formula){
  structure(
    list(
      formula = formula
    ),
    class = "endogenmodel"
  )
}

#' Build a model specification
#'
#' Creates a model specification that will become part of a simulation system.
#' The spec is a lightweight object storing the model type, formula, and
#' arguments. Actual fitting happens later via [fit_model()] (called
#' internally by [fit_system()]).
#'
#' @section Model types and required arguments:
#'
#' \describe{
#'   \item{`"deterministic"`}{Two-sided formula `outcome ~ I(expr)`. The RHS
#'     must be wrapped in `I()`. Evaluated at each simulated time step `t`.
#'     No fitting; no extra arguments.}
#'   \item{`"parametric_distribution"`}{One-sided formula `~ var`. Pass
#'     `distribution = "norm"` (or any distribution name accepted by
#'     [fitdistrplus::fitdist()]). Extra arguments (`start`, `method`,
#'     `lower`, `upper`, …) are forwarded to `fitdist()`.}
#'   \item{`"linear"`}{Two-sided formula. Optional `boot ∈ {"resid","wild"}`
#'     selects residual or wild bootstrap; omit `boot` for plain OLS.}
#'   \item{`"glm"`}{Two-sided formula. `family = stats::gaussian()` by
#'     default; pass any `stats::family` (e.g. `stats::quasibinomial()`).
#'     Optional `boot` as for `"linear"`.}
#'   \item{`"exogen"`}{One-sided formula `~var`. The variable must already
#'     be present in `data` for every row of the forecast horizon — the
#'     model just copies those values into the simulation grid.}
#'   \item{`"univariate_fable"`}{Two-sided fable formula (e.g.
#'     `y ~ error("A") + trend("N") + season("N")`). Pass `method = "ets"`
#'     or `method = "arima"`. Requires the `fable`/`fabletools`/`tsibble`
#'     packages.}
#'   \item{`"heterolm"`}{Two-sided mean formula plus `variance = ~ ...`
#'     (one-sided log-variance formula, defaults to `~ 1`). Requires the
#'     `heterolm` package.}
#'   \item{`"spatial_lag"`}{Two-sided formula `sl_y ~ lag(y)` (use `lag()`
#'     to avoid a circular dependency on a same-period outcome). Pass `nb`,
#'     `wt`, and `unit_ids` from [st_weights_from_sf()] or `sfdep`
#'     directly. Optional `island_default` for units with no neighbours.}
#'   \item{`"glmmTMB"`}{Two-sided formula, including lme4-style random-effects
#'     bars `(1 + lag(x) | group)` and glmmTMB covariance-structure wrappers
#'     (e.g. `ar1(times + 0 | group)`, `us(…|g)`). Pass `family =` (default
#'     `stats::gaussian()`), `dispformula = ~…` (default `~1`), `ziformula =
#'     ~…` (default `~0`), and optionally `control =
#'     glmmTMB::glmmTMBControl()`. Grouping factors and cov-struct coordinate
#'     columns (e.g. `region`, `times`) are read at every forecast step, so —
#'     like any predictor — they must be produced by some model: add an
#'     `exogen` (e.g. `build_model("exogen", formula = ~region)`) to carry the
#'     column forward, or group by a panel key (`unit`/`time`), which is always
#'     present. Temporal covariance structures (`ar1`/`ou`/…) are forecast
#'     multi-step by predicting the whole forecast-so-far block at each step so
#'     glmmTMB applies the correct `phi^k` decay; the coordinate must be carried
#'     into the horizon and be contiguous and unit-spaced. Response-scale
#'     predictive draws are implemented for `gaussian`, `poisson`, `binomial`,
#'     `Gamma`, `nbinom1`, `nbinom2`, `beta`, `betabinomial`, `t`, `lognormal`,
#'     `skewnormal`, and
#'     `truncated_poisson`/`truncated_nbinom1`/`truncated_nbinom2`; `tweedie`
#'     needs the `tweedie` package; other families fall back to the conditional
#'     mean. Requires the `glmmTMB` package.}
#'   \item{`"gamlss"`}{Two-sided `formula` for the location parameter `mu`
#'     (may include `pb()`/`cs()`/`lo()` smoothers and `random()`/`ra()`/`re()`
#'     grouping terms). Optional `sigma.formula`, `nu.formula`, `tau.formula`
#'     (one-sided, default `~1`). `family =` a `gamlss.family` object (default
#'     `gamlss.dist::NO()`). Optional `control = gamlss::gamlss.control(...)`.
#'     Grouping factors inside `random()`/`ra()`/`re()` are read at every
#'     forecast step, so they must be produced by some model — add an `exogen`
#'     (e.g. `build_model("exogen", formula = ~region)`) to carry the grouping
#'     column forward, or group by a panel key. Requires the `gamlss` package.}
#' }
#'
#' @param type One of `"deterministic"`, `"parametric_distribution"`,
#'   `"linear"`, `"glm"`, `"exogen"`, `"univariate_fable"`, `"heterolm"`,
#'   `"spatial_lag"`, `"glmmTMB"`, or `"gamlss"`.
#' @param formula An R formula. See the model-type section for the expected
#'   shape per type.
#' @param ... Model-specific arguments. See the model-type section.
#'
#' @return An `endogenr_spec` object (a list with `$type`, `$formula`, `$args`).
#' @seealso [setup_system()], [fit_system()], [simulate_system()], [fit_model()]
#' @family build
#' @export
#'
#' @examples
#' df <- endogenr::example_data
#' train <- df[df$year >= 1970 & df$year < 2010, ]
#' c1 <- yjbest ~ lag(zoo::rollsumr(yjbest, k = 5, fill = NA)) + lag(log(gdppc))
#' model_system <- list(
#'   build_model("deterministic", formula = gdppc ~ I(abs(lag(gdppc)*(1+gdppc_grwt)))),
#'   build_model("deterministic", formula = gdp ~ I(abs(gdppc*population))),
#'   build_model("parametric_distribution", formula = ~gdppc_grwt, distribution = "t_ls",
#'     start = list(df = 1, mu = mean(train$gdppc_grwt), sigma = sd(train$gdppc_grwt))),
#'   build_model("linear", formula = c1, boot = "resid"),
#'   build_model("exogen", formula = ~psecprop),
#'   build_model("exogen", formula = ~population)
#' )
build_model <- function(type, formula, ...) {
  valid_types <- c("deterministic", "parametric_distribution", "linear", "glm",
                   "exogen", "univariate_fable", "heterolm", "spatial_lag",
                   "glmmTMB", "gamlss")
  if (!type %in% valid_types) {
    stop("Unknown model type: ", type)
  }

  dots <- list(...)

  spec <- structure(
    list(type = type, formula = formula, args = dots),
    class = c(paste0(type, "_spec"), "endogenr_spec")
  )
  spec
}


#' Fit a model from a specification
#'
#' Generic function that dispatches to type-specific fitting methods based on
#' Generic function that dispatches to type-specific fitting methods based on
#' the spec's class. Called internally by [fit_system()].
#'
#' @param spec An `endogenr_spec` object from [build_model()].
#' @param ... Arguments passed to the type-specific method (typically `data`,
#'   `ctx`, `subset`).
#'
#' @return A fitted endogenmodel object.
#' @family build
#' @export
fit_model <- function(spec, ...) UseMethod("fit_model")
