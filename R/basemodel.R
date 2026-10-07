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
#'   \item{`"cross_section"`}{Two-sided formula `outcome ~ I(expr)` with a
#'     single `I()` term. `expr` is evaluated across all units in the
#'     simulation at each simulated time step `t`, separately within each
#'     simulation draw. A scalar result (e.g.
#'     `quantile(x, 0.9, na.rm = TRUE)`, `mean(x)`) is copied to every unit; a
#'     result with one value per unit (e.g. `rank(x)`, `x - mean(x)`) is kept
#'     per unit. Time-series functions (`lag`, `diff`, rolling, cumulative) are
#'     not allowed; reference a lagged column produced by another model (e.g.
#'     `build_model("deterministic", x_l1 ~ I(lag(x)))`). Training-period
#'     values of the outcome column must already be in `data` (e.g.
#'     `data[, front := quantile(x, 0.9, na.rm = TRUE), by = year]`). No
#'     fitting; no extra arguments.}
#'   \item{`"parametric_distribution"`}{One-sided formula `~ var`. Pass
#'     `distribution = "norm"` (or any distribution name accepted by
#'     [fitdistrplus::fitdist()]). Extra arguments (`start`, `method`,
#'     `lower`, `upper`, …) are forwarded to `fitdist()`. `start` may be a
#'     named list **or a `function(x)`** returning one — the function is
#'     evaluated on the training-window data at fit time, so it stays correct
#'     under [run_experiments()] / sliding-window refits. For the bundled
#'     location-scale Student-t (`distribution = "t_ls"`, see [t_ls]) starting
#'     values are derived automatically (median/MAD location-scale,
#'     kurtosis-matched `df`) when none are given, and the `t_ls` d/p/q/r
#'     functions are exported and made visible to `fitdist()` automatically
#'     even when endogenr is not attached. Optional
#'     `param_uncertainty = TRUE` propagates MLE parameter uncertainty by
#'     drawing one parameter vector per simulation from the asymptotic
#'     MVN(estimate, vcov); default `FALSE` conditions on the point
#'     estimates.}
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
#' @param type One of `"deterministic"`, `"cross_section"`,
#'   `"parametric_distribution"`, `"linear"`, `"glm"`, `"exogen"`,
#'   `"univariate_fable"`, `"heterolm"`, `"spatial_lag"`, `"glmmTMB"`, or
#'   `"gamlss"`.
#' @param formula An R formula. See the model-type section for the expected
#'   shape per type.
#' @param ... Model-specific arguments. See the model-type section.
#' @param bounds Optional length-2 numeric `c(lower, upper)`. When set, this
#'   model's simulated outcome is clamped to `[lower, upper]` at every forecast
#'   step of [simulate_system()], before the value feeds later steps. Use it to
#'   stabilize autoregressive feedback (e.g. `bounds = c(-1, 1)` for a growth
#'   rate). Infinities are allowed for a one-sided limit (`c(0, Inf)`). A draw
#'   that returns non-finite after clamping is reset to a finite in-range value
#'   (the midpoint of two finite bounds, or the finite bound of a one-sided
#'   limit), so a bounded outcome never propagates `NaN`/`Inf`. Default `NULL`
#'   applies no clamping.
#'
#' @return An `endogenr_spec` object (a list with `$type`, `$formula`, `$args`,
#'   and `$bounds`).
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
build_model <- function(type, formula, ..., bounds = NULL) {
  valid_types <- c("deterministic", "cross_section", "parametric_distribution",
                   "linear", "glm", "exogen", "univariate_fable", "heterolm",
                   "spatial_lag", "glmmTMB", "gamlss")
  if (!type %in% valid_types) {
    stop("Unknown model type: ", type)
  }

  if (!is.null(bounds)) {
    if (!is.numeric(bounds) || length(bounds) != 2L ||
        anyNA(bounds) || bounds[1L] > bounds[2L]) {
      stop("`bounds` must be a length-2 numeric c(lower, upper) with lower <= upper.",
           call. = FALSE)
    }
  }

  if (identical(type, "cross_section")) .check_cross_section_formula(formula)

  dots <- list(...)

  spec <- structure(
    list(type = type, formula = formula, args = dots, bounds = bounds),
    class = c(paste0(type, "_spec"), "endogenr_spec")
  )
  spec
}


#' Fit a model from a specification
#'
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
