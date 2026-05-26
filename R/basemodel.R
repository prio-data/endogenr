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
#' internally by [setup_simulator()]).
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
#' }
#'
#' @param type One of "deterministic", "parametric_distribution", "linear",
#'   "glm", "exogen", "univariate_fable", "heterolm", or "spatial_lag".
#' @param formula An R formula. See the model-type section for the expected
#'   shape per type.
#' @param ... Model-specific arguments. See the model-type section.
#'
#' @return An `endogenr_spec` object (a list with `$type`, `$formula`, `$args`).
#' @seealso [setup_simulator()], [simulate_endogenr()], [fit_model()]
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
                   "exogen", "univariate_fable", "heterolm", "spatial_lag")
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
#' the spec's class. Called internally by [setup_simulator()] and
#' [inner_simulation()].
#'
#' @param spec An `endogenr_spec` object from [build_model()].
#' @param ... Arguments passed to the type-specific method (typically `data`,
#'   `ctx`, `subset`).
#'
#' @return A fitted endogenmodel object.
#' @family build
#' @export
fit_model <- function(spec, ...) UseMethod("fit_model")
