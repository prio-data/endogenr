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

#' Build a model specification
#'
#' Creates a model specification that will become part of a simulation system.
#' The spec is a lightweight object storing the model type, formula, and
#' arguments. Actual fitting happens later via [fit_model()] (called
#' internally by [setup_simulator()]).
#'
#' @param type One of "deterministic", "parametric_distribution", "linear",
#'   "glm", "exogen", "univariate_fable", "heterolm", or "spatial_lag".
#' @param formula An R formula. See details for each model type.
#' @param ... Model-specific arguments (e.g. `boot`, `distribution`, `family`,
#'   `variance`, `nb`, `wt`, `unit_ids`).
#'
#' @return An `endogenr_spec` object (a list with `$type`, `$formula`, `$args`).
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
#' @export
fit_model <- function(spec, ...) UseMethod("fit_model")
