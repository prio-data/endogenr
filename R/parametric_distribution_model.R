#' Student-t distribution density function
#'
#' @param x Numeric vector.
#' @param df Degrees of freedom.
#' @param mu Location parameter.
#' @param sigma Scale parameter.
#'
#' @return Numeric vector of density values.
#' @keywords internal
dt_ls <- function(x, df = 1, mu = 0, sigma = 1) {
  1 / sigma * stats::dt((x - mu) / sigma, df)
}

#' Student-t distribution distribution function
#'
#' @param q Numeric vector of quantiles.
#' @param df Degrees of freedom.
#' @param mu Location parameter.
#' @param sigma Scale parameter.
#'
#' @return Numeric vector of probabilities.
#' @keywords internal
pt_ls <- function(q, df = 1, mu = 0, sigma = 1) {
  stats::pt((q - mu) / sigma, df)
}

#' Student-t distribution quantile function
#'
#' @param p Numeric vector of probabilities.
#' @param df Degrees of freedom.
#' @param mu Location parameter.
#' @param sigma Scale parameter.
#'
#' @return Numeric vector of quantiles.
#' @keywords internal
qt_ls <- function(p, df = 1, mu = 0, sigma = 1) {
  stats::qt(p, df) * sigma + mu
}

#' Student-t distribution random generation
#'
#' @param n Number of observations.
#' @param df Degrees of freedom.
#' @param mu Location parameter.
#' @param sigma Scale parameter.
#'
#' @return Numeric vector of random deviates.
#' @keywords internal
rt_ls <- function(n, df = 1, mu = 0, sigma = 1) {
  stats::rt(n, df) * sigma + mu
}

#' Sample from a fitted distribution
#'
#' Generates random samples from the distribution fitted by fitdistrplus::fitdist().
#'
#' @param fitobj A fitdist object.
#' @param n Integer. Number of samples to draw.
#'
#' @return Numeric vector of random samples.
#' @keywords internal
.sample_from_fitdist <- function(fitobj, n) {
  dname <- fitobj$distname
  est <- as.list(fitobj$estimate)

  if (dname == "t_ls") {
    return(rt_ls(n, df = est$df, mu = est$mu, sigma = est$sigma))
  }
  if (dname == "nbinom") {
    prob <- est$size / (est$size + est$mu)
    return(stats::rnbinom(n, size = est$size, prob = prob))
  }

  # Standard distributions: rnorm, rcauchy, rgamma, rpois, etc.
  rname <- paste0("r", dname)
  rfun <- tryCatch(match.fun(rname), error = function(e) NULL)

  if (is.null(rfun)) {
    # Try actuar for gumbel etc.
    if (requireNamespace("actuar", quietly = TRUE)) {
      rfun <- tryCatch(utils::getFromNamespace(rname, "actuar"), error = function(e) NULL)
    }
  }

  if (is.null(rfun)) {
    stop(sprintf("Cannot find random generation function '%s' for distribution '%s'", rname, dname))
  }

  do.call(rfun, c(list(n = n), est))
}

#' Fits a parametric distribution using fitdistrplus::fitdist
#'
#' @param model An endogenmodel with `$outcome`, `$distribution`, and optional `$fit_args`.
#' @param data A data.frame or data.table containing the outcome column.
#'
#' @return A fitdist object.
#' @keywords internal
fit_parametric_distribution_model <- function(model, data) {
  if (is.null(model$fit_args)) {
    fitted <- fitdistrplus::fitdist(data[[model$outcome]], model$distribution)
  } else {
    args <- model$fit_args
    args$data <- data[[model$outcome]]
    args$distr <- model$distribution
    fitted <- do.call(fitdistrplus::fitdist, args)
  }
  fitted
}

#' Parametric distribution model
#'
#' This model is static across time, and therefore independent/exogenous.
#' Uses fitdistrplus::fitdist() to fit a distribution to the pooled training data.
#'
#' @param formula A one-sided formula (e.g. `~gdppc_grwt`).
#' @param distribution Character. Distribution name compatible with fitdistrplus.
#' @param data A data.table or data.frame.
#' @param ctx A panel_context object.
#' @param ... Additional arguments passed to fitdistrplus::fitdist().
#'
#' @return An endogenmodel of class `parametric_distribution`.
#' @export
#' @exportS3Method
fit_model.parametric_distribution_spec <- function(spec, data = NULL, ctx = NULL, ...) {
  # Extract fitdist-specific args, excluding 'distribution'
  extra_args <- spec$args[!names(spec$args) %in% "distribution"]
  do.call(parametric_distribution_model, c(
    list(formula = spec$formula, distribution = spec$args$distribution,
         data = data, ctx = ctx),
    extra_args
  ))
}

#' @keywords internal
parametric_distribution_model <- function(formula = NULL, distribution = NULL, data = NULL, ctx = NULL, ...) {
  model <- new_endogenmodel(formula)
  model$distribution <- distribution
  model$fit_args <- rlang::list2(...)

  class(model) <- c("parametric_distribution", class(model))
  model$independent <- TRUE
  model$outcome <- parse_formula(model)$outcome
  model$fitted <- fit_parametric_distribution_model(model, data)

  return(model)
}

#' Predict function for parametric_distribution models
#'
#' Generates random samples from the fitted distribution for all forecast rows.
#'
#' @param model A parametric_distribution endogenmodel.
#' @param data A data.table (the simulation grid).
#' @param ctx A panel_context object.
#' @param test_start Integer. Start of the forecast period.
#' @param horizon Integer. Number of forecast steps.
#' @param inner_sims Integer. Number of inner simulations.
#' @param ... Ignored.
#'
#' @return A data.table with columns: unit, sim, time, outcome.
#' @export
predict.parametric_distribution <- function(model, data, ctx, test_start, horizon, inner_sims, ...) {
  unit_col <- ctx_unit(ctx)
  time_col <- ctx_time(ctx)
  all_keys <- ctx_keys(ctx)

  # Filter to forecast period
  pred_data <- data[data[[time_col]] >= test_start & data[[time_col]] <= (test_start + horizon - 1)]

  # Keep only relevant columns
  result_cols <- c(all_keys, time_col, model$outcome)
  result <- pred_data[, ..result_cols]

  # Generate samples from the fitted distribution
  n <- nrow(result)
  data.table::set(result, j = model$outcome, value = .sample_from_fitdist(model$fitted, n))

  return(result)
}
