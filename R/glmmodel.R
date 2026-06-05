#' A bootstrapped GLM
#'
#' Supports residual based bootstrapping and wild bootstrapping on the link scale.
#'
#' @param formula A two-sided R formula.
#' @param data A data.frame or data.table.
#' @param family A `family` object (e.g. `stats::gaussian()`).
#' @param type Bootstrap type: `"resid"` or `"wild"`.
#'
#' @return A fitted `glm` object.
#' @keywords internal
bootstrapglm <- function(formula, data, family, type){
  data <- na.omit(data)
  fitted <- stats::glm(formula, data, family = family)

  # Work on the link scale for resampling
  eta <- predict(fitted, type = "link")
  resid_working <- residuals(fitted, type = "working")

  if (type == "resid") {
    resampled_residuals <- base::sample(resid_working, size = length(resid_working), replace = TRUE)
  } else if (type == "wild") {
    resampled_residuals <- resid_working * stats::rnorm(length(resid_working))
  } else {
    stop("Unknown bootstrap type")
  }

  outcome <- fitted$terms |> rlang::f_lhs() |> as.character()
  data[[outcome]] <- family$linkinv(eta + resampled_residuals)
  stats::glm(formula, data, family = family)
}


#' @exportS3Method
fit_model.glm_spec <- function(spec, data = NULL, ctx = NULL, subset = NULL, ...) {
  family <- if (!is.null(spec$args$family)) spec$args$family else stats::gaussian()
  glmmodel(
    formula = spec$formula,
    family = family,
    boot = spec$args$boot,
    data = data, ctx = ctx, subset = subset
  )
}

#' GLM model
#'
#' @param formula A two-sided R formula.
#' @param family A family object (e.g. quasibinomial(), gaussian(), poisson())
#' @param boot Optional bootstrap type: `"resid"`, `"wild"`, or `NULL`.
#' @param data A data.table or data.frame.
#' @param ctx A panel_context object.
#' @param subset Optional list with start/end for training window.
#' @param ... Additional arguments stored on the model spec.
#'
#' @return An endogenmodel of class `glm_endogenr`.
#' @keywords internal
glmmodel <- function(formula = NULL, family = stats::gaussian(), boot = NULL, data = NULL, ctx = NULL, subset = NULL, ...){
  model <- new_endogenmodel(formula)
  model$boot <- boot
  model$family <- family
  model$fit_args <- rlang::list2(...)
  model$independent <- FALSE

  panel_frame <- create_panel_frame(model$formula, data, ctx)
  model$naive_formula <- panel_frame$naive_formula
  model$data <- panel_frame$data
  model$timevar <- ctx_time(ctx)
  model$subset <- subset

  # Cache materialization helpers for fast predict (reuse from create_panel_frame)
  model$mat_formula <- panel_frame$mat_formula
  model$col_mapping <- panel_frame$col_mapping

  class(model) <- c("glm_endogenr", class(model))

  model$fit <- function(formula, data, family, boot, subset, timevar){
    if(!is.null(subset)){
      data <- data[data[[timevar]] >= subset$start & data[[timevar]] <= subset$end]
    }

    if(!is.null(boot)){
      fitted <- bootstrapglm(formula, data, family = family, type = boot)
    } else {
      fitted <- stats::glm(formula, data, family = family)
    }
  }

  model$fitted <- model$fit(model$naive_formula, model$data, model$family, model$boot, model$subset, model$timevar)

  model$coefs <- broom::tidy(model$fitted)
  model$gof <- broom::glance(model$fitted)

  model$outcome <- parse_formula(model)$outcome
  model$required_history <- .required_history(model$formula)
  # Response-scale dispersion (estimated for gaussian/Gamma/quasi-* families;
  # fixed at 1 for poisson/binomial). Cached so getpi_glm() need not call the
  # expensive summary() on every predict step.
  model$dispersion <- summary(model$fitted)$dispersion

  return(model)
}

# Family-specific response-scale draw given the fitted mean(s) `mu` and the
# estimated `dispersion`. Returns a numeric vector the same length as `mu`.
# `family` is the glm family name string. Unsupported families return `mu`
# unchanged (parameter-uncertainty only) and warn once per session.
.glm_response_draw <- function(mu, family, dispersion) {
  n <- length(mu)
  switch(family,
    "gaussian" = stats::rnorm(n, mean = mu, sd = sqrt(max(dispersion, 0))),
    "poisson"  = stats::rpois(n, lambda = pmax(mu, 0)),
    "quasipoisson"  = .glm_draw_qpois(mu, dispersion),
    "Gamma"         = .glm_draw_gamma(mu, dispersion),
    "binomial"      = .glm_draw_prop(mu, dispersion),
    "quasibinomial" = .glm_draw_prop(mu, dispersion),
    {
      .glm_warn_unsupported(family)
      mu
    }
  )
}

# Quasi-Poisson: overdispersed counts via a negative-binomial approximation
# with mean `mu` and variance `dispersion * mu` (size = mu / (dispersion - 1)).
# Falls back to Poisson when there is no overdispersion.
.glm_draw_qpois <- function(mu, dispersion) {
  mu <- pmax(mu, 0)
  if (!is.finite(dispersion) || dispersion <= 1) {
    return(stats::rpois(length(mu), mu))
  }
  out <- numeric(length(mu))
  pos <- mu > 0
  if (any(pos)) {
    out[pos] <- stats::rnbinom(sum(pos), size = mu[pos] / (dispersion - 1),
                               mu = mu[pos])
  }
  out
}

# Gamma: mean `mu`, variance `dispersion * mu^2` via shape = 1/dispersion,
# scale = mu * dispersion.
.glm_draw_gamma <- function(mu, dispersion) {
  mu <- pmax(mu, .Machine$double.eps)
  if (!is.finite(dispersion) || dispersion <= 0) return(mu)
  stats::rgamma(length(mu), shape = 1 / dispersion, scale = mu * dispersion)
}

# Binomial / quasibinomial PROPORTION outcomes (no trial count in the grid):
# draw from a Beta with mean `mu` and precision derived from `dispersion`. The
# Beta variance is capped at the Bernoulli value mu(1 - mu); overdispersion
# beyond that (quasibinomial dispersion > 1) cannot be represented without a
# trial count, so it collapses to ~Bernoulli draws. See the "Known issues"
# section of ?endogenr for the proportion assumption this encodes.
.glm_draw_prop <- function(mu, dispersion) {
  mu <- pmin(pmax(mu, 0), 1)
  phi <- if (is.finite(dispersion) && dispersion > 0) {
    max(1 / dispersion - 1, 1e-6)
  } else 1e-6
  a <- pmax(mu * phi, 1e-9)
  b <- pmax((1 - mu) * phi, 1e-9)
  out <- stats::rbeta(length(mu), shape1 = a, shape2 = b)
  out[mu <= 0] <- 0
  out[mu >= 1] <- 1
  out
}

# One-time-per-session warning for families with no response-scale sampler.
.glm_warn_unsupported <- function(family) {
  key <- paste0(".endogenr_glm_warned_", family)
  if (!isTRUE(getOption(key))) {
    warning("GLM family '", family, "' has no response-scale predictive draw; ",
            "falling back to a link-scale (parameter-uncertainty-only) draw, ",
            "which under-disperses. ", call. = FALSE)
    do.call(options, stats::setNames(list(TRUE), key))
  }
  invisible(NULL)
}

#' Get the predictive distribution from a GLM
#'
#' Samples `nsamples` points from the predictive distribution on the response
#' scale. Parameter uncertainty enters on the link scale (a t-distributed draw
#' around the linear predictor using `se.fit`); the response is then sampled
#' from the family's distribution at the resulting mean, using the estimated
#' `dispersion`. This restores `lm` parity for gaussian GLMs (link draw +
#' residual scale) and yields realistic counts/positive/proportion draws for
#' the other supported families.
#'
#' @param glmpred prediction object from predict.glm with se.fit = TRUE and type = "link"
#' @param family a family object
#' @param df residual degrees of freedom
#' @param dispersion Estimated response-scale dispersion (1 for poisson/binomial).
#' @param nsamples Integer. Number of draws from the predictive distribution.
#'
#' @return A numeric vector (or matrix if `nsamples > 1`) of samples on the response scale.
#' @keywords internal
getpi_glm <- function(glmpred, family, df, dispersion = 1, nsamples = 1){
  # Parameter uncertainty: draw the linear predictor on the link scale, then
  # map to the conditional mean.
  eta <- glmpred$fit + outer(glmpred$se.fit, stats::rt(nsamples, df))
  mu  <- family$linkinv(eta)
  # Response-scale draw at each mean (adds the family's predictive dispersion).
  draw <- .glm_response_draw(as.vector(mu), family$family, dispersion)
  if (nsamples == 1) draw else matrix(draw, nrow = nrow(mu))
}

#' Predict function for a GLM model
#'
#' @param model A `glm_endogenr` endogenmodel.
#' @param data A data.table.
#' @param t Time step to predict.
#' @param ctx A panel_context object.
#' @param what Either "pi" or "expectation".
#' @param ... Ignored, accepted for S3 generic consistency.
#'
#' @return A data.table with key + index + outcome columns.
#' @family simulation
#' @export
predict.glm_endogenr <- function(model, data, t, ctx, what = "pi", ...) {
  idx <- ctx_time(ctx)
  grp <- ctx_unit(ctx)
  all_keys <- ctx_keys(ctx)

  # Subset to the per-unit history window the RHS actually needs
  data <- .history_subset(data, idx, t, model$required_history)

  # Materialize formula (use cached helpers to skip update/inject/clean_names)
  data <- materialize_formula(model$formula, data, ctx,
                              .mat_formula = model$mat_formula,
                              .col_mapping = model$col_mapping)

  # Filter for specific time point
  data <- data[data[[idx]] == t]

  # Predict on link scale with standard errors
  pred <- predict(model$fitted, newdata = data, type = "link", se.fit = TRUE)

  # Build result data.table
  result_cols <- c(all_keys, idx, model$outcome)
  result <- data[, ..result_cols]

  if (what == "expectation") {
    data.table::set(result, j = model$outcome, value = model$family$linkinv(pred$fit))
  } else if (what == "pi") {
    data.table::set(result, j = model$outcome,
                    value = getpi_glm(pred, model$family, model$fitted$df.residual,
                                      dispersion = model$dispersion))
  } else {
    stop("`what` must be either `pi` or `expectation`")
  }

  return(result)
}
