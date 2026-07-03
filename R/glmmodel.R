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

  data[[".boot_y"]] <- family$linkinv(eta + resampled_residuals)
  refit <- stats::update(formula, .boot_y ~ .)
  withCallingHandlers(
    stats::glm(refit, data, family = family),
    warning = function(w) {
      if (grepl("non-integer", conditionMessage(w))) invokeRestart("muffleWarning")
    }
  )
}

#' Stage-2 pooled GLM fit on materialized data
#'
#' File-level analogue of `.lm_stage2_fit()` — see that helper for why the
#' fit step must not be a closure stored on the model object.
#'
#' @param formula The rewritten, aliased pooled formula.
#' @param data The materialized data.table.
#' @param family A `stats::family` object.
#' @param boot Bootstrap type or `NULL`.
#' @param subset Optional list with `start`/`end` training window.
#' @param timevar Character. Time column name.
#' @return A fitted `glm` object.
#' @keywords internal
.glm_stage2_fit <- function(formula, data, family, boot, subset, timevar) {
  # Restrict the fit data to the model's own columns (plus timevar for the
  # window filter below), so na.omit in the bootstrap helpers drops only
  # rows missing a model term — matching plain glm()'s estimation sample.
  fit_cols <- unique(c(intersect(all.vars(formula), names(data)), timevar))
  data <- data[, ..fit_cols]
  if (!is.null(subset)) {
    data <- .dt_rows(data, data[[timevar]] >= subset$start &
                             data[[timevar]] <= subset$end)
  }
  if (!is.null(boot)) {
    bootstrapglm(formula, data, family = family, type = boot)
  } else {
    stats::glm(formula, data, family = family)
  }
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
glmmodel <- function(formula = NULL, family = stats::gaussian(), boot = NULL,
                     data = NULL, ctx = NULL, subset = NULL, ...) {
  model <- new_endogenmodel(formula)
  model$boot      <- boot
  model$family    <- family
  model$fit_args  <- rlang::list2(...)
  model$independent <- FALSE

  grp_keys <- ctx_keys(ctx)
  timevar  <- ctx_time(ctx)

  # Stage 1: per-unit ts materialisation (same two-stage approach as linearmodel).
  pm        <- panel_materialize(model$formula, data,
                                 groupvar = grp_keys, timevar = timevar)
  alias_map <- .pt_make_aliases(pm$map)
  .pt_apply_aliases(pm$data, alias_map)
  fit_formula <- .pt_alias_formula(pm$formula, alias_map)

  model$ts_map       <- pm$map
  model$pt_alias_map <- alias_map
  model$fit_formula  <- fit_formula
  model$data         <- pm$data
  model$timevar      <- timevar
  model$groupvar     <- grp_keys
  model$subset       <- subset

  class(model) <- c("glm_endogenr", class(model))

  # Stage 2: pooled GLM fit (file-level helper — see .lm_stage2_fit for why
  # this must not be a closure stored on the model).
  model$fitted <- .glm_stage2_fit(fit_formula, pm$data, family, boot,
                                  subset, timevar)

  model$coefs   <- broom::tidy(model$fitted)
  model$gof     <- broom::glance(model$fitted)
  model$outcome <- parse_formula(model)$outcome
  model$required_history <- .required_history(model$formula)
  # Response-scale dispersion (estimated for gaussian/Gamma/quasi-* families;
  # fixed at 1 for poisson/binomial). Cached so .glm_predictive_draws() need
  # not call the expensive summary() on every predict step.
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

#' Predictive draws from a fitted glm prediction object
#'
#' The draw kernel behind [draw_predictive()] for `glm` models. Parameter
#' uncertainty enters on the link scale (a t-distributed draw around the
#' linear predictor using `se.fit` — lm parity, deliberately `rt` not
#' `rnorm`); the response is then sampled from the family's distribution at
#' the resulting mean, using the estimated `dispersion`. This restores `lm`
#' parity for gaussian GLMs (link draw + residual scale) and yields realistic
#' counts/positive/proportion draws for the other supported families.
#'
#' `param_scope` controls the parameter-deviate granularity. In the
#' row-expansion path (`"row"`) each row is one sim instance, so draw an
#' INDEPENDENT t per row — matching the `linear` kernel and the per-row
#' normal draw in `predict.glmmTMB_endogenr`. (A single shared draw would
#' freeze the parameter-uncertainty component across every unit and inner sim
#' at a time step, under-dispersing the ensemble.) Under `"draw"` each column
#' block keeps one t draw shared across rows: each column is one coefficient
#' realisation.
#'
#' @param glmpred prediction object from predict.glm with se.fit = TRUE and type = "link"
#' @param family a family object
#' @param df residual degrees of freedom
#' @param dispersion Estimated response-scale dispersion (1 for poisson/binomial).
#' @param n_param Integer >= 0. Number of parameter draws (0 = point estimates).
#' @param n_innov Integer >= 0. Innovation draws per parameter draw (0 = mean).
#' @param param_scope `"draw"` (shared per column block) or `"row"`.
#'
#' @return A numeric matrix of samples on the response scale,
#'   `length(glmpred$fit)` x `max(n_param,1) * max(n_innov,1)`.
#' @keywords internal
.glm_predictive_draws <- function(glmpred, family, df, dispersion = 1,
                                  n_param = 1L, n_innov = 1L, param_scope = "draw") {
  n <- length(glmpred$fit)
  P <- max(n_param, 1L)
  K <- max(n_innov, 1L)
  eta <- if (n_param == 0L) {
    matrix(glmpred$fit, n, P)
  } else if (param_scope == "row") {
    glmpred$fit + matrix(stats::rt(n * P, df), n, P) * glmpred$se.fit
  } else {
    glmpred$fit + outer(glmpred$se.fit, stats::rt(P, df))
  }
  mu <- family$linkinv(eta)
  out <- matrix(NA_real_, n, P * K)
  for (j in seq_len(P)) {
    for (k in seq_len(K)) {
      out[, (j - 1L) * K + k] <-
        if (n_innov == 0L) mu[, j]
        else .glm_response_draw(as.vector(mu[, j]), family$family, dispersion)
    }
  }
  out
}

#' @rdname draw_predictive
#' @export
draw_predictive.glm_endogenr <- function(model, newdata, n_param = 1L, n_innov = 1L,
                                         param_scope = c("draw", "row"), ...) {
  param_scope <- match.arg(param_scope)
  .check_draw_counts(n_param, n_innov)
  pred <- predict(model$fitted, newdata = newdata, type = "link", se.fit = TRUE)
  .glm_predictive_draws(pred, model$family, model$fitted$df.residual,
                        dispersion = model$dispersion,
                        n_param = n_param, n_innov = n_innov, param_scope = param_scope)
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
  idx      <- ctx_time(ctx)
  all_keys <- ctx_keys(ctx)

  # Subset to the per-unit history window the RHS actually needs.
  data <- .history_subset(data, idx, t, model$required_history)

  # Re-materialise ts columns per (unit, sim) group, then apply aliases.
  env <- rlang::f_env(model$formula)
  mat <- .apply_ts_map(model$ts_map, data, all_keys, idx, env = env, copy = FALSE)
  .pt_apply_aliases(mat, model$pt_alias_map)

  # Filter to the prediction time step.
  mat <- .dt_rows(mat, mat[[idx]] == t)

  # Link-scale prediction and draws happen inside draw_predictive.glm_endogenr.
  result_cols <- c(all_keys, idx, model$outcome)
  result      <- mat[, ..result_cols]

  if (what == "expectation") {
    data.table::set(result, j = model$outcome,
                    value = as.vector(draw_predictive(model, mat, n_param = 0L, n_innov = 0L)))
  } else if (what == "pi") {
    data.table::set(result, j = model$outcome,
                    value = as.vector(draw_predictive(model, mat, n_param = 1L, n_innov = 1L,
                                                      param_scope = "row")))
  } else {
    stop("`what` must be either `pi` or `expectation`")
  }

  return(result)
}
