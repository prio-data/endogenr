#' A bootstrapped linear model
#'
#' Supports residual based bootstrapping and wild bootstrapping.
#'
#' @param formula A two-sided R formula.
#' @param data A data.frame or data.table.
#' @param type Bootstrap type: `"resid"` or `"wild"`.
#'
#' @return A fitted `lm` object.
#' @keywords internal
bootstraplm <- function(formula, data, type){
  data <- na.omit(data)
  fitted <- stats::lm(formula, data)
  resid <- residuals(fitted)

  if (type == "resid") {
    resampled_residuals <- base::sample(resid, size = length(resid), replace = TRUE)
  } else if (type == "wild") {
    resampled_residuals <- resid * stats::rnorm(length(resid))
  } else {
    stop("Unknown bootstrap type")
  }

  data[[".boot_y"]] <- fitted$fitted.values + resampled_residuals
  stats::lm(stats::update(formula, .boot_y ~ .), data)
}

#' Stage-2 pooled OLS fit on materialized data
#'
#' A file-level function (NOT a closure stored on the model): a closure's
#' environment would drag the whole constructor frame — input data, the
#' materialized copy, and the model itself — into every serialized model
#' object, multiplying the payload shipped to parallel workers by an order
#' of magnitude.
#'
#' @param formula The rewritten, aliased pooled formula.
#' @param data The materialized data.table.
#' @param boot Bootstrap type or `NULL`.
#' @param subset Optional list with `start`/`end` training window.
#' @param timevar Character. Time column name.
#' @return A fitted `lm` object.
#' @keywords internal
.lm_stage2_fit <- function(formula, data, boot, subset, timevar) {
  # Restrict the fit data to the model's own columns (plus timevar for the
  # window filter below), so na.omit in the bootstrap helpers drops only
  # rows missing a model term — matching plain lm()'s estimation sample.
  fit_cols <- unique(c(intersect(all.vars(formula), names(data)), timevar))
  data <- data[, ..fit_cols]
  if (!is.null(subset)) {
    data <- .dt_rows(data, data[[timevar]] >= subset$start &
                             data[[timevar]] <= subset$end)
  }
  if (!is.null(boot)) {
    bootstraplm(formula, data, type = boot)
  } else {
    stats::lm(formula, data)
  }
}



#' @exportS3Method
fit_model.linear_spec <- function(spec, data = NULL, ctx = NULL, subset = NULL, ...) {
  linearmodel(
    formula = spec$formula,
    boot = spec$args$boot,
    data = data, ctx = ctx, subset = subset,
    outcome = spec$args$outcome
  )
}

#' Linear model
#'
#' @param formula A two-sided R formula.
#' @param boot Optional bootstrap type: `"resid"`, `"wild"`, or `NULL`.
#' @param data A data.table or data.frame.
#' @param ctx A panel_context object.
#' @param subset Optional list with start/end for training window.
#' @param ... Additional arguments stored on the model spec.
#'
#' @return An endogenmodel of class `linear`.
#' @keywords internal
linearmodel <- function(formula = NULL, boot = NULL, data = NULL, ctx = NULL,
                        subset = NULL, ...) {
  model <- new_endogenmodel(formula)
  model$boot      <- boot
  model$fit_args  <- rlang::list2(...)
  model$independent <- FALSE

  grp_keys <- ctx_keys(ctx)
  timevar  <- ctx_time(ctx)

  # Stage 1: per-unit ts materialisation. Every maximal time-series
  # sub-expression in the formula is evaluated within group, time-ordered,
  # and stored as a `.pt#` synthetic column. The formula is rewritten to
  # reference those columns so Stage 2 sees only plain variables.
  pm        <- panel_materialize(model$formula, data,
                                 groupvar = grp_keys, timevar = timevar)

  # Build human-readable aliases for the synthetic columns and rename them
  # in-place: `.pt1` -> `lag_x`, `.pt2` -> `lag_log_gdppc`, etc. This keeps
  # coefficient names identical to what base `lm()` would show with a
  # janitor-cleaned column name.
  alias_map <- .pt_make_aliases(pm$map)
  .pt_apply_aliases(pm$data, alias_map)

  # Rewrite the pooled formula to use the alias names.
  fit_formula <- .pt_alias_formula(pm$formula, alias_map)

  model$ts_map      <- pm$map
  model$pt_alias_map <- alias_map
  model$fit_formula  <- fit_formula
  model$data         <- pm$data
  model$timevar      <- timevar
  model$groupvar     <- grp_keys
  model$subset       <- subset

  class(model) <- c("linear", class(model))

  # Stage 2: pooled fit. factor contrasts / poly / spline / interaction bases
  # are resolved across all units here; predict.lm stores predvars/xlevels for
  # coherent basis reconstruction at predict time.
  model$fitted  <- .lm_stage2_fit(fit_formula, pm$data, boot, subset, timevar)
  model$coefs   <- broom::tidy(model$fitted)
  model$gof     <- broom::glance(model$fitted)
  model$outcome <- parse_formula(model)$outcome
  model$required_history <- .required_history(model$formula)

  return(model)
}

#' Predictive draws from a fitted lm prediction object
#'
#' The draw kernel behind [draw_predictive()] for `linear` models. Separates
#' parameter uncertainty from innovation uncertainty: per parameter draw, the
#' residual scale is drawn as `s^2 ~ scale^2 * df / chisq(df)` and the mean is
#' shifted by `se.fit * z * (s / scale)`; innovations are `N(0, s)` per row.
#' The per-row marginal is `fit + sqrt(se.fit^2 + scale^2) * t_df` — identical
#' in distribution to the historical fused single-t draw (a normal /
#' sqrt(chisq/df) mixture).
#'
#' `param_scope` controls the parameter-deviate granularity. In the
#' row-expansion path (`"row"`) each row is one sim instance, so draw an
#' INDEPENDENT deviate per row — a single shared draw would freeze the
#' parameter-uncertainty component across every unit and inner sim at a time
#' step, under-dispersing the ensemble. Under `"draw"` each column block is
#' one coefficient realisation: one deviate shared across rows — a rank-1
#' approximation of the design-based MVN coefficient draw (exact per-row
#' marginals, approximate cross-row dependence).
#'
#' @param lmpred A prediction object from `predict.lm()` with `se.fit = TRUE`.
#' @param n_param Integer >= 0. Number of parameter draws (0 = point estimates).
#' @param n_innov Integer >= 0. Innovation draws per parameter draw (0 = mean).
#' @param param_scope `"draw"` (shared per column block) or `"row"`.
#'
#' @return Numeric matrix, `length(lmpred$fit)` x `max(n_param,1) * max(n_innov,1)`.
#' @keywords internal
.lm_predictive_draws <- function(lmpred, n_param = 1L, n_innov = 1L, param_scope = "draw") {
  fit   <- lmpred$fit
  se    <- lmpred$se.fit
  scale <- lmpred$residual.scale
  df    <- lmpred$df
  n <- length(fit)
  P <- max(n_param, 1L)
  K <- max(n_innov, 1L)
  out <- matrix(NA_real_, n, P * K)
  for (j in seq_len(P)) {
    if (n_param == 0L) {
      mu <- fit
      s  <- rep(scale, n)
    } else if (param_scope == "row") {
      s  <- scale * sqrt(df / stats::rchisq(n, df))
      z  <- stats::rnorm(n)
      mu <- fit + se * z * (s / scale)
    } else {
      s  <- rep(scale * sqrt(df / stats::rchisq(1L, df)), n)
      z  <- rep(stats::rnorm(1L), n)
      mu <- fit + se * z * (s / scale)
    }
    for (k in seq_len(K)) {
      out[, (j - 1L) * K + k] <- if (n_innov == 0L) mu else mu + stats::rnorm(n, 0, s)
    }
  }
  out
}

#' @rdname draw_predictive
#' @export
draw_predictive.linear <- function(model, newdata, n_param = 1L, n_innov = 1L,
                                   param_scope = c("draw", "row"), ...) {
  param_scope <- match.arg(param_scope)
  .check_draw_counts(n_param, n_innov)
  pred <- predict(model$fitted, newdata = newdata, se.fit = TRUE)
  .lm_predictive_draws(pred, n_param, n_innov, param_scope)
}

#' Predict function for a linear model
#'
#' @param model A `linear` endogenmodel.
#' @param data A data.table.
#' @param t Time step to predict.
#' @param ctx A panel_context object.
#' @param what Either "pi" or "expectation".
#' @param ... Ignored, accepted for S3 generic consistency.
#'
#' @return A data.table with key + index + outcome columns.
#' @family simulation
#' @export
predict.linear <- function(model, data, t, ctx, what = "pi", ...) {
  idx      <- ctx_time(ctx)
  all_keys <- ctx_keys(ctx)

  # Subset to the per-unit history window the RHS actually needs.
  data <- .history_subset(data, idx, t, model$required_history)

  # Stage 1 (predict path): re-materialise the per-unit ts columns from the
  # stored ts_map, then rename `.pt#` to readable aliases so newdata column
  # names match those in the fitted model.
  env <- rlang::f_env(model$formula)
  mat <- .apply_ts_map(model$ts_map, data, all_keys, idx, env = env, copy = FALSE)
  .pt_apply_aliases(mat, model$pt_alias_map)

  # Filter to the prediction time step.
  mat <- .dt_rows(mat, mat[[idx]] == t)

  # Stage 2: predict.lm (inside draw_predictive.linear) reproduces the pooled
  # design (factor contrasts, poly/spline bases, interactions) via its stored
  # predvars/xlevels.
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
