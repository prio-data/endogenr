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

  outcome <- fitted$terms |> rlang::f_lhs() |> as.character()
  data[[outcome]] <- fitted$fitted.values + resampled_residuals
  stats::lm(formula, data)
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
  old <- intersect(names(alias_map), names(pm$data))
  if (length(old) > 0L) data.table::setnames(pm$data, old, alias_map[old])

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
  model$fit <- function(formula, data, boot, subset, timevar) {
    if (!is.null(subset)) {
      data <- data[data[[timevar]] >= subset$start & data[[timevar]] <= subset$end]
    }
    if (!is.null(boot)) {
      bootstraplm(formula, data, type = boot)
    } else {
      stats::lm(formula, data)
    }
  }

  model$fitted  <- model$fit(model$fit_formula, model$data, model$boot,
                             model$subset, model$timevar)
  model$coefs   <- broom::tidy(model$fitted)
  model$gof     <- broom::glance(model$fitted)
  model$outcome <- parse_formula(model)$outcome
  model$required_history <- .required_history(model$formula)

  return(model)
}

#' Get the standard error for prediction
#'
#' @param lmpred A prediction object from `predict.lm()` with `se.fit = TRUE`.
#'
#' @return Numeric vector. Per-row standard error including residual scale.
#' @keywords internal
get_sepi <- function(lmpred){
  se <- lmpred$se.fit
  scale <- lmpred$residual.scale
  sqrt(se^2 + scale^2)
}

#' Selects a column per row in a matrix
#'
#' @param mat A numeric matrix.
#' @param column_ids Integer vector. One column index per row of `mat`.
#'
#' @return A numeric vector with one value per row.
#' @keywords internal
select_col_per_row <- function(mat, column_ids){
  cidx <- cbind(1:nrow(mat), column_ids)
  mat[cidx]
}

#' Get the predictive distribution from a linear model
#'
#' Samples nsamples points from the predictive distribution.
#'
#' @param lmpred A prediction object from `predict.lm()` with `se.fit = TRUE`.
#' @param nsamples Integer. Number of draws (1 for row-expansion path).
#'
#' @return Numeric vector (or matrix if `nsamples > 1`) of predictive draws.
#' @keywords internal
getpi <- function(lmpred, nsamples = 1){
  sepi <- get_sepi(lmpred)
  if (nsamples == 1) {
    # Row-expansion architecture: each row is already one sim instance, so draw
    # one independent t-value per row (length(sepi) draws).
    as.vector(lmpred$fit + sepi * stats::rt(length(sepi), lmpred$df))
  } else {
    # Multi-sample path (e.g. longhorizon): one row per unit, nsamples draws
    # each — outer() produces an n x nsamples matrix.
    lmpred$fit + outer(sepi, stats::rt(nsamples, lmpred$df))
  }
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
  old <- intersect(names(model$pt_alias_map), names(mat))
  if (length(old) > 0L) data.table::setnames(mat, old, model$pt_alias_map[old])

  # Filter to the prediction time step.
  mat <- mat[mat[[idx]] == t]

  # Stage 2: predict.lm reproduces the pooled design (factor contrasts,
  # poly/spline bases, interactions) via its stored predvars/xlevels.
  pred <- predict(model$fitted, newdata = mat, se.fit = TRUE)

  result_cols <- c(all_keys, idx, model$outcome)
  result      <- mat[, ..result_cols]

  if (what == "expectation") {
    data.table::set(result, j = model$outcome, value = pred$fit)
  } else if (what == "pi") {
    data.table::set(result, j = model$outcome, value = getpi(pred))
  } else {
    stop("`what` must be either `pi` or `expectation`")
  }

  return(result)
}
