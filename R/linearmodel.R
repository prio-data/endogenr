#' A bootstrapped linear model
#'
#' Supports residual based bootstrapping and wild bootstrapping.
#'
#' @param formula
#' @param data
#' @param type
#'
#' @return
#' @keywords internal
#'
#' @examples
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
#' @param formula
#' @param boot
#' @param data A data.table or data.frame.
#' @param ctx A panel_context object.
#' @param subset Optional list with start/end for training window.
#' @param ...
#'
#' @return
#' @keywords internal
#'
#' @examples
linearmodel <- function(formula = NULL, boot = NULL, data = NULL, ctx = NULL, subset = NULL, ...){
  model <- new_endogenmodel(formula)
  model$boot <- boot
  model$fit_args <- rlang::list2(...)
  model$independent <- FALSE

  panel_frame <- create_panel_frame(model$formula, data, ctx)
  model$naive_formula <- panel_frame$naive_formula
  model$data <- panel_frame$data
  model$timevar <- ctx_time(ctx)
  model$subset <- subset

  # Cache materialization helpers for fast predict
  mat_cache <- .build_mat_cache(model$formula, data, ctx)
  model$mat_formula <- mat_cache$mat_formula
  model$col_mapping <- mat_cache$col_mapping

  class(model) <- c("linear", class(model))

  model$fit <- function(formula, data, boot, subset, timevar){
    if(!is.null(subset)){
      data <- data[data[[timevar]] >= subset$start & data[[timevar]] <= subset$end]
    }

    if(!is.null(boot)){
      fitted <- bootstraplm(formula, data, type = boot)
    } else{
      fitted <- stats::lm(formula, data)
    }
  }

  model$fitted <- model$fit(model$naive_formula, model$data, model$boot, model$subset, model$timevar)

  model$coefs <- broom::tidy(model$fitted)
  model$gof <- broom::glance(model$fitted)

  model$outcome <- parse_formula(model)$outcome
  model$max_history <- .max_lag_depth(model$formula)
  if (model$max_history < 1L) model$max_history <- 1L

  return(model)
}

#' Get the standard error for prediction
#'
#' @param lmpred
#'
#' @return
#' @keywords internal
#'
#' @examples
get_sepi <- function(lmpred){
  se <- lmpred$se.fit
  scale <- lmpred$residual.scale
  sqrt(se^2 + scale^2)
}

#' Selects a column per row in a matrix
#'
#' @param mat
#' @param column_ids
#'
#' @return
#' @keywords internal
#'
#' @examples
select_col_per_row <- function(mat, column_ids){
  cidx <- cbind(1:nrow(mat), column_ids)
  mat[cidx]
}

#' Get the predictive distribution from a linear model
#'
#' Samples nsamples points from the predictive distribution.
#'
#' @param lmpred
#' @param nsamples
#'
#' @return
#' @keywords internal
#'
#' @examples
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
#' @param model
#' @param data A data.table.
#' @param t Time step to predict.
#' @param ctx A panel_context object.
#' @param what Either "pi" or "expectation".
#' @param ...
#'
#' @return A data.table with key + index + outcome columns.
#' @export
predict.linear <- function(model, data, t, ctx, what = "pi", ...) {
  idx <- ctx_time(ctx)
  grp <- ctx_unit(ctx)
  all_keys <- ctx_keys(ctx)
  max_history <- model$max_history

  # Subset to relevant history window
  data <- data[data[[idx]] <= t & data[[idx]] > (t - max_history - 1)]

  # Materialize formula (use cached helpers to skip update/inject/clean_names)
  data <- materialize_formula(model$formula, data, ctx,
                              .mat_formula = model$mat_formula,
                              .col_mapping = model$col_mapping)

  # Filter for specific time point
  data <- data[data[[idx]] == t]

  # Make predictions
  pred <- predict(model$fitted, newdata = data, se.fit = TRUE)

  # Build result data.table
  result_cols <- c(all_keys, idx, model$outcome)
  result <- data[, ..result_cols]

  if (what == "expectation") {
    data.table::set(result, j = model$outcome, value = pred$fit)
  } else if (what == "pi") {
    data.table::set(result, j = model$outcome, value = getpi(pred))
  } else {
    stop("`what` must be either `pi` or `expectation`")
  }

  return(result)
}
