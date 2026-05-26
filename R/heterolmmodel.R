#' Heteroscedastic linear model
#'
#' Fits a heteroscedastic linear model using [heterolm::hetero()], where both
#' the mean and variance are modelled as functions of covariates. The variance
#' model uses a log-link, so observation-specific prediction intervals are
#' wider where predicted variance is higher.
#'
#' Use [build_model()] with `type = "heterolm"`. Pass the mean equation as the
#' main `formula`, and the (one-sided) log-variance equation through `variance`
#' (defaults to `~ 1`).
#'
#' @param spec A `heterolm_spec` object from [build_model()].
#' @param data A data.table or data.frame.
#' @param ctx A panel_context object.
#' @param subset Optional list with \code{start} and \code{end} for subsetting training data.
#' @param ... Additional arguments forwarded to [heterolm::hetero()].
#'
#' @return An endogenmodel of class \code{heterolm}.
#' @family simulation
#' @export
#' @exportS3Method
fit_model.heterolm_spec <- function(spec, data = NULL, ctx = NULL, subset = NULL, ...) {
  extra_args <- spec$args[!names(spec$args) %in% "variance"]
  do.call(heterolmmodel, c(
    list(formula = spec$formula, variance = spec$args$variance,
         data = data, ctx = ctx, subset = subset),
    extra_args
  ))
}

#' @keywords internal
heterolmmodel <- function(formula = NULL, variance = NULL, data = NULL, ctx = NULL, subset = NULL, ...) {
  if (!requireNamespace("heterolm", quietly = TRUE)) {
    stop("Package 'heterolm' is required for heterolm models. Install it from GitHub.")
  }

  model <- new_endogenmodel(formula)
  model$independent <- FALSE
  model$variance_formula <- if (is.null(variance)) 1 else variance
  model$fit_args <- rlang::list2(...)

  # Get panel metadata from context
  grp <- ctx_unit(ctx)
  idx <- ctx_time(ctx)
  all_keys <- ctx_keys(ctx)
  model$timevar <- idx
  model$subset <- subset

  # Get outcome name
  outcome_name <- base::all.vars(rlang::f_lhs(formula))

  # Build combined formula with all terms from both mean and variance
  mean_terms <- labels(stats::terms(formula))
  var_terms <- labels(stats::terms(model$variance_formula))
  all_terms <- unique(c(mean_terms, var_terms))
  combined_formula <- stats::reformulate(all_terms, response = outcome_name)
  model$combined_formula <- combined_formula

  # Run create_panel_frame on combined formula for the full data
  combined_cpf <- create_panel_frame(combined_formula, data, ctx)

  # Run create_panel_frame on mean-only formula to identify mean columns
  mean_cpf <- create_panel_frame(formula, data, ctx)

  # Determine naive column names for mean and variance
  mean_naive_names <- setdiff(names(mean_cpf$data), c(all_keys, idx, outcome_name))
  var_naive_names <- setdiff(names(combined_cpf$data), names(mean_cpf$data))

  # If variance is intercept-only, use "1"
  if (length(var_naive_names) == 0) {
    var_naive_names <- "1"
  }

  model$naive_mean_names <- mean_naive_names
  model$naive_var_names <- var_naive_names

  # Build heterolm formula: outcome ~ mean_vars | var_vars
  hetero_formula_str <- paste0(
    outcome_name, " ~ ",
    paste(mean_naive_names, collapse = " + "),
    " | ",
    paste(var_naive_names, collapse = " + ")
  )
  model$naive_hetero_formula <- stats::as.formula(hetero_formula_str)

  # Store data
  model$data <- combined_cpf$data

  # Cache materialization helpers for fast predict (using combined formula)
  mat_cache <- .build_mat_cache(combined_formula, data, ctx)
  model$mat_formula <- mat_cache$mat_formula
  model$col_mapping <- mat_cache$col_mapping

  # Apply subset if provided
  fit_data <- model$data
  if (!is.null(subset)) {
    fit_data <- fit_data[fit_data[[idx]] >= subset$start & fit_data[[idx]] <= subset$end]
  }

  # Fit the model
  fit_args <- c(
    list(formula = model$naive_hetero_formula, data = as.data.frame(fit_data), panel.id = NULL),
    model$fit_args
  )
  model$fitted <- do.call(heterolm::hetero, fit_args)

  model$outcome <- outcome_name
  model$max_history <- max(.max_lag_depth(model$formula),
                           .max_lag_depth(model$variance_formula))
  if (model$max_history < 1L) model$max_history <- 1L

  class(model) <- c("heterolm", class(model))
  return(model)
}


#' Predict function for a heteroscedastic linear model
#'
#' Samples from N(mu_i, sigma_i) per observation, where sigma_i comes from the
#' heterolm variance model.
#'
#' @param model A heterolm model object.
#' @param data A data.table.
#' @param t The time point to predict for.
#' @param ctx A panel_context object.
#' @param what Either \code{"pi"} (prediction interval sample) or \code{"expectation"}.
#' @param ... Not used.
#'
#' @return A data.table with key + index + outcome columns.
#' @family simulation
#' @export
predict.heterolm <- function(model, data, t, ctx, what = "pi", ...) {
  idx <- ctx_time(ctx)
  all_keys <- ctx_keys(ctx)
  max_history <- model$max_history

  # Subset to relevant history window
  data <- data[data[[idx]] <= t & data[[idx]] > (t - max_history - 1)]

  # Materialize combined formula (use cached helpers to skip update/inject/clean_names)
  data <- materialize_formula(model$combined_formula, data, ctx,
                              .mat_formula = model$mat_formula,
                              .col_mapping = model$col_mapping)

  # Filter for specific time point
  data <- data[data[[idx]] == t]

  # Make predictions using heterolm (returns list with mu and sigma)
  pred <- predict(model$fitted, newdata = as.data.frame(data), type = "response")

  # Build result data.table with only necessary columns
  result_cols <- c(all_keys, idx, model$outcome)
  result <- data[, ..result_cols]

  # Update outcome column based on prediction type
  if (what == "expectation") {
    data.table::set(result, j = model$outcome, value = pred$mu)
  } else if (what == "pi") {
    data.table::set(result, j = model$outcome, value = stats::rnorm(nrow(result), pred$mu, pred$sigma))
  } else {
    stop("`what` must be either `pi` or `expectation`")
  }

  return(result)
}
