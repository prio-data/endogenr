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
  # Default to an intercept-only log-variance. Keep this a one-sided formula
  # (not the literal `1`) so terms()/.max_lag_depth()/.edges_from_formula() all
  # accept it; `labels(terms(~1))` is empty, which the naive-name logic below
  # already maps to the intercept "1".
  model$variance_formula <- if (is.null(variance)) ~1 else variance
  model$fit_args <- rlang::list2(...)

  # Get panel metadata from context
  idx <- ctx_time(ctx)
  model$timevar <- idx
  model$subset <- subset

  # Get outcome name
  outcome_name <- base::all.vars(rlang::f_lhs(formula))

  # Build combined formula with all terms from both mean and variance
  mean_terms <- labels(stats::terms(formula))
  var_terms  <- labels(stats::terms(model$variance_formula))
  all_terms  <- unique(c(mean_terms, var_terms))
  combined_formula <- stats::reformulate(all_terms, response = outcome_name)
  model$combined_formula <- combined_formula

  # Materialize combined formula once; reuse cache for predict
  combined_cpf <- create_panel_frame(combined_formula, data, ctx)
  cm <- combined_cpf$col_mapping

  # Derive mean and variance term labels preserving interaction structure.
  # Store the raw cleaned labels (before model.matrix expansion) for the
  # predict path to re-use.
  mean_terms_clean <- .clean_term_labels(formula, cm)
  var_terms_clean  <- .clean_term_labels(model$variance_formula, cm)
  model$hetero_mean_clean_labels <- mean_terms_clean
  model$hetero_var_clean_labels  <- var_terms_clean

  # Clean outcome name (usually unchanged, but apply mapping for consistency)
  outcome_clean <- if (outcome_name %in% names(cm)) cm[[outcome_name]] else outcome_name

  # heterolm::hetero() uses column-name lookups, not R formula algebra, so
  # interactions (g:x) and factor dummies must be pre-computed as flat columns.
  mat_data <- combined_cpf$data
  mean_exp <- .hetero_expand_terms(mean_terms_clean, mat_data)
  var_exp  <- .hetero_expand_terms(var_terms_clean,  mat_data)
  # mat_data was modified in-place; update model$data below

  mean_col_names <- if (length(mean_exp$col_names) == 0L) "1" else mean_exp$col_names
  var_col_names  <- if (length(var_exp$col_names)  == 0L) "1" else var_exp$col_names

  # Build heterolm formula: outcome ~ mean_vars | var_vars
  hetero_formula_str <- paste0(
    outcome_clean, " ~ ",
    paste(mean_col_names, collapse = " + "),
    " | ",
    paste(var_col_names, collapse = " + ")
  )
  model$naive_hetero_formula <- stats::as.formula(hetero_formula_str)

  # Store data (with any newly-expanded columns) and cache
  model$data        <- mat_data
  model$mat_formula <- combined_cpf$mat_formula
  model$col_mapping <- cm

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
  model$required_history <- max(.required_history(model$formula),
                                .required_history(model$variance_formula))

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

  # Subset to the per-unit history window the mean/variance RHS actually needs
  data <- .history_subset(data, idx, t, model$required_history)

  # Materialize combined formula (use cached helpers to skip update/inject/clean_names)
  data <- materialize_formula(model$combined_formula, data, ctx,
                              .mat_formula = model$mat_formula,
                              .col_mapping = model$col_mapping)

  # Filter for specific time point
  data <- data[data[[idx]] == t]

  # Re-expand interaction/factor columns that heterolm::hetero() needs as flat
  # column-name lookups (mirrors what was done at fit time).
  .hetero_expand_terms(model$hetero_mean_clean_labels, data)
  .hetero_expand_terms(model$hetero_var_clean_labels,  data)

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
