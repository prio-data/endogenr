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
heterolmmodel <- function(formula = NULL, variance = NULL, data = NULL,
                          ctx = NULL, subset = NULL, ...) {
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
  grp_keys <- ctx_keys(ctx)
  timevar  <- ctx_time(ctx)
  model$timevar <- timevar
  model$subset  <- subset

  # Get outcome name
  outcome_name <- base::all.vars(rlang::f_lhs(formula))

  # Build a combined formula covering all terms from both mean and variance
  # formulas. This is materialised once so ts columns are shared.
  mean_terms <- labels(stats::terms(formula))
  var_terms  <- labels(stats::terms(model$variance_formula))
  all_terms  <- unique(c(mean_terms, var_terms))
  combined_formula <- stats::reformulate(all_terms, response = outcome_name,
                                         env = rlang::f_env(formula))

  # Stage 1: per-unit ts materialisation of the combined formula. `poly(dem,2)`
  # and other design operators are NOT extracted here — they stay in the formula
  # and are evaluated pooled by model.matrix() in Stage 2.
  pm        <- panel_materialize(combined_formula, data,
                                 groupvar = grp_keys, timevar = timevar)
  alias_map <- .pt_make_aliases(pm$map)
  old <- intersect(names(alias_map), names(pm$data))
  if (length(old) > 0L) data.table::setnames(pm$data, old, alias_map[old])

  model$ts_map       <- pm$map
  model$pt_alias_map <- alias_map

  # Rewrite the mean and variance formulas through the combined ts_map, then
  # apply aliases. `.rewrite_from_map()` re-uses the same map built for the
  # combined formula so no ts expression is extracted twice.
  mean_rw <- .rewrite_from_map(formula,                  pm$map)
  var_rw  <- .rewrite_from_map(model$variance_formula,   pm$map)
  mean_fit_formula <- .pt_alias_formula(mean_rw, alias_map)
  var_fit_formula  <- .pt_alias_formula(var_rw,  alias_map)

  # Stage 2: model.matrix expansion for heterolm column-name lookups.
  # heterolm::hetero() does not use R's formula algebra — interaction terms
  # (a:b) and factor dummy codes must be pre-computed as flat columns.
  mat_data <- pm$data
  mean_labels <- labels(stats::terms(mean_fit_formula))
  var_labels  <- labels(stats::terms(var_fit_formula))

  # Clean outcome name (apply alias map for consistency)
  outcome_clean <- if (outcome_name %in% names(alias_map)) {
    alias_map[[outcome_name]]
  } else {
    outcome_name
  }

  mean_exp <- .hetero_expand_terms(mean_labels, mat_data)
  var_exp  <- .hetero_expand_terms(var_labels,  mat_data)

  # Store terms + xlevels for coherent predict-time reconstruction
  model$hetero_mean_terms   <- mean_exp$terms_obj
  model$hetero_mean_xlevels <- mean_exp$xlevels
  model$hetero_var_terms    <- var_exp$terms_obj
  model$hetero_var_xlevels  <- var_exp$xlevels

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

  # Store data (with newly-expanded columns) and the combined formula for
  # predict-time ts rematerialisation
  model$data             <- mat_data
  model$combined_formula <- combined_formula

  # Apply subset if provided
  fit_data <- model$data
  if (!is.null(subset)) {
    fit_data <- fit_data[fit_data[[timevar]] >= subset$start &
                           fit_data[[timevar]] <= subset$end]
  }

  # Fit the model
  fit_args <- c(
    list(formula = model$naive_hetero_formula, data = as.data.frame(fit_data),
         panel.id = NULL),
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
  idx      <- ctx_time(ctx)
  all_keys <- ctx_keys(ctx)

  # Subset to the per-unit history window the mean/variance RHS actually needs.
  data <- .history_subset(data, idx, t, model$required_history)

  # Re-materialise ts columns per (unit, sim) group.
  env <- rlang::f_env(model$combined_formula)
  mat <- .apply_ts_map(model$ts_map, data, all_keys, idx, env = env)
  old <- intersect(names(model$pt_alias_map), names(mat))
  if (length(old) > 0L) data.table::setnames(mat, old, model$pt_alias_map[old])

  # Filter to the prediction time step.
  mat <- mat[mat[[idx]] == t]

  # Re-expand interaction/factor columns using the stored training-time terms
  # objects (coherent basis: same poly scaling, factor contrasts, etc. as fit).
  .hetero_expand_from_terms(model$hetero_mean_terms, model$hetero_mean_xlevels, mat)
  .hetero_expand_from_terms(model$hetero_var_terms,  model$hetero_var_xlevels,  mat)

  # Make predictions using heterolm (returns list with mu and sigma)
  pred <- predict(model$fitted, newdata = as.data.frame(mat), type = "response")

  # Build result data.table with only necessary columns
  result_cols <- c(all_keys, idx, model$outcome)
  result <- mat[, ..result_cols]

  # Update outcome column based on prediction type
  if (what == "expectation") {
    data.table::set(result, j = model$outcome, value = pred$mu)
  } else if (what == "pi") {
    data.table::set(result, j = model$outcome,
                    value = stats::rnorm(nrow(result), pred$mu, pred$sigma))
  } else {
    stop("`what` must be either `pi` or `expectation`")
  }

  return(result)
}
