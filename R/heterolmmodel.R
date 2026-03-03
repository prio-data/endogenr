#' Heteroscedastic linear model
#'
#' Fits a heteroscedastic linear model using [heterolm::hetero()], where both the mean
#' and variance are modeled as functions of covariates. The variance model uses a
#' log-link, so observation-specific prediction intervals are wider where predicted
#' variance is higher.
#'
#' @param formula A two-sided formula for the mean equation (e.g. \code{y ~ lag(x1) + lag(x2)}).
#' @param variance A one-sided formula for the log-variance equation (e.g. \code{~ lag(z1) + lag(z2)}).
#'   Defaults to \code{~ 1} (homoscedastic).
#' @param data A tsibble.
#' @param subset Optional list with \code{start} and \code{end} for subsetting training data.
#' @param ... Additional arguments passed to [heterolm::hetero()].
#'
#' @return An endogenmodel of class \code{heterolm}.
#' @export
heterolmmodel <- function(formula = NULL, variance = NULL, data = NULL, subset = NULL, ...) {
  if (!requireNamespace("heterolm", quietly = TRUE)) {
    stop("Package 'heterolm' is required for heterolm models. Install it from GitHub.")
  }

  model <- new_endogenmodel(formula)
  model$independent <- FALSE
  model$variance_formula <- if (is.null(variance)) ~ 1 else variance
  model$fit_args <- rlang::list2(...)

  # Get tsibble metadata
  grp <- tsibble::key_vars(data)
  idx <- tsibble::index_var(data)
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
  combined_cpf <- create_panel_frame(combined_formula, data)

  # Run create_panel_frame on mean-only formula to identify mean columns
  mean_cpf <- create_panel_frame(formula, data)

  # Determine naive column names for mean and variance
  mean_naive_names <- setdiff(names(mean_cpf$data), c(grp, idx, outcome_name))
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

  # Apply subset if provided
  fit_data <- model$data
  if (!is.null(subset)) {
    fit_data <- fit_data |> dplyr::filter(!!rlang::sym(idx) >= subset$start, !!rlang::sym(idx) <= subset$end)
  }

  # Fit the model
  fit_args <- c(
    list(formula = model$naive_hetero_formula, data = as.data.frame(fit_data), panel.id = NULL),
    model$fit_args
  )
  model$fitted <- do.call(heterolm::hetero, fit_args)

  model$outcome <- outcome_name

  class(model) <- c("heterolm", class(model))
  return(model)
}


#' Predict function for a heteroscedastic linear model
#'
#' Samples from N(mu_i, sigma_i) per observation, where sigma_i comes from the
#' heterolm variance model.
#'
#' @param model A heterolm model object.
#' @param data A tsibble with simulation data.
#' @param t The time point to predict for.
#' @param what Either \code{"pi"} (prediction interval sample) or \code{"expectation"}.
#' @param ... Not used.
#'
#' @return A tsibble with key + index + outcome column.
#' @export
predict.heterolm <- function(model, data, t, what = "pi", ...) {
  # Get index and key variables
  idx <- tsibble::index_var(data)
  grp <- tsibble::key_vars(data)

  # Find max_history from both formulas
  mean_nums <- stringr::str_extract_all(as.character(model$formula)[[3]], "[0-9]+")[[1]] |> as.numeric()
  var_char <- as.character(model$variance_formula)
  var_nums <- stringr::str_extract_all(var_char[length(var_char)], "[0-9]+")[[1]] |> as.numeric()
  all_nums <- c(mean_nums, var_nums)
  max_history <- if (length(all_nums) > 0) max(all_nums) else 1

  data <- data |>
    dplyr::filter(!!rlang::sym(idx) <= t, !!rlang::sym(idx) > (t - max_history - 1))

  # Create panel frame using the combined formula
  data <- create_panel_frame(model$combined_formula, data)$data

  # Filter for specific time point
  data <- data |>
    dplyr::filter(!!rlang::sym(idx) == t)

  # Make predictions using heterolm (returns data.table with mu and sigma)
  pred <- predict(model$fitted, newdata = as.data.frame(data), type = "response")

  # Create result tsibble with only necessary columns
  result <- data |>
    dplyr::select(!!!rlang::syms(c(grp, idx, model$outcome))) |>
    tsibble::as_tsibble(
      key = grp,
      index = idx
    )

  # Update outcome column based on prediction type
  if (what == "expectation") {
    result <- result |>
      dplyr::mutate(!!rlang::sym(model$outcome) := pred$mu)
  } else if (what == "pi") {
    result <- result |>
      dplyr::mutate(
        !!rlang::sym(model$outcome) := stats::rnorm(dplyr::n(), pred$mu, pred$sigma)
      )
  } else {
    stop("`what` must be either `pi` or `expectation`")
  }

  return(result)
}
