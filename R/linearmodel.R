#' A bootstrapped linear model
#'
#' Supports residual based bootstrapping and wild bootstrapping.
#'
#' @param formula
#' @param data
#' @param type
#'
#' @return
#' @export
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



#' Linear model
#'
#' @param formula
#' @param boot
#' @param data
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
linearmodel <- function(formula = NULL, boot = NULL, data = NULL, subset = NULL, ...){
  model <- new_endogenmodel(formula)
  model$boot <- boot
  model$fit_args <- rlang::list2(...)
  model$independent <- FALSE

  panel_frame <- create_panel_frame(model$formula, data)
  model$naive_formula <- panel_frame$naive_formula
  model$data <- panel_frame$data
  model$timevar <- tsibble::index_var(data)
  model$subset <- subset
  class(model) <- c("linear", class(model))

  model$fit <- function(formula, data, boot, subset, timevar){
    if(!is.null(subset)){
      data <- data |> dplyr::filter(!!rlang::sym(timevar) >= subset$start, !!rlang::sym(timevar) <= subset$end)
    }

    if(!is.null(boot)){
      fitted <- bootstraplm(formula, data, type = boot)
    } else(
      fitted = stats::lm(formula, data)
    )
  }

  model$fitted <- model$fit(model$naive_formula, model$data, model$boot, model$subset, model$timevar)

  model$coefs <- broom::tidy(model$fitted)
  model$gof <- broom::glance(model$fitted)

  model$outcome <- parse_formula(model)$outcome

  return(model)
}

#' Get the standard error for prediction
#'
#' @param lmpred
#'
#' @return
#' @export
#'
#' @examples
get_sepi <- function(lmpred){
  #lmpred <- predict(lmfit, se.fit = T)
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
#' @export
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
#' @export
#'
#' @examples
getpi <- function(lmpred, nsamples = 1){
  sepi <- get_sepi(lmpred)
  pi <- lmpred$fit + outer(sepi, stats::rt(nsamples, lmpred$df))
  if(nsamples == 1){
    pi <- as.vector(pi)
  }
  return(pi)
}

#' Predict function for a linear model
#'
#' @param model
#' @param data
#' @param t
#' @param what
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
predict.linear <- function(model, data, t, what = "pi", ...) {
  # Get index and key variables
  idx <- tsibble::index_var(data)
  grp <- tsibble::key_vars(data)

  # Find any numbers in the RHS and use that as the max number of past time periods to include.
  max_history <- stringr::str_extract_all(model$formula |> as.character(), "[0-9]+")[[3]] |> as.numeric() |> max()

  data <- data |>
    dplyr::filter(!!rlang::sym(idx) <= t, !!rlang::sym(idx) > (t-max_history-1))

  # Create panel frame using the tsibble version
  data <- create_panel_frame(model$formula, data)$data

  # Filter for specific time point
  data <- data |>
    dplyr::filter(!!rlang::sym(idx) == t)

  # Make predictions
  if (!is.null(model$boot)) {
    pred <- predict(model$fitted, newdata = data, se.fit = TRUE)
  } else {
    pred <- predict(model$fitted, newdata = data, se.fit = TRUE)
  }

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
      dplyr::mutate(!!rlang::sym(model$outcome) := pred$fit)
  } else if (what == "pi") {
    result <- result |>
      dplyr::mutate(
        !!rlang::sym(model$outcome) := getpi(pred)
      )
  } else{
    stop("`what´ must be either `pi´ or `expectation´")
  }

  return(result)
}
