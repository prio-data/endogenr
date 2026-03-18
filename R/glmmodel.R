#' A bootstrapped GLM
#'
#' Supports residual based bootstrapping and wild bootstrapping on the link scale.
#'
#' @param formula
#' @param data
#' @param family
#' @param type
#'
#' @return
#' @export
#'
#' @examples
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


#' GLM model
#'
#' @param formula
#' @param family A family object (e.g. quasibinomial(), gaussian(), poisson())
#' @param boot
#' @param data
#' @param subset
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
glmmodel <- function(formula = NULL, family = stats::gaussian(), boot = NULL, data = NULL, subset = NULL, ...){
  model <- new_endogenmodel(formula)
  model$boot <- boot
  model$family <- family
  model$fit_args <- rlang::list2(...)
  model$independent <- FALSE

  panel_frame <- create_panel_frame(model$formula, data)
  model$naive_formula <- panel_frame$naive_formula
  model$data <- panel_frame$data
  model$timevar <- tsibble::index_var(data)
  model$subset <- subset
  class(model) <- c("glm_endogenr", class(model))

  model$fit <- function(formula, data, family, boot, subset, timevar){
    if(!is.null(subset)){
      data <- data |> dplyr::filter(!!rlang::sym(timevar) >= subset$start, !!rlang::sym(timevar) <= subset$end)
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

  return(model)
}

#' Get the predictive distribution from a GLM
#'
#' Samples nsamples points from the predictive distribution on the response scale.
#' Sampling is done on the link scale and then transformed via the inverse link function.
#'
#' @param glmpred prediction object from predict.glm with se.fit = TRUE and type = "link"
#' @param family a family object
#' @param df residual degrees of freedom
#' @param nsamples
#'
#' @return
#' @export
#'
#' @examples
getpi_glm <- function(glmpred, family, df, nsamples = 1){
  # Sample on the link scale using t-distribution, then apply inverse link
  eta_samples <- glmpred$fit + outer(glmpred$se.fit, stats::rt(nsamples, df))
  pi <- family$linkinv(eta_samples)
  if(nsamples == 1){
    pi <- as.vector(pi)
  }
  return(pi)
}

#' Predict function for a GLM model
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
predict.glm_endogenr <- function(model, data, t, what = "pi", ...) {
  # Get index and key variables
  idx <- tsibble::index_var(data)
  grp <- tsibble::key_vars(data)

  # Find any numbers in the RHS and use that as the max number of past time periods to include.
  # Default to 1 when the RHS contains no numeric literals.
  nums <- stringr::str_extract_all(model$formula |> as.character(), "[0-9]+")[[3]] |> as.numeric()
  max_history <- if (length(nums) > 0) max(nums) else 1L

  data <- data |>
    dplyr::filter(!!rlang::sym(idx) <= t, !!rlang::sym(idx) > (t-max_history-1))

  # Create panel frame using the tsibble version
  data <- create_panel_frame(model$formula, data)$data

  # Filter for specific time point
  data <- data |>
    dplyr::filter(!!rlang::sym(idx) == t)

  # Predict on link scale with standard errors
  pred <- predict(model$fitted, newdata = data, type = "link", se.fit = TRUE)

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
      dplyr::mutate(!!rlang::sym(model$outcome) := model$family$linkinv(pred$fit))
  } else if (what == "pi") {
    result <- result |>
      dplyr::mutate(
        !!rlang::sym(model$outcome) := getpi_glm(pred, model$family, model$fitted$df.residual)
      )
  } else{
    stop("`what` must be either `pi` or `expectation`")
  }

  return(result)
}
