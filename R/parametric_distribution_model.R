#' Creates a distributional distribution object from a fitdistrplus result
#'
#' @param fitobj
#'
#' @return
#' @export
#'
#' @examples
create_distribution_object = function(fitobj) {
  get_nbinom <- function(size, mu){
    prob <- size / (size + mu)
    distributional::dist_negative_binomial(size = size, prob = prob)
  }

  switch(fitobj$distname,
         "norm" = distributional::dist_normal(
           mean = fitobj$estimate["mean"],
           sd = fitobj$estimate["sd"]
         ),
         "cauchy" = distributional::dist_cauchy(
           location = fitobj$estimate["location"],
           scale = fitobj$estimate["scale"]
         ),
         "gumbel" = distributional::dist_gumbel(
           alpha = fitobj$estimate["alpha"],
           scale = fitobj$estimate["scale"]
         ),
         "gamma" = distributional::dist_gamma(
           shape = fitobj$estimate["shape"],
           rate = fitobj$estimate["rate"]
         ),
         "t_ls" = distributional::dist_student_t(
           df = fitobj$estimate["df"],
           mu = fitobj$estimate["mu"],
           sigma = fitobj$estimate["sigma"]
         ),
         "nbinom" = get_nbinom(fit$estimate["size"], fitobj$estimate["mu"]),
         "pois" = distributional::dist_poisson(
           lambda = fitobj$estimate["lambda"]
         ),
         stop("Unsupported distribution")
  )
}

#' Student-t distribution density function
#'
#' @param x
#' @param df
#' @param mu
#' @param sigma
#'
#' @return
#' @export
#'
#' @examples
dt_ls <- function(x, df=1, mu=0, sigma=1){
  1/sigma * stats::dt((x - mu)/sigma, df)
}

#' Student-t distribution distribution function
#'
#' @param x
#' @param df
#' @param mu
#' @param sigma
#'
#' @return
#' @export
#'
#' @examples
pt_ls <- function(q, df=1, mu=0, sigma=1){
  stats::pt((q - mu)/sigma, df)
}

#' Student-t distribution quantile function
#'
#' @param x
#' @param df
#' @param mu
#' @param sigma
#'
#' @return
#' @export
#'
#' @examples
qt_ls <- function(p, df=1, mu=0, sigma=1){
  stats::qt(p, df)*sigma + mu
}

#' Student-t distribution random generation
#'
#' @param x
#' @param df
#' @param mu
#' @param sigma
#'
#' @return
#' @export
#'
#' @examples
rt_ls <- function(n, df=1, mu=0, sigma=1){
  stats::rt(n, df)*sigma + mu
}

#' Fits a parametric distribution using fitdistrplus::fitdist
#'
#' @param model
#' @param data
#'
#' @return
#' @export
#'
#' @examples
fit_parametric_distribution_model <- function(model, data){
  if (is.null(model$fit_args)) {
    fitted <- fitdistrplus::fitdist(data[[model$outcome]], model$distribution)
  } else {
    args <- model$fit_args
    args$data <- data[[model$outcome]]
    args$distr <- model$distribution
    fitted <- do.call(fitdistrplus::fitdist, args)
  }
  create_distribution_object(fitted)
}

#' Parametric distribution model
#'
#' This model is static across time, and therefore independent/exogenous
#'
#' @param formula
#' @param distribution
#' @param data
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
parametric_distribution_model <- function(formula = NULL, distribution = NULL, data = NULL, ...){
  model <- new_endogenmodel(formula)
  model$distribution <- distribution
  model$fit_args <- rlang::list2(...)

  class(model) <- c("parametric_distribution", class(model))
  model$independent <- TRUE
  model$outcome <- parse_formula(model)$outcome
  model$fitted <- fit_parametric_distribution_model(model, data)

  return(model)
}

#' Predict function for parametric_distribution models
#'
#' @param model
#' @param data
#' @param test_start
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
predict.parametric_distribution <- function(model, data, test_start, ...) {
  # Get index and key variables from tsibble
  idx <- tsibble::index_var(data)
  grp <- tsibble::key_vars(data)

  # Create prediction data frame with only necessary columns
  pred_data <- data |>
    dplyr::filter(!!rlang::sym(idx) >= test_start) |>
    dplyr::select(
      !!!rlang::syms(grp),
      !!rlang::sym(idx),
      !!rlang::sym(model$outcome)
    )
  # Generate new outcomes using the fitted distribution
  pred_data <- pred_data |>
    dplyr::mutate(
      !!rlang::sym(model$outcome) := distributional::generate(
        model$fitted,
        dplyr::n()
      )$mean
    ) %>%
    # Ensure tsibble structure
    tsibble::as_tsibble(
      key = grp,
      index = idx
    )

  # Return a fresh copy
  return(pred_data)
}
