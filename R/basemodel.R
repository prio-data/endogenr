#' A baseline structure for a model
#'
#' @param formula
#'
#' @return
#' @export
#'
#' @examples
new_endogenmodel <- function(formula){
  structure(
    list(
      formula = formula
    ),
    class = "endogenmodel"
  )
}

#' Build a model
#'
#' @param type
#' @param formula
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
build_model <- function(type, formula, ...) {
  dots <- list(...)

  f <- switch(type,
              "deterministic" = purrr::partial(
                deterministicmodel,
                formula = formula
              ),
              "parametric_distribution" = purrr::partial(
                parametric_distribution_model,
                formula = formula,
                distribution = dots$distribution,
                start = dots$start,
                method = ifelse(is.null(dots$method), "mle", dots$method),
                discrete = ifelse(is.null(dots$discrete), FALSE, dots$discrete),
                fix.arg = dots$fix.arg
              ),
              "linear" = purrr::partial(
                linearmodel,
                formula = formula,
                outcome = dots$outcome,
                boot = dots$boot
              ),
              "exogen" = purrr::partial(
                exogenmodel,
                formula = formula
              ),
              "univariate_fable" = purrr::partial(
                univariate_fable_model,
                formula = formula,
                method = dots$method
              ),
              stop("Unknown model type: ", type)
  )

  class(f) <- c(class(f), type)
  return(f)
}
