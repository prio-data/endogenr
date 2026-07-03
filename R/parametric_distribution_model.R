#' Location-scale Student-t distribution
#'
#' Density, distribution function, quantile function and random generation
#' for the location-scale Student-t distribution with `df` degrees of
#' freedom, location `mu`, and scale `sigma`.
#'
#' These functions are exported so that [fitdistrplus::fitdist()] can resolve
#' them by name (`"t_ls"`) from the search path — the same visibility contract
#' `actuar` uses for its distributions. When endogenr is not attached, the
#' package makes them temporarily visible to `fitdist()` automatically during
#' fitting.
#'
#' @param x,q Numeric vector of quantiles.
#' @param p Numeric vector of probabilities.
#' @param n Number of observations.
#' @param df Degrees of freedom.
#' @param mu Location parameter.
#' @param sigma Scale parameter.
#'
#' @return `dt_ls` gives the density, `pt_ls` the distribution function,
#'   `qt_ls` the quantile function, and `rt_ls` random deviates.
#' @examples
#' x <- rt_ls(500, df = 5, mu = 2, sigma = 1.5)
#' dt_ls(0, df = 5, mu = 2, sigma = 1.5)
#' if (requireNamespace("fitdistrplus", quietly = TRUE)) {
#'   fitdistrplus::fitdist(x, "t_ls",
#'                         start = list(df = 10, mu = mean(x), sigma = sd(x)))
#' }
#' @name t_ls
NULL

#' @rdname t_ls
#' @export
dt_ls <- function(x, df = 1, mu = 0, sigma = 1) {
  1 / sigma * stats::dt((x - mu) / sigma, df)
}

#' @rdname t_ls
#' @export
pt_ls <- function(q, df = 1, mu = 0, sigma = 1) {
  stats::pt((q - mu) / sigma, df)
}

#' @rdname t_ls
#' @export
qt_ls <- function(p, df = 1, mu = 0, sigma = 1) {
  stats::qt(p, df) * sigma + mu
}

#' @rdname t_ls
#' @export
rt_ls <- function(n, df = 1, mu = 0, sigma = 1) {
  stats::rt(n, df) * sigma + mu
}

#' Make package distribution functions visible to fitdistrplus
#'
#' Evaluates `expr` with this package's d/p/q/r functions for `distribution`
#' temporarily visible from the global search path. `fitdistrplus::fitdist()`
#' resolves distribution functions by name from its own namespace chain
#' (`namespace:fitdistrplus -> imports -> base -> globalenv -> search path`),
#' so endogenr-internal functions are invisible to it unless the package is
#' attached. This shim injects only functions that exist in the endogenr
#' namespace and are not already visible, restores any shadowed non-function
#' binding, and always cleans up on exit — safe under `future` workers (each
#' process has its own globalenv, and injections are removed before return).
#'
#' @param distribution Character distribution name (e.g. `"t_ls"`).
#' @param expr Expression to evaluate (typically the `fitdist()` call).
#'
#' @return The value of `expr`.
#' @keywords internal
.with_dist_visible <- function(distribution, expr) {
  fnames <- paste0(c("d", "p", "q", "r"), distribution)
  ns <- parent.env(environment())  # endogenr namespace
  inject <- character(0)
  shadowed <- list()
  for (f in fnames) {
    have_ours <- exists(f, envir = ns, inherits = FALSE)
    visible   <- exists(f, envir = globalenv(), mode = "function")  # globalenv + search path
    if (have_ours && !visible) {
      if (exists(f, envir = globalenv(), inherits = FALSE)) {
        shadowed[[f]] <- get(f, envir = globalenv(), inherits = FALSE)  # non-function binding
      }
      assign(f, get(f, envir = ns, inherits = FALSE), envir = globalenv())
      inject <- c(inject, f)
    }
  }
  on.exit({
    rm(list = inject, envir = globalenv())
    for (f in names(shadowed)) assign(f, shadowed[[f]], envir = globalenv())
  }, add = TRUE)
  force(expr)
}

#' Automatic starting values for distributions fitdist cannot self-start
#'
#' Returns method-of-moments starting values for distributions outside
#' `fitdistrplus:::startargdefault`'s registry, computed on the actual
#' training-window data so they are window-aware. Returns `NULL` when
#' fitdist's own defaults cover the distribution.
#'
#' @param distribution Character distribution name.
#' @param x Numeric vector of training observations.
#'
#' @return A named list of starting values, or `NULL`.
#' @keywords internal
.dist_start <- function(distribution, x) {
  if (distribution == "t_ls") {
    x <- x[is.finite(x)]
    m <- stats::median(x)
    s <- stats::mad(x)
    if (!is.finite(s) || s <= 0) s <- stats::sd(x)
    z <- (x - mean(x)) / stats::sd(x)
    k <- mean(z^4, na.rm = TRUE) - 3            # excess kurtosis
    # moment-matched df: excess kurtosis of t is 6 / (df - 4)
    df <- if (is.finite(k) && k > 0.1) 6 / k + 4 else 30
    return(list(df = min(max(df, 2.5), 100), mu = m, sigma = s))
  }
  NULL
}

#' Sample from a fitted distribution
#'
#' Generates random samples from the distribution fitted by
#' fitdistrplus::fitdist(), optionally at parameter values other than the
#' MLE point estimates (used for parameter-uncertainty draws).
#'
#' @param fitobj A fitdist object.
#' @param n Integer. Number of samples to draw.
#' @param est Named list of parameter values; defaults to the point estimates.
#'
#' @return Numeric vector of random samples.
#' @keywords internal
.sample_from_fitdist <- function(fitobj, n, est = as.list(fitobj$estimate)) {
  dname <- fitobj$distname

  if (dname == "nbinom") {
    prob <- est$size / (est$size + est$mu)
    return(stats::rnbinom(n, size = est$size, prob = prob))
  }

  # Standard distributions (rnorm, rcauchy, rgamma, rpois, ...) and package
  # distributions such as rt_ls, which match.fun resolves inside the
  # endogenr namespace.
  rname <- paste0("r", dname)
  rfun <- tryCatch(match.fun(rname), error = function(e) NULL)

  if (is.null(rfun)) {
    # Try actuar for gumbel etc.
    if (requireNamespace("actuar", quietly = TRUE)) {
      rfun <- tryCatch(utils::getFromNamespace(rname, "actuar"), error = function(e) NULL)
    }
  }

  if (is.null(rfun)) {
    stop(sprintf("Cannot find random generation function '%s' for distribution '%s'", rname, dname))
  }

  do.call(rfun, c(list(n = n), est))
}

#' Fits a parametric distribution using fitdistrplus::fitdist
#'
#' Builds the `fitdist()` argument list from the model spec, derives automatic
#' starting values when none are given (see [.dist_start()]), and evaluates the
#' fit with the package's d/p/q/r functions visible (see [.with_dist_visible()]).
#' `start` may be a named list or a `function(x)` evaluated by fitdist on the
#' training data.
#'
#' @param model An endogenmodel with `$outcome`, `$distribution`, and optional `$fit_args`.
#' @param data A data.frame or data.table containing the outcome column.
#'
#' @return A fitdist object.
#' @keywords internal
fit_parametric_distribution_model <- function(model, data) {
  args <- if (is.null(model$fit_args)) list() else model$fit_args
  args$data  <- data[[model$outcome]]
  args$distr <- model$distribution
  if (is.null(args$start)) {
    args$start <- .dist_start(model$distribution, args$data)
  }

  tryCatch(
    .with_dist_visible(model$distribution, do.call(fitdistrplus::fitdist, args)),
    error = function(e) {
      if (grepl("starting values", conditionMessage(e), fixed = TRUE)) {
        rlang::abort(
          paste0("fitdist() could not derive starting values for distribution '",
                 model$distribution, "'. Pass start = list(...) or ",
                 "start = function(x) list(...) to build_model(); a function is ",
                 "evaluated on the training-window data at fit time."),
          parent = e
        )
      }
      stop(e)
    }
  )
}

#' Parametric distribution model
#'
#' This model is static across time, and therefore independent/exogenous.
#' Uses [fitdistrplus::fitdist()] to fit a distribution to the pooled training
#' data.
#'
#' Use [build_model()] with `type = "parametric_distribution"`, a one-sided
#' formula (e.g. `~gdppc_grwt`), and `distribution = "norm"` (or any
#' distribution name accepted by `fitdist()`).
#'
#' @param spec A `parametric_distribution_spec` object from [build_model()].
#' @param data A data.table or data.frame containing the outcome column.
#' @param ctx A panel_context object.
#' @param ... Additional arguments forwarded to [fitdistrplus::fitdist()].
#'
#' @return An endogenmodel of class `parametric_distribution`.
#' @family simulation
#' @export
#' @exportS3Method
fit_model.parametric_distribution_spec <- function(spec, data = NULL, ctx = NULL, ...) {
  # Extract fitdist-specific args, excluding endogenr-level options
  extra_args <- spec$args[!names(spec$args) %in% c("distribution", "param_uncertainty")]
  do.call(parametric_distribution_model, c(
    list(formula = spec$formula, distribution = spec$args$distribution,
         data = data, ctx = ctx,
         param_uncertainty = isTRUE(spec$args$param_uncertainty)),
    extra_args
  ))
}

#' @keywords internal
parametric_distribution_model <- function(formula = NULL, distribution = NULL, data = NULL,
                                          ctx = NULL, param_uncertainty = FALSE, ...) {
  model <- new_endogenmodel(formula)
  model$distribution <- distribution
  model$fit_args <- rlang::list2(...)
  model$param_uncertainty <- isTRUE(param_uncertainty)

  class(model) <- c("parametric_distribution", class(model))
  model$independent <- TRUE
  model$outcome <- parse_formula(model)$outcome
  model$fitted <- fit_parametric_distribution_model(model, data)

  return(model)
}

#' @rdname draw_predictive
#' @export
draw_predictive.parametric_distribution <- function(model, newdata, n_param = 1L, n_innov = 1L, ...) {
  .check_draw_counts(n_param, n_innov)
  if (n_innov == 0L) {
    stop("n_innov = 0 (conditional mean) is not available for parametric_distribution models",
         call. = FALSE)
  }

  fitobj    <- model$fitted
  n         <- nrow(newdata)
  P         <- max(n_param, 1L)
  K         <- max(n_innov, 1L)
  point_est <- as.list(fitobj$estimate)

  vc <- fitobj$vcov
  use_vcov <- n_param > 0L
  if (use_vcov && (is.null(vc) || !all(is.finite(vc)))) {
    .warn_once(".endogenr_fitdist_no_vcov",
               paste0("fitdist object for distribution '", fitobj$distname,
                      "' carries no usable vcov; parameter draws fall back to ",
                      "the point estimates."))
    use_vcov <- FALSE
  }

  out <- matrix(NA_real_, n, P * K)
  for (j in seq_len(P)) {
    est_j <- point_est
    if (use_vcov) {
      # One MVN(estimate, vcov) parameter draw per column block (asymptotic MLE
      # sampling distribution), shared across all rows of the block.
      est_j <- stats::setNames(as.list(drop(.cf_rmv(1L, fitobj$estimate, vc))),
                               names(fitobj$estimate))
    }
    for (k in seq_len(K)) {
      v <- tryCatch(.sample_from_fitdist(fitobj, n, est = est_j),
                    error = function(e) NULL)
      if (is.null(v) || !all(is.finite(v))) {
        # MVN draws can leave the parameter space (e.g. negative df); fall back
        # to the point estimates for this column.
        if (use_vcov) {
          .warn_once(".endogenr_fitdist_bad_param_draw",
                     paste0("A parameter draw for distribution '", fitobj$distname,
                            "' produced invalid samples; falling back to the ",
                            "point estimates for affected draws."))
        }
        v <- .sample_from_fitdist(fitobj, n, est = point_est)
      }
      out[, (j - 1L) * K + k] <- v
    }
  }
  out
}

#' Predict function for parametric_distribution models
#'
#' Generates random samples from the fitted distribution for all forecast rows
#' via [draw_predictive()]. With `param_uncertainty = TRUE` on the spec, one
#' parameter vector is drawn per simulation from the MLE's asymptotic
#' MVN(estimate, vcov) and shared across that simulation's rows and times;
#' the default conditions on the point estimates.
#'
#' @param model A parametric_distribution endogenmodel.
#' @param data A data.table (the simulation grid).
#' @param ctx A panel_context object.
#' @param test_start Integer. Start of the forecast period.
#' @param horizon Integer. Number of forecast steps.
#' @param inner_sims Integer. Number of inner simulations.
#' @param ... Ignored.
#'
#' @return A data.table with columns: unit, sim, time, outcome.
#' @family simulation
#' @export
predict.parametric_distribution <- function(model, data, ctx, test_start, horizon, inner_sims, ...) {
  unit_col <- ctx_unit(ctx)
  time_col <- ctx_time(ctx)
  all_keys <- ctx_keys(ctx)

  # Filter to forecast period
  pred_data <- .dt_rows(data, data[[time_col]] >= test_start &
                                data[[time_col]] <= (test_start + horizon - 1))

  # Keep only relevant columns
  result_cols <- c(all_keys, time_col, model$outcome)
  result <- pred_data[, ..result_cols]

  if (isTRUE(model$param_uncertainty)) {
    sim_col <- ctx_sim(ctx)
    if (is.null(sim_col)) {
      data.table::set(result, j = model$outcome,
                      value = as.vector(draw_predictive(model, result,
                                                        n_param = 1L, n_innov = 1L)))
    } else {
      result[, (model$outcome) := as.vector(draw_predictive(model, .SD,
                                                            n_param = 1L, n_innov = 1L)),
             by = c(sim_col)]
    }
  } else {
    data.table::set(result, j = model$outcome,
                    value = as.vector(draw_predictive(model, result,
                                                      n_param = 0L, n_innov = 1L)))
  }

  return(result)
}
