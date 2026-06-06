# Coefficient-trajectory forecasting -------------------------------------------
#
# An in-package, dependency-light analogue of time-varying parameter (TVP) /
# state-space regression. We estimate a deterministic coefficient *path* from
# expanding (or rolling) training windows, then project each coefficient forward
# as a Gaussian random walk (optionally with drift), starting from the forecast
# origin's estimation uncertainty (the fit's vcov) and fanning out by the
# empirical evolution covariance of the path. The result is "the distribution of
# the regression parameters that changes over time" that the dynamic simulator's
# random-window bootstrap only approximates pointwise.
#
# This is deliberately the empirical version. The principled estimators are
# state-space TVP models: dlm / KFAS (Kalman, no MCMC) and shrinkTVP / walker
# (Bayesian, shrinkage toward a static coefficient); see forecast_coefficients()
# docs.


# Internal: draw n rows from MVN(mu, Sigma) (or a heavier-tailed multivariate-t
# with `df` degrees of freedom). Robust to a non-PSD Sigma via an eigen clip, so
# pairwise/empirical covariances are safe to pass.
.cf_rmv <- function(n, mu, Sigma, df = NULL) {
  p <- length(mu)
  if (p == 0L) return(matrix(numeric(0), nrow = n, ncol = 0L))
  e <- eigen(Sigma, symmetric = TRUE)
  vals <- pmax(e$values, 0)
  # symmetric square root: V diag(sqrt(vals)) V'
  root <- e$vectors %*% (t(e$vectors) * sqrt(vals))
  z <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
  out <- z %*% root
  if (!is.null(df) && is.finite(df) && df > 0) {
    out <- out * sqrt(df / stats::rchisq(n, df))   # row i scaled by scale[i]
  }
  sweep(out, 2L, mu, "+")
}


# Internal: the grid of window-end anchors tau, always including the forecast
# origin (test_start - 1) and never crossing it (no leakage).
.cf_anchor_grid <- function(train_start, test_start, window, min_train, width, step) {
  origin <- test_start - 1L
  first  <- if (window == "rolling") train_start + width - 1L
            else train_start + min_train - 1L
  if (first > origin) {
    stop("No valid training windows: the first anchor (", first, ") is past the ",
         "forecast origin (", origin, "). Lower `min_train`/`width` or extend ",
         "the training range.", call. = FALSE)
  }
  taus <- seq.int(first, origin, by = as.integer(step))
  if (taus[length(taus)] != origin) taus <- c(taus, origin)
  taus
}


# Internal: estimate the coefficient path for ONE coefficient-bearing spec.
# Bootstrap is forced off so each window yields the OLS/MLE point estimate (a
# clean series), and the origin fit is kept for its coefficients and vcov.
# Returns NULL for specs that carry no tidy coefficients (e.g. heterolm).
.cf_path_one <- function(spec, object, taus, window, width) {
  spec0 <- spec
  spec0$args$boot <- NULL                 # point path, not a bootstrap realisation

  train_start <- object$train_start
  rows <- vector("list", length(taus))
  origin_fit <- NULL

  for (i in seq_along(taus)) {
    tau <- taus[i]
    sub <- if (window == "rolling") list(start = tau - width + 1L, end = tau)
           else list(start = train_start, end = tau)
    m <- .fit_spec(spec0, object, subset = sub)
    if (is.null(m$coefs)) return(NULL)    # not a coefficient-bearing model
    ct <- data.table::as.data.table(m$coefs)
    rows[[i]] <- data.table::data.table(
      outcome   = m$outcome,
      term      = as.character(ct$term),
      .tau      = tau,
      estimate  = as.numeric(ct$estimate),
      std.error = as.numeric(ct$std.error)
    )
    if (i == length(taus)) origin_fit <- m$fitted
  }

  list(path = data.table::rbindlist(rows), origin_fit = origin_fit)
}


# Internal: project ONE outcome's coefficient path forward and summarise.
# `path` is the tidy path for a single outcome; `origin_fit` its last-window fit.
.cf_forecast_one <- function(outcome, path, origin_fit, test_start, horizon,
                             method, nsim, student_t, step) {
  cterms <- names(stats::coef(origin_fit))         # canonical term set (origin)
  cterms <- cterms[!is.na(stats::coef(origin_fit))]
  p <- length(cterms)
  b_m <- stats::coef(origin_fit)[cterms]
  Sigma0 <- as.matrix(stats::vcov(origin_fit))[cterms, cterms, drop = FALSE]

  # Path matrix B: anchors (rows) x canonical terms (cols).
  taus <- sort(unique(path$.tau))
  m <- length(taus)
  if (m < 3L) {
    stop("Need at least 3 window anchors to estimate a coefficient trajectory ",
         "for outcome '", outcome, "' (got ", m, "). Lower `min_train`/`width`, ",
         "reduce `step`, or extend the training range.", call. = FALSE)
  }
  B <- matrix(NA_real_, nrow = m, ncol = p, dimnames = list(NULL, cterms))
  tau_idx <- match(path$.tau, taus)
  col_idx <- match(path$term, cterms)
  keep <- !is.na(col_idx)
  B[cbind(tau_idx[keep], col_idx[keep])] <- path$estimate[keep]

  # Per-year evolution covariance from the path's first differences. Anchors are
  # `step` years apart, so a one-step diff spans `step` years: divide by `step`.
  dB <- diff(B)
  sigvar <- apply(dB, 2L, function(x) {
    x <- x[!is.na(x)]; if (length(x) >= 2L) stats::var(x) else 0
  })
  Sigma_w <- stats::cov(dB, use = "pairwise.complete.obs")
  Sigma_w[is.na(Sigma_w)] <- 0
  diag(Sigma_w) <- sigvar
  Sigma_w <- Sigma_w / step

  # Drift: 0 for a pure random walk; the per-year OLS slope of each term's path
  # for "drift".
  drift <- rep(0, p)
  names(drift) <- cterms
  if (method == "drift") {
    for (j in seq_len(p)) {
      y <- B[, j]; ok <- !is.na(y)
      if (sum(ok) >= 2L) drift[j] <- unname(stats::coef(stats::lm(y[ok] ~ taus[ok]))[2L])
    }
    drift[is.na(drift)] <- 0
  }

  df <- if (isTRUE(student_t)) origin_fit$df.residual else NULL

  # Coherent random-walk trajectories: origin draw carries estimation
  # uncertainty (Sigma0); each yearly step adds an MVN(drift, Sigma_w) innovation
  # so the marginal at step h is N(b_m + h*drift, Sigma0 + h*Sigma_w).
  cur <- .cf_rmv(nsim, b_m, Sigma0, df)
  out_rows <- vector("list", horizon)
  for (k in seq_len(horizon)) {
    cur <- cur + .cf_rmv(nsim, drift, Sigma_w, df)
    t_k <- test_start + k - 1L
    qs  <- apply(cur, 2L, stats::quantile, probs = c(0.05, 0.5, 0.95), names = FALSE)
    out_rows[[k]] <- data.table::data.table(
      outcome = outcome, term = cterms, .time = t_k, .type = "forecast", .h = k,
      .mean = colMeans(cur), .sd = apply(cur, 2L, stats::sd),
      .q05 = qs[1L, ], .q50 = qs[2L, ], .q95 = qs[3L, ],
      .draws = lapply(seq_len(p), function(j) cur[, j])
    )
  }
  fc <- data.table::rbindlist(out_rows)

  # Observed path rows (Normal display band from the per-window std.error).
  obs <- data.table::copy(path[path$term %in% cterms])
  z <- stats::qnorm(0.95)
  obs <- obs[, .(outcome = outcome, term = term, .time = .tau, .type = "observed",
                 .h = NA_integer_, .mean = estimate, .sd = std.error,
                 .q05 = estimate - z * std.error, .q50 = estimate,
                 .q95 = estimate + z * std.error,
                 .draws = vector("list", .N))]
  list(obs = obs, fc = fc)
}


#' Forecast where regression coefficients are heading (time-varying parameters)
#'
#' Estimates a deterministic coefficient *path* for each `linear`/`glm` spec in a
#' fitted (or set-up) system by refitting on a grid of expanding or rolling
#' training windows, then projects every coefficient forward over the forecast
#' horizon as a Gaussian random walk (optionally with drift). The forward
#' distribution starts from the forecast origin's estimation uncertainty (the
#' fit's variance-covariance matrix) and fans out by the empirical evolution
#' covariance of the path, preserving cross-term correlations.
#'
#' @details
#' This is an in-package, dependency-light analogue of **time-varying parameter
#' (TVP) / state-space (dynamic linear) regression**: `y_t = X_t b_t + e_t` with
#' `b_t = b_{t-1} + w_t`. The dynamic simulator's random-window bootstrap
#' (`min_window` in [setup_system()]) approximates the same object pointwise but
#' indexes it by a random window (start *and* end), confounding recency with
#' sample size. [forecast_coefficients()] instead uses a deterministic grid
#' anchored at the window *end*, so the coefficient series is a clean function of
#' time that can be projected forward.
#'
#' The forecast is deliberately the empirical/Gaussian version. The principled
#' estimators are state-space TVP models: `dlm` / `KFAS` (Kalman filter, no
#' MCMC), `shrinkTVP` / `walker` (Bayesian, with shrinkage toward a static
#' coefficient that guards against TVP overfitting), and `bvarsv` (TVP-VAR,
#' Primiceri 2005). Reach for those when you need filtered/smoothed states,
#' shrinkage, or formal posterior inference.
#'
#' Bootstrap is forced off for the path fits (so each window yields the OLS/MLE
#' point estimate, not a single bootstrap realisation); the spec's own `boot`
#' setting is irrelevant here. The coefficient path is **pooled across units**,
#' matching the pooled `linear`/`glm` model. Only coefficient-bearing specs
#' (`linear`/`glm`) are forecast; others (including `heterolm`, which carries no
#' tidy coefficients) are skipped, mirroring [get_coefficients()].
#'
#' @param object An `endogenr_fitted_system` from [fit_system()] or an
#'   `endogenr_system_setup` from [setup_system()]. Either supplies the specs,
#'   training data, panel context, and `train_start`/`test_start`/`horizon`.
#' @param horizon Integer or `NULL`. Number of steps ahead to forecast the
#'   coefficients (`test_start ... test_start + horizon - 1`). Defaults to the
#'   system's `horizon`.
#' @param window One of `"expanding"` (default; window `[train_start, tau]`) or
#'   `"rolling"` (window `[tau - width + 1, tau]`).
#' @param width Integer or `NULL`. Rolling-window length (ignored for
#'   `"expanding"`). Defaults to `min_train`.
#' @param min_train Integer or `NULL`. Minimum training length: the first
#'   expanding window (and the default rolling `width`). Defaults to the system's
#'   `min_window`, or `10` when that is unset.
#' @param step Integer. Spacing (in time units) between window-end anchors.
#'   Default `1`.
#' @param method One of `"rw"` (default; random walk, no trend, fanning
#'   uncertainty) or `"drift"` (random walk plus the per-year OLS slope of the
#'   path, so the central path continues the recent trend).
#' @param nsim Integer. Number of coefficient-trajectory draws used to summarise
#'   the forward distribution. Default `1000`.
#' @param student_t Logical. When `TRUE`, innovations are drawn from a
#'   heavier-tailed multivariate-t using the origin fit's residual degrees of
#'   freedom instead of a Gaussian. Default `FALSE`.
#' @param draws Logical. When `TRUE`, forecast rows carry a `.draws` list-column
#'   of the per-time coefficient draws (length `nsim`). Default `FALSE` drops it
#'   to keep the table lean.
#'
#' @return A `data.table`, one row per `(outcome, term, .time)`, with columns:
#'   `outcome`, `term`, `.time`, `.type` (`"observed"` for the fitted path or
#'   `"forecast"` for the projection), `.h` (steps ahead; `NA` for observed),
#'   `.mean`, `.sd`, `.q05`, `.q50`, `.q95`, and (when `draws = TRUE`) `.draws`.
#'   The observed band uses a Normal approximation from each window's standard
#'   error; the forecast band is the empirical quantiles of the trajectory draws.
#'   Carries `panel_unit`/`panel_time` attributes and the chosen `window`/`method`.
#' @seealso [plot_coefficient_forecast()], [get_coefficients()], [setup_system()],
#'   [fit_system()]
#' @family postprocess
#' @export
forecast_coefficients <- function(object, horizon = NULL,
                                  window = c("expanding", "rolling"),
                                  width = NULL, min_train = NULL, step = 1L,
                                  method = c("rw", "drift"), nsim = 1000L,
                                  student_t = FALSE, draws = FALSE) {
  if (!inherits(object, "endogenr_system_setup")) {
    stop("`object` must be an endogenr_system_setup (setup_system()) or an ",
         "endogenr_fitted_system (fit_system()).", call. = FALSE)
  }
  window <- match.arg(window)
  method <- match.arg(method)
  step   <- as.integer(step)
  nsim   <- as.integer(nsim)
  if (is.na(step) || step < 1L) stop("`step` must be a positive integer.", call. = FALSE)
  if (is.na(nsim) || nsim < 2L) stop("`nsim` must be an integer >= 2.", call. = FALSE)

  if (is.null(horizon)) horizon <- object$horizon
  horizon <- as.integer(horizon)
  if (is.na(horizon) || horizon < 1L) stop("`horizon` must be a positive integer.", call. = FALSE)

  if (is.null(min_train)) min_train <- if (!is.null(object$min_window)) object$min_window else 10L
  min_train <- as.integer(min_train)
  if (is.null(width)) width <- min_train
  width <- as.integer(width)

  taus <- .cf_anchor_grid(object$train_start, object$test_start, window,
                          min_train, width, step)

  specs <- object$specs
  candidate <- which(vapply(specs, function(s) s$type %in% c("linear", "glm"),
                            logical(1)))
  if (length(candidate) == 0L) {
    stop("No `linear`/`glm` specs to forecast: only those carry coefficients.",
         call. = FALSE)
  }

  parts <- list()
  for (k in candidate) {
    res <- .cf_path_one(specs[[k]], object, taus, window, width)
    if (is.null(res)) next                       # no tidy coefficients
    outcome <- res$path$outcome[1L]
    one <- .cf_forecast_one(outcome, res$path, res$origin_fit, object$test_start,
                            horizon, method, nsim, student_t, step)
    parts[[length(parts) + 1L]] <- one$obs
    parts[[length(parts) + 1L]] <- one$fc
  }
  if (length(parts) == 0L) {
    stop("No coefficient-bearing models produced a forecast.", call. = FALSE)
  }

  out <- data.table::rbindlist(parts, use.names = TRUE)
  if (!isTRUE(draws)) out[, .draws := NULL]
  data.table::setorderv(out, c("outcome", "term", ".time"))

  data.table::setattr(out, "panel_unit", ctx_unit(object$fit_ctx))
  data.table::setattr(out, "panel_time", ctx_time(object$fit_ctx))
  data.table::setattr(out, "cf_window", window)
  data.table::setattr(out, "cf_method", method)
  data.table::setattr(out, "cf_test_start", object$test_start)
  class(out) <- c("coef_forecast", class(out))
  out[]
}
