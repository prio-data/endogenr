#' Exponential decay since event
#'
#' Computes an exponentially decaying signal from the most recent occurrence
#' of an event. Intended for use inside model formulas in endogenr.
#'
#' At the time of the event, the value is 1. It then decays as
#' \code{exp(-lambda * time_since_event)}. Before any event has occurred,
#' the value equals \code{exp(-lambda * max_years)} (i.e., as if the last
#' event happened \code{max_years} ago). The output is floored at this value
#' so that it never drops below the left-censored default.
#'
#' @param event A numeric vector. An event is considered to occur when
#'   \code{event > 0}.
#' @param lambda Decay rate. Larger values mean faster decay.
#'   Half-life in time units is \code{log(2) / lambda}.
#' @param max_years The assumed time since last event for units with no
#'   observed event. Also used as the floor for the decay.
#'
#' @return A numeric vector the same length as \code{event}.
#' @family formula_helpers
#' @export
#'
#' @examples
#' # In a formula:
#' # build_model("linear",
#' #   formula = y ~ lag(decay_since_event(conflict, 0.3)) + lag(log(gdppc)),
#' #   boot = "resid")
#'
#' x <- c(0, 0, 1, 0, 0, 0, 1, 0, 0)
#' decay_since_event(x, lambda = 0.5)
decay_since_event <- function(event, lambda = 0.5, max_years = 50) {
  floor <- exp(-lambda * max_years)
  n <- length(event)
  out <- rep(floor, n)
  time_since <- NA_real_
  for (i in 1:n) {
    if (!is.na(event[i]) && event[i] > 0) {
      time_since <- 0
    }
    out[i] <- if (!is.na(time_since)) max(exp(-lambda * time_since), floor) else floor
    if (!is.na(time_since)) time_since <- time_since + 1
  }
  out
}

#' Exponential decay time since event
#'
#' The inverse of \code{\link{decay_since_event}}: 0 at the time of the event,
#' rising toward 1 as time passes. Useful for modelling a "peace dividend" or
#' recovery effect.
#'
#' @inheritParams decay_since_event
#'
#' @return A numeric vector the same length as \code{event}.
#' @family formula_helpers
#' @export
#'
#' @examples
#' x <- c(0, 0, 1, 0, 0, 0, 1, 0, 0)
#' time_since_event(x, lambda = 0.5)
time_since_event <- function(event, lambda = 0.5, max_years = 50) {
  ceiling <- 1 - exp(-lambda * max_years)
  n <- length(event)
  out <- rep(ceiling, n)
  time_since <- NA_real_
  for (i in 1:n) {
    if (!is.na(event[i]) && event[i] > 0) {
      time_since <- 0
    }
    out[i] <- if (!is.na(time_since)) min(1 - exp(-lambda * time_since), ceiling) else ceiling
    if (!is.na(time_since)) time_since <- time_since + 1
  }
  out
}

#' Intensity-weighted exponential decay since event
#'
#' Like \code{\link{decay_since_event}}, but the peak value equals the
#' event intensity (e.g., fatalities) rather than 1. Multiple events
#' accumulate: each past event contributes its own decaying intensity.
#'
#' @param event A numeric vector. An event occurs when \code{event > 0}.
#' @param intensity A numeric vector of event intensities (e.g., fatalities).
#'   Only used when \code{event > 0}.
#' @param lambda Decay rate.
#' @param max_years Used to set the floor per-event contribution.
#'
#' @return A numeric vector the same length as \code{event}.
#' @family formula_helpers
#' @export
#'
#' @examples
#' evt <- c(0, 1, 0, 0, 1, 0)
#' fat <- c(0, 100, 0, 0, 500, 0)
#' intensity_decay(evt, fat, lambda = 0.3)
intensity_decay <- function(event, intensity, lambda = 0.5, max_years = 50) {
  n <- length(event)
  out <- rep(0, n)
  for (i in 1:n) {
    total <- 0
    for (j in 1:i) {
      if (!is.na(event[j]) && event[j] > 0) {
        dt <- i - j
        total <- total + intensity[j] * max(exp(-lambda * dt), exp(-lambda * max_years))
      }
    }
    out[i] <- total
  }
  out
}

#' Trailing (expanding or windowed) mean within a time-ordered group
#'
#' Computes, position-by-position over an already time-ordered within-group
#' vector `x`, the trailing mean of the last `window` observations.
#' Intended for use inside model formulas in endogenr, particularly for
#' Mundlak-style between-effects that must update dynamically during simulation.
#'
#' At position \eqn{i}, the value is the mean of
#' \eqn{x[\max(1, i - \text{window} + 1) : i]}, NA-omitted.
#' Returns `NA_real_` until at least `min_obs` non-NA values are available.
#'
#' @param x A numeric vector, time-ordered within one panel unit.
#' @param window `Inf` (expanding/cumulative mean, the default) or a positive
#'   integer specifying the trailing window width.
#' @param min_obs Positive integer. Minimum number of non-NA observations
#'   required before a non-NA value is returned. Default `1L`.
#'
#' @return A numeric vector the same length as `x`.
#' @family formula_helpers
#' @export
#'
#' @examples
#' hist_mean(c(1, 2, 3, 4))                   # c(1, 1.5, 2, 2.5)
#' hist_mean(c(1, 2, 3, 4), window = 2)       # c(1, 1.5, 2.5, 3.5)
#' hist_mean(c(1, NA, 3))                     # c(1, 1, 2)
hist_mean <- function(x, window = Inf, min_obs = 1L) {
  # Validate arguments
  if (!(is.infinite(window) && window > 0) &&
      !(is.numeric(window) && length(window) == 1L && is.finite(window) &&
        window == as.integer(window) && window >= 1L)) {
    stop("`window` must be Inf or a positive integer.", call. = FALSE)
  }
  if (!is.numeric(min_obs) || length(min_obs) != 1L ||
      !is.finite(min_obs) || min_obs != as.integer(min_obs) || min_obs < 1L) {
    stop("`min_obs` must be a positive integer.", call. = FALSE)
  }
  min_obs <- as.integer(min_obs)
  n <- length(x)
  if (n == 0L) return(numeric(0L))

  if (is.infinite(window)) {
    # Expanding mean: fast cumulative path
    cs  <- cumsum(ifelse(is.na(x), 0, x))
    cn  <- cumsum(!is.na(x))
    out <- cs / cn
    out[cn < min_obs] <- NA_real_
    return(out)
  }

  # Finite window: trailing loop (series are short per unit)
  window <- as.integer(window)
  out <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    lo   <- max(1L, i - window + 1L)
    vals <- x[lo:i]
    vals <- vals[!is.na(vals)]
    if (length(vals) >= min_obs) {
      out[i] <- mean(vals)
    }
  }
  out
}

#' Build Mundlak between-effect deterministic specs
#'
#' Constructs a list of `deterministic` model specs that update the
#' Mundlak-style within-group running means (`m_*` columns) each simulation
#' step via [hist_mean()]. Append the result with `c()` into a model system.
#'
#' Each source column and each `m_*` output column must be present in the
#' input data before [setup_system()]. Purely-derived output columns may be
#' NA-filled; source columns must be model-produced or panel keys
#' (checked by [validate_system_closure()]).
#'
#' @param vars Either a character vector of source column names (output names
#'   are `paste0(prefix, vars)`), or a **named** character vector where names
#'   are output column names and values are source expression strings (e.g.
#'   `c(m_grwt_l1 = "grwt_l1", m_c25 = "c25")`).
#' @param window Scalar (`Inf` or a positive integer applied to all variables)
#'   or a named numeric vector keyed by output column name (per-variable
#'   windows). `Inf` yields an expanding (cumulative) mean.
#' @param prefix Character prefix prepended to source names when `vars` is
#'   unnamed. Default `"m_"`.
#' @param min_obs Positive integer passed to [hist_mean()]. Default `1L`.
#' @param bounds Optional two-element numeric vector `c(lower, upper)` applied
#'   to all output columns, or `NULL` (no clamping).
#'
#' @return A list of `deterministic` model specs ready to `c()` into a system.
#' @family build
#' @export
#'
#' @examples
#' \dontrun{
#' sys <- setup_system(
#'   c(
#'     exogen(~x),
#'     mundlak_means("x"),
#'     build_model("linear", formula = y ~ lag(y) + m_x)
#'   ),
#'   data = dt, train_start = 1, test_start = 25, horizon = 5,
#'   groupvar = "unit", timevar = "time"
#' )
#' }
mundlak_means <- function(vars, window = Inf, prefix = "m_", min_obs = 1L,
                           bounds = NULL) {
  # Normalise vars to a named character vector (output -> source expression)
  if (is.null(names(vars))) {
    if (!is.character(vars) || length(vars) == 0L) {
      stop("`vars` must be a non-empty character vector.", call. = FALSE)
    }
    src_exprs <- vars
    out_names <- paste0(prefix, vars)
    vars <- stats::setNames(src_exprs, out_names)
  } else {
    if (length(vars) == 0L) stop("`vars` must be non-empty.", call. = FALSE)
    if (!is.character(vars)) stop("`vars` must be a character vector.", call. = FALSE)
    out_names <- names(vars)
  }

  # Validate output uniqueness
  if (anyDuplicated(out_names)) {
    stop("output names in `vars` must be unique.", call. = FALSE)
  }

  # Normalise window to a named vector keyed by output name
  if (length(window) == 1L) {
    window <- stats::setNames(rep(window, length(out_names)), out_names)
  } else {
    if (!all(out_names %in% names(window))) {
      stop("`window` vector must be named by every output column name.", call. = FALSE)
    }
    window <- window[out_names]
  }

  # Validate each window entry
  for (nm in out_names) {
    w <- window[[nm]]
    ok <- (is.infinite(w) && w > 0) ||
          (is.numeric(w) && length(w) == 1L && is.finite(w) &&
           w == as.integer(w) && w >= 1L)
    if (!ok) stop("`window` entries must be Inf or a positive integer.", call. = FALSE)
  }

  # Build one deterministic spec per output column
  specs <- vector("list", length(out_names))
  for (i in seq_along(out_names)) {
    out  <- out_names[i]
    src  <- vars[[i]]
    w    <- window[[out]]
    w_str <- if (is.infinite(w)) "Inf" else format(w)
    rhs  <- sprintf("I(hist_mean(%s, window = %s, min_obs = %dL))", src, w_str, as.integer(min_obs))
    f    <- stats::reformulate(rhs, response = out)
    environment(f) <- parent.frame()
    specs[[i]] <- build_model("deterministic", formula = f, bounds = bounds)
  }
  specs
}
