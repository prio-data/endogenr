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
