# Statistical helpers for behavioural tests ----------------------------------
#
# All operate on a list of per-cell draw vectors (`draws`) aligned with a
# `truth` vector. NA-safe: cells with no draws or NA truth are dropped.

# Empirical coverage of the central `level`% predictive interval: the fraction
# of cells whose truth falls inside [q_{a/2}, q_{1-a/2}].
empirical_coverage <- function(draws, truth, level = 90) {
  a <- (1 - level / 100)
  hit <- vapply(seq_along(truth), function(i) {
    d <- draws[[i]]
    d <- d[!is.na(d)]
    if (length(d) == 0L || is.na(truth[i])) return(NA)
    lo <- stats::quantile(d, a / 2, names = FALSE)
    hi <- stats::quantile(d, 1 - a / 2, names = FALSE)
    truth[i] >= lo && truth[i] <= hi
  }, logical(1))
  mean(hit, na.rm = TRUE)
}

# Probability-integral-transform values: the fraction of draws at or below the
# truth, per cell. Should be ~Uniform(0, 1) for a calibrated forecast.
pit_values <- function(draws, truth) {
  out <- vapply(seq_along(truth), function(i) {
    d <- draws[[i]]
    d <- d[!is.na(d)]
    if (length(d) == 0L || is.na(truth[i])) return(NA_real_)
    mean(d <= truth[i])
  }, numeric(1))
  out[!is.na(out)]
}

# Mean predictive-interval width (q_{1-a/2} - q_{a/2}) over cells.
mean_interval_width <- function(draws, level = 90) {
  a <- (1 - level / 100)
  w <- vapply(draws, function(d) {
    d <- d[!is.na(d)]
    if (length(d) == 0L) return(NA_real_)
    unname(diff(stats::quantile(d, c(a / 2, 1 - a / 2))))
  }, numeric(1))
  mean(w, na.rm = TRUE)
}

# CRPS skill score of forecast `a` relative to baseline `b`: 1 - mean(a)/mean(b).
# Positive => `a` is the better (lower-CRPS) forecast.
crps_skill <- function(a, b) {
  1 - mean(a, na.rm = TRUE) / mean(b, na.rm = TRUE)
}

# Run `f(seed)` over `seeds` (an integer count -> 1:n, or an explicit vector)
# and return the results. Scalar results are simplified to a numeric vector.
replicate_seeds <- function(seeds, f) {
  if (length(seeds) == 1L && seeds >= 1) seeds <- seq_len(seeds)
  res <- lapply(seeds, function(s) {
    set.seed(s)
    f(s)
  })
  if (all(vapply(res, function(x) is.numeric(x) && length(x) == 1L, logical(1)))) {
    return(unlist(res))
  }
  res
}
