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

  # Finite window: right-aligned trailing mean via windowed cumulative sums.
  window <- as.integer(window)
  cs  <- cumsum(ifelse(is.na(x), 0, x))
  cn  <- cumsum(!is.na(x))
  lo  <- pmax(0L, seq_len(n) - window)        # #elements strictly before each window
  num <- cs - c(0, cs)[lo + 1L]
  den <- cn - c(0, cn)[lo + 1L]
  out <- num / den
  out[den < min_obs] <- NA_real_              # also converts empty-window 0/0 NaN -> NA
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

#' Build time-dummy unit means for two-way Mundlak correction on unbalanced panels
#'
#' Constructs the \eqn{\bar{f}_{r,i}} columns — the fraction of unit
#' \eqn{i}'s complete-case training periods that fall in time period \eqn{r}
#' — which, together with \code{factor(timevar)} and the covariate unit means
#' (\code{m_*}), restore TWFE-equivalent within-slopes on an unbalanced panel
#' (Wooldridge 2025, §10.2).
#'
#' The columns are data-derived constants: fixed over the training window and
#' \strong{set to zero in all forecast rows} so that the \code{exogen} carrier
#' introduces no structural forecast shift. Covariate means (\code{m_*}) remain
#' untouched.
#'
#' @section Exactness caveat:
#' The \code{fbar_*} correction yields coefficients numerically identical to
#' TWFE only when the covariate unit means (\code{m_*}) in \code{formula} are
#' also \strong{constant complete-case means} over the same estimation sample —
#' not the rolling \code{hist_mean()} means produced by \code{\link{mundlak_means}}.
#' Rolling means yield approximate equivalence. Compute constant means before
#' calling this function and pass the resulting data.
#'
#' @section Forecast behaviour:
#' Because \code{fbar_*} columns are ordinary data columns carried forward by
#' the returned \code{exogen} spec, the zero-in-forecast default is set directly
#' in the augmented \code{data}. A caller who wants the training-window constant
#' propagated instead simply overwrites the forecast rows
#' (\code{data[time >= test_start, (tw$terms) := <per-unit const>]}) before
#' passing \code{tw$data} to \code{\link{setup_system}}; no engine change is
#' needed.
#'
#' @param data A data.table or data.frame containing all units (training and
#'   forecast rows). Forecast rows (where \code{timevar >= test_start}) must
#'   exist for every simulated unit; the helper sets their \code{fbar_*}
#'   columns to 0.
#' @param formula The \strong{base} model formula \emph{without} any
#'   \code{fbar_*} terms. It defines the estimation sample via the same
#'   \code{na.omit} logic used by \code{\link[=linearmodel]{linearmodel}}.
#'   Any \code{m_*} columns it references must already be present in
#'   \code{data}.
#' @param groupvar Character. Name of the unit (group) identifier column.
#' @param timevar Character. Name of the time identifier column.
#' @param test_start Scalar numeric. First time period of the forecast window;
#'   rows with \code{timevar >= test_start} receive \code{fbar_* = 0}.
#' @param prefix Character prefix for the generated column names.
#'   Default \code{"fbar_"}.
#' @param drop_ref Logical. If \code{TRUE} (default), the column corresponding
#'   to the earliest training year is dropped to remove the sum-to-1
#'   collinearity with the intercept, matching \code{factor()}'s default
#'   reference-level encoding. Set to \code{FALSE} only when an explicit
#'   reference is managed by the caller.
#'
#' @return A named list with four elements:
#' \describe{
#'   \item{\code{data}}{The input \code{data} augmented with the \code{fbar_*}
#'     columns. Training rows carry the per-unit fraction; forecast rows carry
#'     0.}
#'   \item{\code{terms}}{Character vector of retained \code{fbar_*} column
#'     names (after the reference-level drop).}
#'   \item{\code{spec}}{An \code{exogen} model spec (from
#'     \code{\link{build_model}}) that carries the \code{fbar_*} columns into
#'     the simulation grid.}
#'   \item{\code{formula}}{The base \code{formula} extended with all
#'     \code{terms} as additive predictors.}
#' }
#'
#' @family formula_helpers
#' @export
#'
#' @examples
#' \dontrun{
#' library(data.table)
#' # Build data with constant complete-case unit means of x
#' d <- sim_panel_ragged(units = 8, n_time = 30, seed = 7,
#'                       enter = c("5" = 8), exit = c("6" = 22))
#' d[, m_x := mean(x), by = unit]   # constant cc mean (no NAs in DGP)
#'
#' core <- y ~ x + m_x + factor(time)
#' tw <- mundlak_time_means(d, core, groupvar = "unit", timevar = "time",
#'                          test_start = 31L)
#'
#' # Wire into a simulation system:
#' system <- c(
#'   list(build_model("linear", formula = tw$formula, boot = "resid")),
#'   mundlak_means("x"),           # m_x running mean for simulation dynamics
#'   list(tw$spec),                # fbar_* carrier (0 in forecast rows)
#'   list(build_model("exogen", formula = ~x))
#' )
#' sys <- setup_system(system, tw$data,
#'                     train_start = 1L, test_start = 31L, horizon = 5L,
#'                     groupvar = "unit", timevar = "time")
#' }
mundlak_time_means <- function(data, formula, groupvar, timevar, test_start,
                                prefix = "fbar_", drop_ref = TRUE) {
  # --- Input guards ---
  if (!is.character(groupvar) || length(groupvar) != 1L)
    stop("`groupvar` must be a length-1 character string.", call. = FALSE)
  if (!is.character(timevar) || length(timevar) != 1L)
    stop("`timevar` must be a length-1 character string.", call. = FALSE)
  if (!is.numeric(test_start) || length(test_start) != 1L || !is.finite(test_start))
    stop("`test_start` must be a finite scalar numeric.", call. = FALSE)
  if (!is.character(prefix) || length(prefix) != 1L)
    stop("`prefix` must be a length-1 character string.", call. = FALSE)
  if (!is.logical(drop_ref) || length(drop_ref) != 1L)
    stop("`drop_ref` must be a length-1 logical.", call. = FALSE)

  # Warn if formula does not contain factor(<timevar>); fbar_* is only
  # TWFE-equivalent alongside the time dummies.
  term_labels <- attr(stats::terms(formula), "term.labels")
  has_time_fe <- any(vapply(term_labels, function(lbl) {
    expr <- tryCatch(str2lang(lbl), error = function(e) NULL)
    if (is.null(expr) || !is.call(expr)) return(FALSE)
    fn <- as.character(expr[[1L]])
    if (!fn %in% c("factor", "stats::factor")) return(FALSE)
    if (length(expr) < 2L) return(FALSE)
    arg1 <- expr[[2L]]
    is.symbol(arg1) && as.character(arg1) == timevar
  }, logical(1L)))
  if (!has_time_fe) {
    warning(
      "`formula` does not contain `factor(", timevar, ")`. ",
      "The `fbar_*` correction achieves TWFE equivalence only when time ",
      "dummies are also included.", call. = FALSE
    )
  }

  # --- 1. Copy to data.table; do not mutate caller ---
  data <- data.table::copy(data.table::as.data.table(data))
  data.table::setkeyv(data, c(groupvar, timevar))

  # --- 2. Training subset ---
  train <- data[data[[timevar]] < test_start]

  # --- 3. Reproduce the linear model's estimation sample exactly.
  #        panel_materialize() is the identical call linearmodel() makes at
  #        R/linearmodel.R:97-98; na.omit over fit_vars matches .lm_stage2_fit.
  pm       <- panel_materialize(formula, train,
                                groupvar = groupvar, timevar = timevar)
  fit_vars <- intersect(all.vars(pm$formula), names(pm$data))
  mask     <- stats::complete.cases(
    as.data.frame(pm$data)[, fit_vars, drop = FALSE]
  )

  # --- 4. Years with at least one complete case ---
  yrs <- sort(unique(train[[timevar]][mask]))
  if (length(yrs) == 0L)
    stop("No complete cases found in the training window.", call. = FALSE)

  # --- 5. Per-unit fractions over complete cases.
  #        Count (unit, year) cells then divide by each unit's T_i.
  cc_train  <- train[mask]
  freq_long <- cc_train[, .(count = .N), by = c(groupvar, timevar)]
  freq      <- data.table::dcast(
    freq_long,
    stats::as.formula(paste(groupvar, "~", timevar)),
    value.var = "count",
    fill      = 0L
  )
  yr_cols  <- as.character(yrs)
  frac_mat <- as.matrix(freq[, ..yr_cols])
  Ti       <- rowSums(frac_mat)
  safe_Ti  <- pmax(Ti, 1L)   # T_i >= 1 by construction; guard is cosmetic
  frac_mat <- frac_mat / safe_Ti

  fbar_dt              <- data.table::as.data.table(frac_mat)
  fbar_dt[[groupvar]]  <- freq[[groupvar]]

  # --- 6. Column names: fbar_<clean_timevar>_<year> ---
  time_clean <- janitor::make_clean_names(timevar)
  cn         <- paste0(prefix, time_clean, "_", yrs)
  cn         <- make.unique(cn, sep = "_")
  data.table::setnames(fbar_dt, yr_cols, cn)

  # --- 7. Broadcast constants to all rows; zero out forecast rows.
  #        Units absent from cc_train receive NA from the join -> coerced to 0.
  data[fbar_dt, (cn) := mget(paste0("i.", cn)), on = groupvar]
  for (col in cn) {
    data[is.na(get(col)), (col) := 0]
  }
  data[data[[timevar]] >= test_start, (cn) := 0]

  # --- 8. Drop reference level (earliest year = cn[1]) to remove sum-to-1
  #        collinearity with the intercept, matching factor()'s default.
  kept    <- if (drop_ref) cn[-1L] else cn
  dropped <- setdiff(cn, kept)
  if (length(dropped) > 0L) data[, (dropped) := NULL]

  # Balanced-panel advisory: all fbar columns constant across units -> aliased.
  if (length(kept) > 0L) {
    all_constant <- all(vapply(kept, function(col) {
      vals <- data[[col]][data[[timevar]] < test_start]
      data.table::uniqueN(vals) <= 1L
    }, logical(1L)))
    if (all_constant) {
      message(
        "All `fbar_*` columns are constant across units: the panel appears ",
        "balanced. `lm()` treats aliased coefficients as NA (equivalent to 0 ",
        "in prediction); this is harmless."
      )
    }
  }

  # --- 9. Build exogen spec to carry fbar_* into the simulation grid ---
  spec <- build_model("exogen", formula = stats::reformulate(kept))

  # --- 10. Build extended model formula ---
  full <- stats::update(
    formula,
    stats::as.formula(paste("~ . +", paste(kept, collapse = " + ")))
  )

  list(data = data, terms = kept, spec = spec, formula = full)
}
