# Scenario parameter framework -----------------------------------------------
#
# Precomputed, model-keyed scenario parameters with fixed-effect heuristic
# helpers. All stochastic parameters are baked at build time via setup_param();
# simulate_system() reads the baked values instead of drawing lazily.

# --- Internal: detect factor(timevar) in a fitted linear model ---------------

# Scan fit_formula for a term of the form factor(<timevar>). Extract the
# estimated year effects from fitted_lm and return metadata, or NULL.
# `vcov` (optional) is the fit's coefficient covariance; when supplied, the
# returned `vcov` is the L x L covariance of the year effects in xlevels order
# (reference row/column and aliased entries = 0).
.detect_time_fe <- function(fit_formula, coefs, xlevels, timevar, vcov = NULL) {
  term_labels <- attr(stats::terms(fit_formula), "term.labels")

  fe_label <- NULL
  for (lbl in term_labels) {
    expr <- tryCatch(str2lang(lbl), error = function(e) NULL)
    if (is.null(expr) || !is.call(expr)) next
    fn <- as.character(expr[[1L]])
    if (!fn %in% c("factor", "stats::factor")) next
    if (length(expr) < 2L) next
    arg1 <- expr[[2L]]
    if (!is.symbol(arg1)) next
    if (as.character(arg1) != timevar) next
    fe_label <- lbl
    break
  }

  if (is.null(fe_label)) return(NULL)

  # Guard: timevar must not appear in any other RHS term
  other_labels <- setdiff(term_labels, fe_label)
  for (lbl in other_labels) {
    vars_in <- tryCatch(all.vars(str2lang(lbl)), error = function(e) character(0L))
    if (timevar %in% vars_in) {
      stop(
        "time fixed-effect scenario handling requires the time variable to enter ",
        "only through `factor(", timevar, ")` ",
        "(found it used elsewhere in the formula).",
        call. = FALSE
      )
    }
  }

  # xlevels for this term — present only when the model has been fitted
  xl <- xlevels[[fe_label]]
  if (is.null(xl)) return(NULL)

  # Build effects vector: baseline (xl[1]) -> 0, remaining levels -> coefficient.
  # Effects are named by level so downstream helpers can look them up by name.
  effects <- stats::setNames(numeric(length(xl)), xl)
  for (lv in xl[-1L]) {
    nm  <- paste0(fe_label, lv)
    val <- coefs[[nm]]
    if (!is.null(val)) effects[[lv]] <- val
  }
  ref_value <- as.numeric(xl[1L])

  V <- NULL
  if (!is.null(vcov)) {
    V    <- matrix(0, length(xl), length(xl), dimnames = list(xl, xl))
    nms  <- paste0(fe_label, xl)
    keep <- which(nms %in% rownames(vcov))
    keep <- keep[keep > 1L]
    if (length(keep) > 0L) {
      V[keep, keep] <- as.matrix(vcov)[nms[keep], nms[keep], drop = FALSE]
    }
    V[!is.finite(V)] <- 0
  }

  list(
    term      = fe_label,
    timevar   = timevar,
    ref_value = ref_value,
    effects   = effects,  # named by time level
    vcov      = V         # L x L over levels, or NULL
  )
}

# --- Internal: detect factor(unitkey) in a fitted linear model ---------------
#
# Scans fit_formula for a term of the form factor(<key>) where <key> is any
# element of unit_keys. Returns unit-FE metadata or NULL. Only the first
# matching key is used (single-key unit FE). Multi-key composite unit FE is
# out of scope: if unit_keys has length > 1, the factor over any one key is
# detected and the rest are left in the design unchanged.
.detect_unit_fe <- function(fit_formula, coefs, xlevels, unit_keys) {
  if (length(unit_keys) == 0L) return(NULL)
  term_labels <- attr(stats::terms(fit_formula), "term.labels")

  fe_label <- NULL
  fe_key   <- NULL
  for (lbl in term_labels) {
    expr <- tryCatch(str2lang(lbl), error = function(e) NULL)
    if (is.null(expr) || !is.call(expr)) next
    fn <- as.character(expr[[1L]])
    if (!fn %in% c("factor", "stats::factor")) next
    if (length(expr) < 2L) next
    arg1 <- expr[[2L]]
    if (!is.symbol(arg1)) next
    key_name <- as.character(arg1)
    if (!key_name %in% unit_keys) next
    fe_label <- lbl
    fe_key   <- key_name
    break
  }

  if (is.null(fe_label)) return(NULL)

  # Guard: unit key must not appear in any other RHS term
  other_labels <- setdiff(term_labels, fe_label)
  for (lbl in other_labels) {
    vars_in <- tryCatch(all.vars(str2lang(lbl)), error = function(e) character(0L))
    if (fe_key %in% vars_in) {
      stop(
        "unit fixed-effect scenario handling requires the unit variable to enter ",
        "only through `factor(", fe_key, ")` ",
        "(found it used elsewhere in the formula).",
        call. = FALSE
      )
    }
  }

  xl <- xlevels[[fe_label]]
  if (is.null(xl)) return(NULL)

  # Build effects named by unit level; baseline (xl[1]) = 0.
  effects <- stats::setNames(numeric(length(xl)), xl)
  for (lv in xl[-1L]) {
    nm  <- paste0(fe_label, lv)
    val <- coefs[[nm]]
    if (!is.null(val)) effects[[lv]] <- val
  }
  ref_value <- xl[1L]  # character unit id (reference level)

  list(
    term      = fe_label,
    unitvar   = fe_key,
    ref_value = ref_value,
    effects   = effects  # named by unit level
  )
}

# --- S3 generic: scenario_terms ---------------------------------------------

#' Retrieve fixed-effect metadata for scenario construction
#'
#' Returns a list with elements `time_fe` and `unit_fe`, each either a
#' metadata list (produced at fit time by `.detect_time_fe` /
#' `.detect_unit_fe`) or `NULL` when the corresponding fixed effect is absent
#' from the formula.
#'
#' This is the per-model-type extension seam consumed by [setup_param()].
#' Currently only `linear` models carry FE metadata; all other types return
#' `NULL` for both elements via the default method.  To add support for a new
#' model type, implement `scenario_terms.<class>` returning the same list
#' shape.
#'
#' @param model A fitted endogenr model object.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{`time_fe`}{Named list with `term`, `timevar`, `ref_value`,
#'     `effects` (named numeric; baseline level = 0), or `NULL`.}
#'   \item{`unit_fe`}{Named list with `term`, `unitvar`, `ref_value`,
#'     `effects` (named numeric; reference unit = 0), or `NULL`.}
#' }
#'
#' @family simulation
#' @export
#'
#' @examples
#' \dontrun{
#' dt  <- sim_panel_common_shock(units = 6L, n_time = 20L, seed = 1L)
#' sys <- setup_system(
#'   list(build_model("linear", formula = y ~ lag(y) + factor(time))),
#'   data = dt, train_start = 1, test_start = 16, horizon = 4,
#'   groupvar = "unit", timevar = "time", inner_sims = 2L
#' )
#' fit <- fit_system(sys, nsim = 2L)
#'
#' # Access time-FE metadata from the representative fitted model
#' model <- fit$fitted_models[[1L]]
#' st    <- scenario_terms(model)
#' st$time_fe$term      # "factor(time)"
#' st$time_fe$effects   # named numeric: year -> estimated offset
#' st$unit_fe           # NULL (no factor(unit) in formula)
#' }
scenario_terms <- function(model) UseMethod("scenario_terms")

#' @rdname scenario_terms
#' @exportS3Method scenario_terms default
scenario_terms.default <- function(model) list(time_fe = NULL, unit_fe = NULL)

#' @rdname scenario_terms
#' @exportS3Method scenario_terms linear
scenario_terms.linear <- function(model) {
  list(time_fe = model$time_fe, unit_fe = model$unit_fe)
}

#' @rdname scenario_terms
#' @exportS3Method scenario_terms endogenr_gamlss
scenario_terms.endogenr_gamlss <- function(model) {
  list(time_fe = model$time_fe, unit_fe = model$unit_fe)
}

# --- Internal: heuristic engine helpers -------------------------------------

# Compute a weight vector in [0, 1] for horizon steps.
# "constant" -> all zeros (no movement toward target / along factor).
# "linear"   -> seq(1/H, 1, 1/H); reaches 1 at the final step.
# "sigmoid"  -> S-curve rescaled to [0, 1], 0 at step 1, 1 at step H.
.fe_weights <- function(path, horizon, mid = NULL, steep = 1) {
  switch(path,
    constant = rep(0, horizon),
    linear   = seq_len(horizon) / horizon,
    sigmoid  = {
      if (is.null(mid)) mid <- (horizon + 1L) / 2
      raw <- stats::plogis((seq_len(horizon) - mid) * steep)
      if (horizon > 1L) {
        (raw - raw[1L]) / (raw[horizon] - raw[1L])
      } else {
        1
      }
    },
    stop("unknown path '", path, "'; must be 'constant', 'linear', or 'sigmoid'",
         call. = FALSE)
  )
}

# Resolve the 'to' argument to a scalar target.
# Accepts: numeric scalar, "zero", "mean", or character unit id(s) in effects.
.resolve_target <- function(to, effects) {
  if (is.numeric(to))          return(as.numeric(to))
  if (identical(to, "zero"))   return(0)
  if (identical(to, "mean"))   return(mean(effects))
  if (is.character(to)) {
    nm    <- names(effects)
    found <- to[to %in% nm]
    if (length(found) == 0L) {
      stop("'to' unit id(s) not found in effect names: ",
           paste(to, collapse = ", "), call. = FALSE)
    }
    return(mean(effects[found]))
  }
  stop("'to' must be numeric, \"zero\", \"mean\", or a character unit id vector",
       call. = FALSE)
}

# n x L matrix of perturbed year-effect series for outer draw i (columns in the
# draw's own level order = chronological xlevels order). One MVN draw per row.
.fe_time_perturbed <- function(block, i, n) {
  mu <- block$effects_by_draw[[i]]
  V  <- if (!is.null(block$vcov_by_draw)) block$vcov_by_draw[[i]] else NULL
  if (is.null(V)) {
    .warn_once(".endogenr_time_fe_no_vcov",
      "Time fixed-effect covariance unavailable; time-FE draws use point estimates without parameter uncertainty.")
    return(matrix(mu, nrow = n, ncol = length(mu), byrow = TRUE,
                  dimnames = list(NULL, names(mu))))
  }
  out <- .cf_rmv(n, mu, V[names(mu), names(mu), drop = FALSE])
  colnames(out) <- names(mu)
  out
}

# --- Exported heuristic helpers: FE blocks -----------------------------------

#' Fixed-effect heuristic: resample from estimated effects
#'
#' Re-bakes an FE block's `values` array by sampling with replacement from the
#' per-draw pool of estimated effects, drawing independently for every
#' `(outer draw, inner sim, forecast step)` combination.  [setup_param()]
#' falls back to this strategy for time fixed effects when fewer than 3 year
#' effects are estimated (the default is [fe_ar()]).
#'
#' @section Parameter uncertainty (time FE):
#' For time fixed effects, each `(outer draw, inner sim)` trajectory first
#' draws one perturbed year-effect series from `MVN(tau_hat, V_tau)` (the
#' fitted estimates and their covariance); the trajectory's resampling pool is
#' that perturbed series.  When the covariance is unavailable, the point
#' estimates are used (with a one-time warning).
#'
#' @section Baked array shapes:
#' \describe{
#'   \item{Time FE (`kind = "time"`)}{`dim = c(nsim, inner_sims, horizon)`.
#'     Each cell is one independent resample from the trajectory's perturbed
#'     year effects.}
#'   \item{Unit FE (`kind = "unit"`)}{`dim = c(n_units, nsim, inner_sims,
#'     horizon)`, `attr(values, "per_trajectory") = TRUE`.  Each cell is an
#'     independent draw from the pooled effects, assigned per unit × trajectory
#'     combination.  The full 4-D array can be memory-intensive for large panels;
#'     consider [fe_persist()] or [fe_converge()] for deterministic alternatives.}
#' }
#'
#' @section Path and target mechanics:
#' All FE helpers share a common framework for gradually moving the baseline
#' toward a target or scaling it over the forecast horizon.  A weight vector
#' `w[h]` in `[0, 1]` is computed from `path`:
#' \describe{
#'   \item{`"constant"` (default)}{`w[h] = 0` for all steps — the baseline is
#'     returned unchanged.}
#'   \item{`"linear"`}{`w[h] = h / H`; reaches 1 at the final step.}
#'   \item{`"sigmoid"`}{S-curve anchored at 0 and 1, inflecting at step `mid`
#'     (default midpoint) with steepness `steep`.}
#' }
#' The weight is applied differently depending on whether `to` is set:
#' \itemize{
#'   \item \strong{`to` is non-`NULL`}: offset at step `h` =
#'     `(1 - w[h]) * baseline + w[h] * target`.
#'   \item \strong{`to = NULL`, `factor != 1`}: offset at step `h` =
#'     `baseline * (1 + (factor - 1) * w[h])`.  A `factor > 1` diverges
#'     (amplifies effects over time); `factor < 1` shrinks them.
#'   \item \strong{All defaults} (`to = NULL`, `factor = 1`,
#'     `path = "constant"`): `w = 0` everywhere → offsets equal the resampled
#'     baseline at every step.
#' }
#'
#' @param block An `endogenr_fe_param` from [setup_param()].
#' @param to Convergence target: `NULL` (no convergence), a numeric scalar,
#'   `"zero"`, `"mean"` (mean of all estimated effects), or a character vector
#'   of unit ids whose mean effect is used.
#' @param factor Multiplicative scale applied along `path` when `to = NULL`
#'   (default `1` = no scaling; `> 1` diverge; `< 1` shrink).
#' @param path Shape of the movement toward `to` or along `factor`:
#'   `"constant"` (default, no movement), `"linear"`, or `"sigmoid"`.
#' @param mid Sigmoid inflection step index (default: midpoint of the horizon).
#' @param steep Sigmoid steepness parameter (default `1`; larger = sharper
#'   transition).
#'
#' @return The same `endogenr_fe_param` block with `values`, `active = TRUE`,
#'   and `heuristic` updated.  Replace the block in the `setup_param()` result
#'   and pass to [simulate_system()].
#' @seealso [fe_distribution()], [fe_fixed()], [fe_persist()], [fe_converge()],
#'   [setup_param()], [simulate_system()]
#' @family simulation
#' @export
#'
#' @examples
#' \dontrun{
#' # -- Shared setup (reused in all FE-helper examples) ----------------------
#' dt  <- sim_panel_common_shock(units = 8L, n_time = 20L, seed = 1L)
#' sys <- setup_system(
#'   list(build_model("linear", formula = y ~ lag(y) + factor(time))),
#'   data = dt, train_start = 1, test_start = 16, horizon = 4,
#'   groupvar = "unit", timevar = "time", inner_sims = 3L
#' )
#' fit <- fit_system(sys, nsim = 5L)
#' set.seed(1); sp <- setup_param(fit)  # default time-FE sampler is fe_ar()
#' blk <- sp$y$time_fe
#' dim(blk$values)       # c(5, 3, 4) = c(nsim, inner_sims, horizon)
#'
#' # Re-bake with a fresh seed (e.g. to change RNG state)
#' set.seed(42)
#' sp$y$time_fe <- fe_resample(blk)
#' res <- simulate_system(fit, scenario_params = sp)
#'
#' # Resample baseline that converges linearly toward 0 over the horizon
#' set.seed(42)
#' sp$y$time_fe <- fe_resample(blk, to = 0, path = "linear")
#' # At step 1 the offset is ~the resampled value; at step 4 it is 0.
#' round(sp$y$time_fe$values[1, 1, ], 3)   # monotonically shrinking
#'
#' # Amplify (diverge) the time-FE variance over the horizon by factor 2
#' set.seed(42)
#' sp$y$time_fe <- fe_resample(blk, factor = 2, path = "linear")
#' }
fe_resample <- function(block, to = NULL, factor = 1, path = "constant",
                        mid = NULL, steep = 1) {
  dims       <- block$dims
  nsim       <- dims$nsim
  inner_sims <- dims$inner_sims
  horizon    <- dims$horizon
  w          <- .fe_weights(path, horizon, mid, steep)
  t_val      <- if (!is.null(to)) .resolve_target(to, block$effects) else NULL

  if (block$kind == "time") {
    vals <- array(NA_real_, dim = c(nsim, inner_sims, horizon))
    for (i in seq_len(nsim)) {
      Tt <- .fe_time_perturbed(block, i, inner_sims)
      # Index columns via sample.int: sample(x) on a length-1 pool means 1:x.
      B  <- matrix(Tt[cbind(rep(seq_len(inner_sims), horizon),
                            sample.int(ncol(Tt), inner_sims * horizon,
                                       replace = TRUE))],
                   nrow = inner_sims, ncol = horizon)
      for (h in seq_len(horizon)) {
        vals[i,, h] <- if (!is.null(t_val))
          (1 - w[h]) * B[, h] + w[h] * t_val
        else
          B[, h] * (1 + (factor - 1) * w[h])
      }
    }
  } else {
    # Unit FE: c(n_units, nsim, inner_sims, horizon) — per_trajectory = TRUE.
    n_units <- length(dims$units)
    vals    <- array(NA_real_, dim = c(n_units, nsim, inner_sims, horizon))
    dimnames(vals)[[1L]] <- as.character(dims$units)
    wfull <- array(rep(w, each = n_units * inner_sims),
                   dim = c(n_units, inner_sims, horizon))
    for (i in seq_len(nsim)) {
      pool  <- block$effects_by_draw[[i]]
      draws <- sample(pool, n_units * inner_sims * horizon, replace = TRUE)
      P <- aperm(array(draws, dim = c(inner_sims, horizon, n_units)), c(3L, 1L, 2L))
      vals[, i, , ] <- if (!is.null(t_val)) (1 - wfull) * P + wfull * t_val
                       else                 P * (1 + (factor - 1) * wfull)
    }
    attr(vals, "per_trajectory") <- TRUE
  }

  block$values    <- vals
  block$active    <- TRUE
  block$heuristic <- list(strategy = "resample", to = to, factor = factor,
                          path = path, mid = mid, steep = steep)
  block
}

#' Time fixed-effect heuristic: AR(1) forecast from the origin
#'
#' Re-bakes a time-FE block's `values` by treating the estimated year effects
#' as a time series and forecasting it with an AR(1) that starts at the last
#' training-year effect.  This is the **default strategy applied to time fixed
#' effects** by [setup_param()] (when at least 3 year effects are estimated).
#'
#' @section Algorithm:
#' For each `(outer draw, inner sim)` trajectory:
#' \enumerate{
#'   \item Draw one perturbed year-effect series
#'     `tau_tilde ~ MVN(tau_hat, V_tau)` from the fitted estimates and their
#'     covariance (point estimates, with a one-time warning, when the
#'     covariance is unavailable).
#'   \item Fit `tau_t = c + rho * tau_{t-1} + e_t` by OLS on that draw;
#'     `rho` is clamped to `[-1, 1]` (a unit root, i.e. random walk with
#'     drift, is the most persistent case).
#'   \item Simulate forward from the last training-year effect, drawing the
#'     innovations by bootstrap from the centred AR residuals.
#'   \item If the training window ends before the forecast origin
#'     (`test_start - 1`, e.g. random-window refits), the recursion is first
#'     stepped through the gap years and only the forecast steps are kept.
#' }
#' The forecast offsets are then moved along `path` toward `to` or scaled by
#' `factor` exactly as in [fe_resample()].  The fitted AR parameters are
#' stored as `nsim x inner_sims` matrices in `heuristic$rho` and
#' `heuristic$intercept`.
#'
#' @inheritParams fe_resample
#' @param block A time-FE `endogenr_fe_param` (`kind = "time"`) from
#'   [setup_param()].
#'
#' @return The same `endogenr_fe_param` block with `values`
#'   (`c(nsim, inner_sims, horizon)`), `active = TRUE`, and `heuristic`
#'   updated.
#' @seealso [fe_resample()], [fe_distribution()], [fe_converge()],
#'   [setup_param()], [simulate_system()]
#' @family simulation
#' @export
#'
#' @examples
#' \dontrun{
#' # (Reuse 'fit' and 'sp' from the fe_resample() example)
#' blk <- sp$y$time_fe
#' set.seed(42)
#' sp$y$time_fe <- fe_ar(blk)
#' summary(as.vector(sp$y$time_fe$heuristic$rho))
#'
#' # AR(1) baseline that converges linearly toward 0 over the horizon
#' sp$y$time_fe <- fe_ar(blk, to = 0, path = "linear")
#' res <- simulate_system(fit, scenario_params = sp)
#' }
fe_ar <- function(block, to = NULL, factor = 1, path = "constant",
                  mid = NULL, steep = 1) {
  if (!identical(block$kind, "time")) {
    stop("fe_ar() supports time fixed effects only.", call. = FALSE)
  }
  if (any(lengths(block$effects_by_draw) < 3L)) {
    stop("fe_ar() needs at least 3 estimated time effects per draw; use fe_resample().",
         call. = FALSE)
  }
  dims       <- block$dims
  nsim       <- dims$nsim
  inner_sims <- dims$inner_sims
  horizon    <- dims$horizon
  w          <- .fe_weights(path, horizon, mid, steep)
  t_val      <- if (!is.null(to)) .resolve_target(to, block$effects) else NULL
  origin     <- dims$test_start - 1L

  vals  <- array(NA_real_, dim = c(nsim, inner_sims, horizon))
  rho_m <- matrix(NA_real_, nrow = nsim, ncol = inner_sims)
  int_m <- matrix(NA_real_, nrow = nsim, ncol = inner_sims)
  rows  <- seq_len(inner_sims)

  for (i in seq_len(nsim)) {
    Tt <- .fe_time_perturbed(block, i, inner_sims)
    L  <- ncol(Tt)

    # Years between the last estimated effect and the forecast origin.
    last <- suppressWarnings(as.numeric(utils::tail(colnames(Tt), 1L)))
    gap  <- origin - last
    if (!is.finite(gap) || gap < 0) gap <- 0L

    # AR(1) by OLS per trajectory; keep centred residuals for the bootstrap.
    E <- matrix(NA_real_, nrow = inner_sims, ncol = L - 1L)
    for (s in rows) {
      x   <- Tt[s, -L]
      y   <- Tt[s, -1L]
      vx  <- stats::var(x)
      rho <- if (vx > 0) stats::cov(x, y) / vx else 0
      rho <- min(max(rho, -1), 1)
      c0  <- mean(y - rho * x)
      e   <- y - c0 - rho * x
      E[s, ]      <- e - mean(e)
      rho_m[i, s] <- rho
      int_m[i, s] <- c0
    }

    B    <- matrix(NA_real_, nrow = inner_sims, ncol = horizon)
    prev <- Tt[, L]
    for (k in seq_len(gap + horizon)) {
      prev <- int_m[i, ] + rho_m[i, ] * prev +
        E[cbind(rows, sample.int(L - 1L, inner_sims, replace = TRUE))]
      if (k > gap) B[, k - gap] <- prev
    }

    for (h in seq_len(horizon)) {
      vals[i,, h] <- if (!is.null(t_val))
        (1 - w[h]) * B[, h] + w[h] * t_val
      else
        B[, h] * (1 + (factor - 1) * w[h])
    }
  }

  block$values    <- vals
  block$active    <- TRUE
  block$heuristic <- list(strategy = "ar1", to = to, factor = factor,
                          path = path, mid = mid, steep = steep,
                          rho = rho_m, intercept = int_m)
  block
}

#' Fixed-effect heuristic: normal distribution
#'
#' Re-bakes an FE block's `values` by drawing i.i.d. Normal offsets,
#' independently for every `(outer draw, inner sim, forecast step)`.  The
#' distribution parameters default to the per-draw mean and sd of the estimated
#' effects, so the spread is calibrated to the model's own historical variation.
#' Fix `mean` and/or `sd` to impose a specific distribution.
#'
#' @section When to use instead of fe_resample():
#' [fe_resample()] draws directly from the discrete pool of estimated year
#' effects; [fe_distribution()] draws from a fitted Normal, which smooths out
#' discreteness and allows the support to extend beyond the observed range.  Use
#' `fe_distribution()` when you want to widen or narrow the spread (via `sd`)
#' or shift the center (via `mean`) relative to historical estimates, or when
#' the pool of estimated levels is too small to resample from reliably.
#'
#' @param block An `endogenr_fe_param` from [setup_param()].
#' @param mean Scalar mean of the Normal distribution (default: per-draw mean
#'   of `effects_by_draw`, so it adapts to each bootstrap draw).
#' @param sd Scalar standard deviation of the Normal distribution (default:
#'   per-draw sd of `effects_by_draw`).  Set to a small value (e.g. `0`) to
#'   produce near-constant offsets while keeping the Normal draw machinery.
#' @inheritParams fe_resample
#'
#' @return The same `endogenr_fe_param` block with `values`, `active = TRUE`,
#'   and `heuristic` updated.
#' @seealso [fe_resample()], [fe_fixed()], [fe_converge()], [setup_param()]
#' @family simulation
#' @export
#'
#' @examples
#' \dontrun{
#' # (Reuse 'fit' and 'sp' from the fe_resample() example)
#' blk <- sp$y$time_fe
#'
#' # Default: Normal with per-draw mean and sd from estimated year effects
#' set.seed(42)
#' sp$y$time_fe <- fe_distribution(blk)
#' res_norm <- simulate_system(fit, scenario_params = sp)
#'
#' # Narrow the spread to 1/10 of the estimated sd (near-deterministic mean)
#' set.seed(42)
#' sp$y$time_fe <- fe_distribution(blk, sd = 0.01)
#'
#' # Fix both mean and sd: offsets drawn from N(0, 0.5) for all trajectories
#' set.seed(42)
#' sp$y$time_fe <- fe_distribution(blk, mean = 0, sd = 0.5)
#'
#' # Combine with path: Normal baseline converging to zero over the horizon
#' set.seed(42)
#' sp$y$time_fe <- fe_distribution(blk, sd = 0.5, to = 0, path = "linear")
#' }
fe_distribution <- function(block, mean = NULL, sd = NULL, to = NULL,
                             factor = 1, path = "constant", mid = NULL, steep = 1) {
  dims       <- block$dims
  nsim       <- dims$nsim
  inner_sims <- dims$inner_sims
  horizon    <- dims$horizon
  w          <- .fe_weights(path, horizon, mid, steep)
  t_val      <- if (!is.null(to)) .resolve_target(to, block$effects) else NULL

  if (block$kind == "time") {
    vals <- array(NA_real_, dim = c(nsim, inner_sims, horizon))
    for (i in seq_len(nsim)) {
      eff_i <- block$effects_by_draw[[i]]
      m_i   <- if (!is.null(mean)) mean else base::mean(eff_i)
      s_i   <- if (!is.null(sd))   sd   else stats::sd(eff_i)
      B     <- matrix(stats::rnorm(inner_sims * horizon, m_i, s_i),
                      nrow = inner_sims, ncol = horizon)
      for (h in seq_len(horizon)) {
        vals[i,, h] <- if (!is.null(t_val))
          (1 - w[h]) * B[, h] + w[h] * t_val
        else
          B[, h] * (1 + (factor - 1) * w[h])
      }
    }
  } else {
    n_units <- length(dims$units)
    vals    <- array(NA_real_, dim = c(n_units, nsim, inner_sims, horizon))
    dimnames(vals)[[1L]] <- as.character(dims$units)
    wfull <- array(rep(w, each = n_units * inner_sims),
                   dim = c(n_units, inner_sims, horizon))
    for (i in seq_len(nsim)) {
      eff_i <- block$effects_by_draw[[i]]
      m_i   <- if (!is.null(mean)) mean else base::mean(eff_i)
      s_i   <- if (!is.null(sd))   sd   else stats::sd(eff_i)
      draws <- stats::rnorm(n_units * inner_sims * horizon, m_i, s_i)
      P <- aperm(array(draws, dim = c(inner_sims, horizon, n_units)), c(3L, 1L, 2L))
      vals[, i, , ] <- if (!is.null(t_val)) (1 - wfull) * P + wfull * t_val
                       else                 P * (1 + (factor - 1) * wfull)
    }
    attr(vals, "per_trajectory") <- TRUE
  }

  block$values    <- vals
  block$active    <- TRUE
  block$heuristic <- list(strategy = "distribution", to = to, factor = factor,
                          path = path, mid = mid, steep = steep,
                          mean = mean, sd = sd)
  block
}

#' Fixed-effect heuristic: fixed constant
#'
#' Re-bakes an FE block's `values` with a single deterministic constant.
#' Every trajectory receives the same offset at every step (unless combined with
#' `to` or `factor` + `path` to ramp it).  This completely removes
#' between-trajectory variance from the FE channel, making it useful as a
#' controlled baseline for sensitivity analysis.
#'
#' The `values` array is always shaped `c(nsim, inner_sims, horizon)` for time
#' FE and `c(n_units, horizon)` for unit FE, preserving the same consumer
#' contract as the stochastic strategies.
#'
#' @section Common uses:
#' \describe{
#'   \item{Null scenario}{`fe_fixed(blk, value = 0)` removes the time-FE
#'     contribution entirely, isolating the effect of other predictors.}
#'   \item{Mean scenario}{`fe_fixed(blk, value = mean(blk$effects))` fixes the
#'     offset at the historical cross-year average for all trajectories.}
#'   \item{Ramp to zero}{`fe_fixed(blk, value = mean(blk$effects), to = 0,
#'     path = "linear")` starts at the historical mean and linearly reaches 0 by
#'     the final forecast step.}
#' }
#'
#' @param block An `endogenr_fe_param` from [setup_param()].
#' @param value Numeric scalar offset applied to all trajectories and steps
#'   (default `0`).
#' @inheritParams fe_resample
#'
#' @return The same `endogenr_fe_param` block with `values`, `active = TRUE`,
#'   and `heuristic` updated.
#' @seealso [fe_resample()], [fe_distribution()], [fe_converge()], [setup_param()]
#' @family simulation
#' @export
#'
#' @examples
#' \dontrun{
#' # (Reuse 'fit' and 'sp' from the fe_resample() example)
#' blk <- sp$y$time_fe
#'
#' # Zero out the time-FE offset for all trajectories (null scenario)
#' sp$y$time_fe <- fe_fixed(blk, value = 0)
#' res_null <- simulate_system(fit, scenario_params = sp)
#'
#' # Fix at the historical mean offset (single deterministic trajectory)
#' sp$y$time_fe <- fe_fixed(blk, value = mean(blk$effects))
#'
#' # Start at the mean and linearly ramp toward 0 (e.g. "fading shock")
#' sp$y$time_fe <- fe_fixed(blk, value = mean(blk$effects),
#'                           to = 0, path = "linear")
#' round(sp$y$time_fe$values[1, 1, ], 3)   # decreasing sequence
#'
#' # Compare ensemble variance: resample vs fixed zero
#' set.seed(1); sp_rs  <- setup_param(fit)   # default fe_ar
#' sp_fix <- setup_param(fit)
#' sp_fix$y$time_fe <- fe_fixed(sp_fix$y$time_fe, value = 0)
#' set.seed(1); res_rs  <- simulate_system(fit)
#' set.seed(1); res_fix <- simulate_system(fit, scenario_params = sp_fix)
#' # res_rs has larger cross-trajectory variance in forecast years
#' }
fe_fixed <- function(block, value = 0, to = NULL, factor = 1, path = "constant",
                     mid = NULL, steep = 1) {
  dims       <- block$dims
  nsim       <- dims$nsim
  inner_sims <- dims$inner_sims
  horizon    <- dims$horizon
  w          <- .fe_weights(path, horizon, mid, steep)
  t_val      <- if (!is.null(to)) .resolve_target(to, block$effects) else NULL

  if (block$kind == "time") {
    vals <- array(NA_real_, dim = c(nsim, inner_sims, horizon))
    for (h in seq_len(horizon)) {
      b_h <- if (!is.null(t_val)) (1 - w[h]) * value + w[h] * t_val
             else                  value * (1 + (factor - 1) * w[h])
      vals[,, h] <- b_h
    }
  } else {
    n_units <- length(dims$units)
    vals    <- matrix(NA_real_, nrow = n_units, ncol = horizon)
    rownames(vals) <- as.character(dims$units)
    for (h in seq_len(horizon)) {
      b_h <- if (!is.null(t_val)) (1 - w[h]) * value + w[h] * t_val
             else                  value * (1 + (factor - 1) * w[h])
      vals[, h] <- b_h
    }
    attr(vals, "per_trajectory") <- FALSE
  }

  block$values    <- vals
  block$active    <- TRUE
  block$heuristic <- list(strategy = "fixed", value = value, to = to,
                          factor = factor, path = path, mid = mid, steep = steep)
  block
}

#' Fixed-effect heuristic: persist historical unit effects (default)
#'
#' The **default strategy for unit fixed effects** applied by [setup_param()].
#' Each simulation unit's offset is set to its historically estimated effect
#' and held constant across the horizon.
#'
#' With all defaults (`to = NULL`, `factor = 1`, `path = "constant"`),
#' `active = FALSE` is set on the block, which means `factor(unit)` stays in
#' the design matrix **completely unchanged** — unit dummies contribute their
#' normal contrast at predict time and no explicit offset array is added.  This
#' is identical to running without any scenario parameter at all.
#'
#' Pass any non-default `to`, `factor`, or `path` to activate the block
#' (`active = TRUE`) and apply a path-shaped adjustment on top of each unit's
#' historical baseline.  For a simple convergence-to-zero use [fe_converge()].
#'
#' @section Active flag and design matrix:
#' When `active = FALSE` (the default), `predict.linear` keeps `factor(unit)`
#' live in the design matrix and adds nothing extra.  When `active = TRUE`, the
#' unit column is neutralised to the reference level in the prediction frame
#' (zeroing out the factor contrast) and the baked `values` matrix — shaped
#' `c(n_units, horizon)` — is added row-by-row.  This lets you override or
#' fade the cross-sectional intercepts explicitly.
#'
#' @param block An `endogenr_fe_param` with `kind == "unit"` from
#'   [setup_param()].
#' @inheritParams fe_resample
#'
#' @return The same `endogenr_fe_param` block with `values` (`c(n_units,
#'   horizon)`, `per_trajectory = FALSE`), `active` (`FALSE` under all defaults,
#'   `TRUE` otherwise), and `heuristic` updated.
#' @seealso [fe_converge()], [fe_fixed()], [setup_param()]
#' @family simulation
#' @export
#'
#' @examples
#' \dontrun{
#' # Build a unit-FE model
#' dt  <- sim_panel_ar1(units = 8L, n_time = 20L, seed = 1L)
#' sys <- setup_system(
#'   list(build_model("exogen", formula = ~x),
#'        build_model("linear", formula = y ~ lag(y) + x + factor(unit))),
#'   data = dt, train_start = 1, test_start = 16, horizon = 4,
#'   groupvar = "unit", timevar = "time", inner_sims = 2L
#' )
#' fit <- fit_system(sys, nsim = 4L)
#' sp  <- setup_param(fit)
#'
#' # Default: active = FALSE — factor(unit) in design, no offset override
#' sp$y$unit_fe$active   # FALSE
#' sp$y$unit_fe$values   # c(n_units, horizon) of historical effects (for reference)
#'
#' # Explicit persist (identity, same as default)
#' sp$y$unit_fe <- fe_persist(sp$y$unit_fe)
#' sp$y$unit_fe$active   # still FALSE
#'
#' # Shrink unit effects to 50% of historical values linearly across the horizon
#' sp$y$unit_fe <- fe_persist(sp$y$unit_fe, factor = 0.5, path = "linear")
#' sp$y$unit_fe$active   # TRUE
#' res <- simulate_system(fit, scenario_params = sp)
#' }
fe_persist <- function(block, to = NULL, factor = 1, path = "constant",
                       mid = NULL, steep = 1) {
  dims       <- block$dims
  n_units    <- length(dims$units)
  horizon    <- dims$horizon
  is_default <- is.null(to) && identical(factor, 1) && identical(path, "constant")

  w        <- .fe_weights(path, horizon, mid, steep)
  t_val    <- if (!is.null(to)) .resolve_target(to, block$effects) else NULL
  unit_nms <- as.character(dims$units)
  base_eff <- block$effects[unit_nms]   # named numeric, length n_units

  vals <- matrix(base_eff, nrow = n_units, ncol = horizon)
  rownames(vals) <- unit_nms

  if (!is_default) {
    for (h in seq_len(horizon)) {
      vals[, h] <- if (!is.null(t_val))
        (1 - w[h]) * base_eff + w[h] * t_val
      else
        base_eff * (1 + (factor - 1) * w[h])
    }
  }

  attr(vals, "per_trajectory") <- FALSE
  block$values    <- vals
  block$active    <- !is_default
  block$heuristic <- list(strategy = "persist", to = to, factor = factor,
                          path = path, mid = mid, steep = steep)
  block
}

#' Fixed-effect heuristic: converge toward a target
#'
#' Convenience wrapper that applies a path-shaped convergence from the block's
#' natural baseline toward `to`.  This is the recommended helper when you want
#' to neutralise fixed effects over the forecast horizon — e.g. letting unit
#' intercepts or time shocks fade toward zero, the cross-unit mean, or a
#' specific historical unit's level.
#'
#' @section Baseline per kind:
#' \describe{
#'   \item{Unit FE (`kind = "unit"`)}{Each unit's baseline is its own
#'     historically estimated effect (same as [fe_persist()]).  The convergence
#'     is deterministic: `values` is `c(n_units, horizon)`,
#'     `per_trajectory = FALSE`.  All units converge at the same rate but from
#'     different starting points, so cross-sectional variation in intercepts
#'     persists early in the horizon and compresses toward the target late.}
#'   \item{Time FE (`kind = "time"`)}{A single effect is drawn per
#'     `(outer draw, inner sim)` from that trajectory's pool — the year effects
#'     perturbed once by `MVN(tau_hat, V_tau)` (see [fe_resample()]) — and held
#'     constant across steps as the baseline (unlike [fe_resample()], which
#'     draws independently at every step).  This keeps each trajectory
#'     internally consistent while still varying across trajectories.  The
#'     drawn baseline then converges to `to` along `path`.}
#' }
#'
#' @section Target specification (`to`):
#' \describe{
#'   \item{Numeric scalar}{Converge to that specific value (e.g. `0`).}
#'   \item{`"zero"`}{Synonym for `to = 0`.}
#'   \item{`"mean"`}{Converge to the mean of all estimated effects
#'     (`mean(block$effects)`).}
#'   \item{Character vector of unit ids}{Converge to the mean effect of those
#'     specific units — useful for anchoring a counterfactual to a reference
#'     group's historical intercept level.}
#' }
#'
#' @param block An `endogenr_fe_param` from [setup_param()].
#' @param to Convergence target (required).  See Target specification above.
#' @param path Shape of the convergence: `"linear"` (default — reaches `to`
#'   exactly at the final forecast step), `"constant"` (jump immediately to
#'   `to` at every step — equivalent to `fe_fixed(value = to)`), or
#'   `"sigmoid"` (S-curve transition).
#' @param mid Sigmoid inflection step (default: midpoint of the horizon).
#' @param steep Sigmoid steepness (default `1`; larger = sharper transition).
#'
#' @return The same `endogenr_fe_param` block with `values`, `active = TRUE`,
#'   and `heuristic` updated.
#' @seealso [fe_persist()], [fe_resample()], [fe_fixed()], [setup_param()],
#'   [simulate_system()]
#' @family simulation
#' @export
#'
#' @examples
#' \dontrun{
#' # ── Time FE: convergence toward zero ──────────────────────────────────────
#' # (Reuse 'fit' and 'sp' from the fe_resample() example)
#' blk <- sp$y$time_fe
#'
#' # Linear convergence: offsets start at the per-trajectory draw and reach 0
#' # at the last forecast step
#' set.seed(42)
#' sp$y$time_fe <- fe_converge(blk, to = 0, path = "linear")
#' mean(abs(sp$y$time_fe$values[,, 1]))  # large
#' mean(abs(sp$y$time_fe$values[,, 4]))  # near 0
#'
#' # Sigmoid convergence (slow start, fast middle, slow end)
#' set.seed(42)
#' sp$y$time_fe <- fe_converge(blk, to = 0, path = "sigmoid", steep = 2)
#'
#' # Converge to the cross-year mean rather than zero
#' set.seed(42)
#' sp$y$time_fe <- fe_converge(blk, to = "mean", path = "linear")
#' res <- simulate_system(fit, scenario_params = sp)
#'
#' # ── Unit FE: fade intercepts toward zero ─────────────────────────────────
#' dt2  <- sim_panel_ar1(units = 8L, n_time = 20L, seed = 2L)
#' sys2 <- setup_system(
#'   list(build_model("exogen", formula = ~x),
#'        build_model("linear", formula = y ~ lag(y) + x + factor(unit))),
#'   data = dt2, train_start = 1, test_start = 16, horizon = 4,
#'   groupvar = "unit", timevar = "time", inner_sims = 2L
#' )
#' fit2 <- fit_system(sys2, nsim = 4L)
#' sp2  <- setup_param(fit2)
#'
#' # Converge each unit's historical intercept to zero over the horizon
#' sp2$y$unit_fe <- fe_converge(sp2$y$unit_fe, to = 0, path = "linear")
#' sp2$y$unit_fe$active  # TRUE
#'
#' # Converge toward the mean cross-sectional intercept
#' sp2$y$unit_fe <- fe_converge(sp2$y$unit_fe, to = "mean")
#'
#' # Converge toward the mean effect of specific reference units
#' ref_units <- as.character(unique(dt2$unit)[1:2])
#' sp2$y$unit_fe <- fe_converge(sp2$y$unit_fe, to = ref_units)
#' res2 <- simulate_system(fit2, scenario_params = sp2)
#' }
fe_converge <- function(block, to, path = "linear", mid = NULL, steep = 1) {
  dims    <- block$dims
  horizon <- dims$horizon
  w       <- .fe_weights(path, horizon, mid, steep)
  t_val   <- .resolve_target(to, block$effects)

  if (block$kind == "unit") {
    n_units  <- length(dims$units)
    unit_nms <- as.character(dims$units)
    base_eff <- block$effects[unit_nms]
    vals     <- matrix(NA_real_, nrow = n_units, ncol = horizon)
    rownames(vals) <- unit_nms
    for (h in seq_len(horizon)) {
      vals[, h] <- (1 - w[h]) * base_eff + w[h] * t_val
    }
    attr(vals, "per_trajectory") <- FALSE
    block$values    <- vals
    block$active    <- TRUE
    block$heuristic <- list(strategy = "converge", to = to, path = path,
                            mid = mid, steep = steep)
  } else {
    # Time FE: draw one value per (outer draw, inner_sim), hold across steps.
    nsim       <- dims$nsim
    inner_sims <- dims$inner_sims
    vals       <- array(NA_real_, dim = c(nsim, inner_sims, horizon))
    for (i in seq_len(nsim)) {
      Tt        <- .fe_time_perturbed(block, i, inner_sims)
      base_draw <- Tt[cbind(seq_len(inner_sims),
                            sample.int(ncol(Tt), inner_sims, replace = TRUE))]
      for (h in seq_len(horizon)) {
        vals[i,, h] <- (1 - w[h]) * base_draw + w[h] * t_val
      }
    }
    block$values    <- vals
    block$active    <- TRUE
    block$heuristic <- list(strategy = "converge", to = to, path = path,
                            mid = mid, steep = steep)
  }
  block
}

#' Apply a coefficient override to a coef block
#'
#' Bakes a per-draw, per-step coefficient override into the `coef` block of a
#' [setup_param()] result, replacing the fitted estimate for the named term
#' with a user-specified value.  The effective value at each forecast step and
#' outer draw is resolved immediately and stored in `block$effective[[term]]`
#' as an `nsim × horizon` matrix, so the override is inspectable before
#' [simulate_system()] is called.
#'
#' @section Override forms:
#' \describe{
#'   \item{Scalar (`value` is a single number)}{The same constant replaces
#'     `beta_hat` for every draw and every forecast step.  Useful for pinning a
#'     coefficient to zero or to a theory-driven value.}
#'   \item{Length-`horizon` numeric vector}{A different value is used at each
#'     forecast step (same across draws).  Useful for specifying a step-by-step
#'     trajectory for the coefficient.}
#'   \item{Function `function(h, beta_hat)`}{Called once per `(draw, step)`
#'     pair with the forecast step index `h` (`1..horizon`) and the per-draw
#'     fitted estimate `beta_hat`.  Use this form to express the override
#'     relative to each draw's own fitted value, e.g.
#'     `function(h, beta_hat) beta_hat * 0.5` halves the coefficient in every
#'     draw independently.  The function must return a single finite number.}
#'   \item{`NULL`}{Clears a previously applied override, restoring the
#'     per-draw fitted estimates.}
#' }
#'
#' @section Supported models:
#' Overrides are only supported for `linear` models (`block$linear == TRUE`).
#' Calling `coef_override()` on a `glm`, `gamlss`, or other non-linear
#' block errors immediately with a clear message.  Check `block$linear` first
#' if the model type is uncertain.
#'
#' @param block An `endogenr_coef_param` from [setup_param()] (`sp$<oc>$coef`).
#' @param term Character.  Name of the coefficient to override; must be in
#'   `colnames(block$beta_by_draw)`.  Use `colnames(sp$<oc>$coef$beta_by_draw)`
#'   or inspect `block$estimates` to see valid names.
#' @param value Override specification: a scalar, a length-`horizon` numeric
#'   vector, a `function(h, beta_hat)`, or `NULL` to clear.
#'
#' @return The same `endogenr_coef_param` block with `overrides[[term]]` and
#'   `effective[[term]]` updated (or removed when `value = NULL`).
#' @seealso [setup_param()], [simulate_system()]
#' @family simulation
#' @export
#'
#' @examples
#' \dontrun{
#' # Build a simple AR(1) + Mundlak means model
#' dt <- sim_panel_ar1(units = 8L, n_time = 20L, seed = 1L)
#' dt[, m_x := hist_mean(x), by = unit]
#' sys <- setup_system(
#'   c(list(build_model("exogen", formula = ~x)),
#'     mundlak_means("x"),
#'     list(build_model("linear", formula = y ~ lag(y) + m_x))),
#'   data = dt, train_start = 1, test_start = 16, horizon = 4,
#'   groupvar = "unit", timevar = "time", inner_sims = 2L
#' )
#' fit <- fit_system(sys, nsim = 4L)
#' sp  <- setup_param(fit)
#'
#' # 1. Scalar override — pin the m_x coefficient to zero for all trajectories
#' sp$y$coef <- coef_override(sp$y$coef, "m_x", 0)
#' sp$y$coef$effective$m_x   # 4 x 4 matrix of zeros
#'
#' # 2. Vector override — ramp the m_x coefficient up across the horizon
#' sp$y$coef <- coef_override(sp$y$coef, "m_x", c(0.2, 0.4, 0.6, 0.8))
#' sp$y$coef$effective$m_x[1, ]  # draw 1: 0.2, 0.4, 0.6, 0.8
#'
#' # 3. Function override — halve each draw's own fitted coefficient
#' sp$y$coef <- coef_override(sp$y$coef, "m_x",
#'                             function(h, beta_hat) beta_hat * 0.5)
#' # effective$m_x differs across rows (draws) but is 0.5 * beta_hat for each
#'
#' # 4. Function override — decay the AR coefficient by h/(H+1) each step
#' sp$y$coef <- coef_override(
#'   sp$y$coef, "lag_y",
#'   function(h, beta_hat) beta_hat * (1 - h / (sp$y$coef$dims$horizon + 1))
#' )
#'
#' # 5. Clear the override (restores per-draw fitted estimates)
#' sp$y$coef <- coef_override(sp$y$coef, "m_x", NULL)
#' is.null(sp$y$coef$effective$m_x)  # TRUE
#'
#' # Simulate with the active override
#' sp$y$coef <- coef_override(sp$y$coef, "m_x", 0)
#' res <- simulate_system(fit, scenario_params = sp)
#' }
coef_override <- function(block, term, value) {
  if (!isTRUE(block$overridable)) {
    stop("coefficient overrides are supported only for `linear` and `gamlss` models",
         call. = FALSE)
  }
  if (!term %in% colnames(block$beta_by_draw)) {
    stop("'", term, "' is not a coefficient in this model. Valid names: ",
         paste(colnames(block$beta_by_draw), collapse = ", "), call. = FALSE)
  }

  # NULL clears the override
  if (is.null(value)) {
    block$overrides[[term]] <- NULL
    block$effective[[term]] <- NULL
    return(block)
  }

  dims    <- block$dims
  nsim    <- dims$nsim
  horizon <- dims$horizon

  if (!is.function(value) &&
      !(is.numeric(value) && length(value) %in% c(1L, horizon))) {
    stop("coefficient override must be a scalar, a length-", horizon,
         " numeric vector, or a function(h, beta_hat)",
         call. = FALSE)
  }

  # Bake effective matrix: nsim x horizon
  eff_mat <- matrix(NA_real_, nrow = nsim, ncol = horizon)
  for (i in seq_len(nsim)) {
    beta_hat_i <- block$beta_by_draw[i, term]
    for (h in seq_len(horizon)) {
      val <- if (is.function(value)) {
        as.numeric(value(h, beta_hat_i))
      } else if (length(value) == 1L) {
        as.numeric(value)
      } else {
        as.numeric(value[h])
      }
      if (!is.finite(val)) {
        stop("coefficient override resolved to a non-finite value at step h = ", h,
             " for draw i = ", i, ".", call. = FALSE)
      }
      eff_mat[i, h] <- val
    }
  }

  block$overrides[[term]] <- value
  block$effective[[term]] <- eff_mat
  block
}

# --- Internal: slice scenario_params for one outer draw ----------------------

# Returns a named-by-outcome list of per-draw scenario slices. Each entry:
#   time_fe = list(var, ref, offsets = c(inner_sims, horizon)) or NULL
#   unit_fe = list(var, ref, per_trajectory, eff) or NULL
#   coef    = list(effective = named list term -> numeric(horizon))
.slice_scenario <- function(scenario_params, i) {
  lapply(scenario_params, function(entry) {
    # Time FE
    time_fe_s <- if (!is.null(entry$time_fe) && isTRUE(entry$time_fe$active)) {
      v <- entry$time_fe$values
      list(
        var     = entry$time_fe$var,
        ref     = entry$time_fe$ref,
        # drop = FALSE keeps c(inner_sims, horizon) when either is 1
        offsets = array(v[i, , , drop = FALSE], dim = dim(v)[2:3])
      )
    } else NULL

    # Unit FE
    unit_fe_s <- if (!is.null(entry$unit_fe) && isTRUE(entry$unit_fe$active)) {
      pt <- isTRUE(attr(entry$unit_fe$values, "per_trajectory"))
      eff <- if (pt) {
        v <- entry$unit_fe$values
        # c(n_units, inner_sims, horizon); rownames are looked up by predict
        array(v[, i, , , drop = FALSE], dim = dim(v)[c(1L, 3L, 4L)],
              dimnames = list(dimnames(v)[[1L]], NULL, NULL))
      } else {
        entry$unit_fe$values            # c(n_units, horizon), same for all draws
      }
      list(
        var            = entry$unit_fe$var,
        ref            = entry$unit_fe$ref,
        per_trajectory = pt,
        eff            = eff
      )
    } else NULL

    # Coef
    coef_s <- if (!is.null(entry$coef)) {
      list(effective = lapply(entry$coef$effective, function(m) m[i, ]))
    } else {
      list(effective = list())
    }

    list(time_fe = time_fe_s, unit_fe = unit_fe_s, coef = coef_s)
  })
}

# --- setup_param: exported --------------------------------------------------

#' Set up precomputed scenario parameters for a fitted endogenr system
#'
#' Builds an `endogenr_scenario_params` object keyed by model outcome (all
#' models in the fitted system, including independent models).  All stochastic
#' parameters — time-FE draw trajectories — are **baked at build time** using
#' the current RNG state, so the exact values entering the simulation are
#' inspectable and modifiable before calling [simulate_system()].
#'
#' @section Default strategies:
#' \describe{
#'   \item{Time FE (`factor(timevar)`)}{[fe_ar()] with `active = TRUE` — per
#'     `(outer draw, inner sim)` trajectory, the year effects are perturbed by
#'     their estimation covariance and forecast by an AR(1) from the last
#'     training-year effect.  Falls back to [fe_resample()] when fewer than 3
#'     year effects are estimated.  Baked as a `c(nsim, inner_sims, horizon)`
#'     array.}
#'   \item{Unit FE (`factor(unitvar)`)}{[fe_persist()] with `active = FALSE` —
#'     `factor(unit)` stays in the design matrix unchanged.  Equivalent to
#'     running without any scenario parameter; the baked `c(n_units, horizon)`
#'     matrix is present for inspection but not applied.}
#' }
#'
#' @section Object structure:
#' The returned object is a named list keyed by outcome, where each entry
#' contains:
#' \describe{
#'   \item{`type`}{Model type string (e.g. `"linear"`, `"glm_endogenr"`).}
#'   \item{`outcome`}{Outcome variable name.}
#'   \item{`independent`}{`TRUE` for models that do not feed back into the
#'     system (e.g. `exogen`).}
#'   \item{`adjustable`}{`TRUE` when at least one of `time_fe`, `unit_fe`, or
#'     `coef` is non-`NULL`.}
#'   \item{`time_fe`}{An `endogenr_fe_param` or `NULL`.  Key sub-fields:
#'     `kind` (`"time"`), `term` (e.g. `"factor(time)"`), `var`, `ref`,
#'     `levels`, `effects` (named numeric, baseline = 0), `effects_by_draw`
#'     (length-`nsim` list of per-draw effects), `vcov_by_draw` (length-`nsim`
#'     list of per-draw effect covariances, `NULL` entries when unavailable),
#'     `dims`, `active`, `heuristic`, `values` (baked array).}
#'   \item{`unit_fe`}{An `endogenr_fe_param` or `NULL`.  Same structure as
#'     `time_fe` but `kind = "unit"` and `values` is `c(n_units, horizon)`.}
#'   \item{`coef`}{An `endogenr_coef_param` or `NULL`.  Key sub-fields:
#'     `estimates` (named numeric, draw-1 β̂), `beta_by_draw` (`nsim × p`
#'     matrix of per-draw fitted coefficients), `linear` (logical),
#'     `overrides` (user spec, for reference), `effective` (named list of
#'     `nsim × horizon` baked override matrices — empty until
#'     [coef_override()] is applied).}
#' }
#' The object also carries a `dims` attribute:
#' `list(nsim, inner_sims, horizon, test_start, units, timevar, unitvar)`.
#'
#' @section Typical workflow:
#' \enumerate{
#'   \item Call `set.seed()` then `setup_param(fit)` to bake stochastic
#'     parameters reproducibly.
#'   \item Inspect the returned object (`print(sp)`, `sp$<oc>$time_fe$values`,
#'     `sp$<oc>$coef$estimates`).
#'   \item Optionally replace FE blocks with a different heuristic:
#'     `sp$<oc>$time_fe <- fe_distribution(sp$<oc>$time_fe, sd = 0.02)`.
#'   \item Optionally bake coefficient overrides:
#'     `sp$<oc>$coef <- coef_override(sp$<oc>$coef, "lag_y", 0.5)`.
#'   \item Pass to `simulate_system(fit, scenario_params = sp)`.
#' }
#' Omitting `scenario_params` (or passing `NULL`) causes `simulate_system()`
#' to call `setup_param()` internally under its own RNG state, which reproduces
#' varied time-FE offsets exactly as before the redesign.
#'
#' @param fitted_system An `endogenr_fitted_system` from [fit_system()].
#'
#' @return An `endogenr_scenario_params` object (S3 class).
#' @seealso [fe_ar()], [fe_resample()], [fe_distribution()], [fe_fixed()],
#'   [fe_persist()], [fe_converge()], [coef_override()], [simulate_system()]
#' @family simulation
#' @export
#'
#' @examples
#' \dontrun{
#' # -- Minimal time-FE example -----------------------------------------------
#' dt  <- sim_panel_common_shock(units = 8L, n_time = 20L, seed = 1L)
#' sys <- setup_system(
#'   list(build_model("linear", formula = y ~ lag(y) + factor(time))),
#'   data = dt, train_start = 1, test_start = 16, horizon = 4,
#'   groupvar = "unit", timevar = "time", inner_sims = 3L
#' )
#' fit <- fit_system(sys, nsim = 5L)
#'
#' # Bake defaults (fe_ar for time FE) under a fixed seed
#' set.seed(1)
#' sp <- setup_param(fit)
#' print(sp)                          # human-readable summary
#' attr(sp, "dims")                   # nsim / inner_sims / horizon / test_start
#'
#' # Inspect baked values
#' dim(sp$y$time_fe$values)           # c(5, 3, 4) = c(nsim, inner_sims, horizon)
#' sp$y$time_fe$heuristic$strategy    # "ar1"
#' sp$y$coef$estimates                # named vector of representative β̂
#'
#' # Modify time-FE strategy: narrow Normal draws
#' set.seed(1)
#' sp <- setup_param(fit)
#' sp$y$time_fe <- fe_distribution(sp$y$time_fe, sd = 0.02)
#' res_narrow <- simulate_system(fit, scenario_params = sp)
#'
#' # Modify time-FE strategy: converge to zero over the horizon
#' set.seed(1)
#' sp <- setup_param(fit)
#' sp$y$time_fe <- fe_converge(sp$y$time_fe, to = 0, path = "linear")
#'
#' # Apply a coefficient override then simulate
#' set.seed(1)
#' sp <- setup_param(fit)
#' sp$y$coef <- coef_override(sp$y$coef, "lag_y",
#'                             function(h, beta_hat) beta_hat * 0.5)
#' res_ov <- simulate_system(fit, scenario_params = sp)
#'
#' # Determinism: same seed => identical baked arrays
#' set.seed(99); a <- setup_param(fit)
#' set.seed(99); b <- setup_param(fit)
#' identical(a$y$time_fe$values, b$y$time_fe$values)  # TRUE
#' }
setup_param <- function(fitted_system) {
  if (!inherits(fitted_system, "endogenr_fitted_system")) {
    stop("`fitted_system` must be the output of fit_system().", call. = FALSE)
  }

  nsim       <- fitted_system$nsim
  if (is.null(nsim)) nsim <- length(fitted_system$fitted_draws)
  inner_sims <- fitted_system$inner_sims
  horizon    <- fitted_system$horizon
  test_start <- fitted_system$test_start
  groupvar   <- fitted_system$groupvar
  timevar    <- fitted_system$timevar
  units      <- as.character(unique(fitted_system$simulation_data[[groupvar]]))

  dims <- list(
    nsim       = as.integer(nsim),
    inner_sims = as.integer(inner_sims),
    horizon    = as.integer(horizon),
    test_start = as.integer(test_start),
    units      = units,
    timevar    = timevar,
    unitvar    = groupvar
  )

  fitted_draws  <- fitted_system$fitted_draws
  fitted_models <- fitted_system$fitted_models

  # Locate the model for a given outcome within one draw's model list.
  find_draw_model <- function(draw, outcome) {
    for (m in draw) {
      if (!is.null(m$outcome) && identical(m$outcome, outcome)) return(m)
    }
    NULL
  }

  entries <- list()

  for (model in fitted_models) {
    oc    <- model$outcome
    type  <- class(model)[1L]
    terms <- scenario_terms(model)
    time_fe_rep <- terms$time_fe
    unit_fe_rep <- terms$unit_fe

    has_coefs   <- !is.null(model$coefs)
    has_time_fe <- !is.null(time_fe_rep)
    has_unit_fe <- !is.null(unit_fe_rep)
    adjustable  <- has_coefs || has_time_fe || has_unit_fe

    # -- Time FE block ---------------------------------------------------------
    time_fe_blk <- if (has_time_fe) {
      ebd <- lapply(fitted_draws, function(draw) {
        m <- find_draw_model(draw, oc)
        if (!is.null(m) && !is.null(m$time_fe)) m$time_fe$effects
        else time_fe_rep$effects
      })
      vbd <- lapply(fitted_draws, function(draw) {
        m <- find_draw_model(draw, oc)
        if (!is.null(m) && !is.null(m$time_fe)) m$time_fe$vcov
        else time_fe_rep$vcov
      })
      blk <- structure(
        list(
          kind            = "time",
          term            = time_fe_rep$term,
          var             = time_fe_rep$timevar,
          ref             = time_fe_rep$ref_value,
          levels          = names(time_fe_rep$effects),
          effects         = time_fe_rep$effects,
          effects_by_draw = ebd,
          vcov_by_draw    = vbd,
          dims            = dims,
          active          = TRUE,
          heuristic       = NULL,
          values          = NULL
        ),
        class = "endogenr_fe_param"
      )
      # default: AR(1) from the origin; resample when the series is too short
      if (all(lengths(ebd) >= 3L)) fe_ar(blk) else fe_resample(blk)
    } else NULL

    # -- Unit FE block ---------------------------------------------------------
    unit_fe_blk <- if (has_unit_fe) {
      ebd_u <- lapply(fitted_draws, function(draw) {
        m <- find_draw_model(draw, oc)
        if (!is.null(m) && !is.null(m$unit_fe)) m$unit_fe$effects
        else unit_fe_rep$effects
      })
      blk <- structure(
        list(
          kind            = "unit",
          term            = unit_fe_rep$term,
          var             = unit_fe_rep$unitvar,
          ref             = unit_fe_rep$ref_value,
          levels          = names(unit_fe_rep$effects),
          effects         = unit_fe_rep$effects,
          effects_by_draw = ebd_u,
          dims            = dims,
          active          = FALSE,
          heuristic       = NULL,
          values          = NULL
        ),
        class = "endogenr_fe_param"
      )
      fe_persist(blk)  # default: persist (active = FALSE, identity)
    } else NULL

    # -- Coef block ------------------------------------------------------------
    coef_blk <- if (has_coefs) {
      all_terms <- model$coefs$term

      # Build nsim x p beta matrix from per-draw fitted coefficients.
      beta_mat <- do.call(rbind, lapply(fitted_draws, function(draw) {
        m <- find_draw_model(draw, oc)
        if (!is.null(m) && !is.null(m$fitted)) {
          cf <- stats::coef(m$fitted)
          cf[all_terms]
        } else {
          stats::setNames(rep(NA_real_, length(all_terms)), all_terms)
        }
      }))
      colnames(beta_mat) <- all_terms

      # Representative estimates from draw 1
      m1 <- find_draw_model(fitted_draws[[1L]], oc)
      estimates <- if (!is.null(m1) && !is.null(m1$fitted)) {
        cf1 <- stats::coef(m1$fitted)
        cf1[all_terms]
      } else {
        stats::setNames(rep(NA_real_, length(all_terms)), all_terms)
      }

      structure(
        list(
          estimates    = estimates,
          beta_by_draw = beta_mat,
          linear       = identical(type, "linear"),
          overridable  = type %in% c("linear", "endogenr_gamlss"),
          dims         = dims,
          overrides    = list(),
          effective    = list()
        ),
        class = "endogenr_coef_param"
      )
    } else NULL

    # Multi-outcome models (e.g. multi-column exogen) produce a vector `oc`;
    # `[[<-` with a vector tries nested indexing, so iterate explicitly.
    for (o in oc) {
      entries[[o]] <- list(
        type        = type,
        outcome     = o,
        independent = isTRUE(model$independent),
        adjustable  = adjustable,
        time_fe     = time_fe_blk,
        unit_fe     = unit_fe_blk,
        coef        = coef_blk
      )
    }
  }

  structure(entries, class = "endogenr_scenario_params", dims = dims)
}

# --- print.endogenr_scenario_params -----------------------------------------

#' @exportS3Method print endogenr_scenario_params
print.endogenr_scenario_params <- function(x, ...) {
  dims <- attr(x, "dims")
  cat("<endogenr_scenario_params>\n")
  if (!is.null(dims)) {
    cat(sprintf(
      "dims: nsim=%d  inner_sims=%d  horizon=%d  test_start=%d  n_units=%d\n",
      dims$nsim, dims$inner_sims, dims$horizon, dims$test_start, length(dims$units)
    ))
  }
  if (length(x) == 0L) {
    cat("  (no parameterised outcomes)\n")
    return(invisible(x))
  }

  for (oc in names(x)) {
    e       <- x[[oc]]
    adj_flg <- if (isTRUE(e$adjustable)) "*" else " "
    cat(sprintf("\n-- outcome: %s  [%s]%s --\n", oc, e$type, adj_flg))

    if (!is.null(e$time_fe)) {
      tfe  <- e$time_fe
      heur <- tfe$heuristic
      smry <- c(mean = mean(tfe$effects), sd = stats::sd(tfe$effects),
                min  = min(tfe$effects),  max = max(tfe$effects))
      cat(sprintf("  time_fe: term='%s', n_levels=%d, strategy='%s'\n",
                  tfe$term, length(tfe$levels),
                  if (!is.null(heur)) heur$strategy else "?"))
      cat(sprintf("    effects: mean=%.3f, sd=%.3f, min=%.3f, max=%.3f\n",
                  smry["mean"], smry["sd"], smry["min"], smry["max"]))
      if (!is.null(tfe$values)) {
        v <- tfe$values
        cat(sprintf("    values[nsim x inner_sims x horizon]: mean=%.3f, sd=%.3f\n",
                    mean(v), stats::sd(v)))
      }
    }

    if (!is.null(e$unit_fe)) {
      ufe  <- e$unit_fe
      heur <- ufe$heuristic
      cat(sprintf("  unit_fe: term='%s', n_units=%d, strategy='%s', active=%s\n",
                  ufe$term, length(ufe$levels),
                  if (!is.null(heur)) heur$strategy else "?",
                  if (isTRUE(ufe$active)) "TRUE" else "FALSE"))
    }

    if (!is.null(e$coef)) {
      cb  <- e$coef
      ovr <- names(cb$effective)
      cat(sprintf("  coef: %d term(s)", ncol(cb$beta_by_draw)))
      if (length(ovr) > 0L) {
        cat(sprintf(", overrides: %s", paste(ovr, collapse = ", ")))
      }
      cat("\n")
    }
  }

  # Workflow hints for the first adjustable outcome
  oc1 <- names(x)[1L]
  if (!is.null(x[[oc1]])) {
    cat("\nHints:\n")
    if (!is.null(x[[oc1]]$time_fe)) {
      cat("  sp$", oc1, "$time_fe <- fe_distribution(sp$", oc1,
          "$time_fe, sd = ..)\n", sep = "")
      cat("  sp$", oc1, "$time_fe <- fe_converge(sp$", oc1,
          "$time_fe, to = 0, path = \"linear\")\n", sep = "")
    }
    if (!is.null(x[[oc1]]$unit_fe)) {
      cat("  sp$", oc1, "$unit_fe <- fe_converge(sp$", oc1,
          "$unit_fe, to = 0)\n", sep = "")
    }
    if (!is.null(x[[oc1]]$coef)) {
      term1 <- colnames(x[[oc1]]$coef$beta_by_draw)[1L]
      cat("  sp$", oc1, "$coef <- coef_override(sp$", oc1,
          "$coef, \"", term1, "\", 0)  ",
          "# scalar, length-horizon vector, or function(h, beta_hat)\n", sep = "")
    }
  }

  invisible(x)
}
