# Scenario parameter framework -----------------------------------------------
#
# Precomputed, model-keyed scenario parameters with fixed-effect heuristic
# helpers. All stochastic parameters are baked at build time via setup_param();
# simulate_system() reads the baked values instead of drawing lazily.

# --- Internal: detect factor(timevar) in a fitted linear model ---------------

# Scan fit_formula for a term of the form factor(<timevar>). Extract the
# estimated year effects from fitted_lm and return metadata, or NULL.
.detect_time_fe <- function(fit_formula, fitted_lm, timevar) {
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
  xl <- fitted_lm$xlevels[[fe_label]]
  if (is.null(xl)) return(NULL)

  # Build effects vector: baseline (xl[1]) -> 0, remaining levels -> coefficient.
  # Effects are named by level so downstream helpers can look them up by name.
  effects <- stats::setNames(numeric(length(xl)), xl)
  for (lv in xl[-1L]) {
    nm  <- paste0(fe_label, lv)
    val <- stats::coef(fitted_lm)[[nm]]
    if (!is.null(val)) effects[[lv]] <- val
  }
  ref_value <- as.numeric(xl[1L])

  list(
    term      = fe_label,
    timevar   = timevar,
    ref_value = ref_value,
    effects   = effects  # named by time level
  )
}

# --- Internal: detect factor(unitkey) in a fitted linear model ---------------
#
# Scans fit_formula for a term of the form factor(<key>) where <key> is any
# element of unit_keys. Returns unit-FE metadata or NULL. Only the first
# matching key is used (single-key unit FE). Multi-key composite unit FE is
# out of scope: if unit_keys has length > 1, the factor over any one key is
# detected and the rest are left in the design unchanged.
.detect_unit_fe <- function(fit_formula, fitted_lm, unit_keys) {
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

  xl <- fitted_lm$xlevels[[fe_label]]
  if (is.null(xl)) return(NULL)

  # Build effects named by unit level; baseline (xl[1]) = 0.
  effects <- stats::setNames(numeric(length(xl)), xl)
  for (lv in xl[-1L]) {
    nm  <- paste0(fe_label, lv)
    val <- stats::coef(fitted_lm)[[nm]]
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
#' metadata list or `NULL`. This is the per-model-type extension seam consumed
#' by [setup_param()]. Currently only `linear` models carry FE metadata; all
#' other types return `NULL` for both elements via the default method. Add
#' support for additional types by implementing `scenario_terms.<class>`.
#'
#' @param model A fitted endogenr model object.
#' @return A list with elements `time_fe` and `unit_fe` (each a metadata list
#'   or `NULL`).
#' @family simulation
#' @export
scenario_terms <- function(model) UseMethod("scenario_terms")

#' @rdname scenario_terms
#' @exportS3Method scenario_terms default
scenario_terms.default <- function(model) list(time_fe = NULL, unit_fe = NULL)

#' @rdname scenario_terms
#' @exportS3Method scenario_terms linear
scenario_terms.linear <- function(model) {
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

# --- Exported heuristic helpers: FE blocks -----------------------------------

#' Fixed-effect heuristic: resample from estimated effects
#'
#' Re-bakes the FE block's `values` by sampling with replacement from the
#' per-draw estimated effects, independently for every `(inner_sim, step)`.
#' This is the **default strategy for time fixed effects**.
#'
#' For unit FE, this produces a full
#' `c(n_units, nsim, inner_sims, horizon)` array (`per_trajectory = TRUE`);
#' be aware of the memory cost for large panels.
#'
#' @param block An `endogenr_fe_param` from [setup_param()].
#' @param to Convergence target after the resample baseline: `NULL` (no
#'   convergence), a numeric scalar, `"zero"`, `"mean"`, or a character vector
#'   of unit ids whose mean effect is the target.
#' @param factor Divergence/shrink multiplier applied along `path` when `to`
#'   is `NULL` (default `1` = no change; `> 1` diverge; `< 1` shrink).
#' @param path Movement path: `"constant"` (default, no movement),
#'   `"linear"`, or `"sigmoid"`.
#' @param mid Sigmoid inflection step (default: midpoint of the horizon).
#' @param steep Sigmoid steepness (default `1`).
#'
#' @return The modified `endogenr_fe_param` block with updated `values`,
#'   `active`, and `heuristic`.
#' @family simulation
#' @export
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
      pool <- block$effects_by_draw[[i]]
      B    <- matrix(sample(pool, inner_sims * horizon, replace = TRUE),
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
    for (i in seq_len(nsim)) {
      pool <- block$effects_by_draw[[i]]
      for (u_idx in seq_len(n_units)) {
        B <- matrix(sample(pool, inner_sims * horizon, replace = TRUE),
                    nrow = inner_sims, ncol = horizon)
        for (h in seq_len(horizon)) {
          vals[u_idx, i,, h] <- if (!is.null(t_val))
            (1 - w[h]) * B[, h] + w[h] * t_val
          else
            B[, h] * (1 + (factor - 1) * w[h])
        }
      }
    }
    attr(vals, "per_trajectory") <- TRUE
  }

  block$values    <- vals
  block$active    <- TRUE
  block$heuristic <- list(strategy = "resample", to = to, factor = factor,
                          path = path, mid = mid, steep = steep)
  block
}

#' Fixed-effect heuristic: normal distribution
#'
#' Re-bakes the FE block's `values` by drawing from a Normal distribution with
#' per-draw mean and sd estimated from `effects_by_draw`, independently for
#' every `(inner_sim, step)`.
#'
#' @param block An `endogenr_fe_param` from [setup_param()].
#' @param mean Optional scalar mean override (default: per-draw mean of
#'   estimated effects).
#' @param sd Optional scalar sd override (default: per-draw sd of estimated
#'   effects).
#' @inheritParams fe_resample
#'
#' @return The modified `endogenr_fe_param` block.
#' @family simulation
#' @export
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
    for (i in seq_len(nsim)) {
      eff_i <- block$effects_by_draw[[i]]
      m_i   <- if (!is.null(mean)) mean else base::mean(eff_i)
      s_i   <- if (!is.null(sd))   sd   else stats::sd(eff_i)
      for (u_idx in seq_len(n_units)) {
        B <- matrix(stats::rnorm(inner_sims * horizon, m_i, s_i),
                    nrow = inner_sims, ncol = horizon)
        for (h in seq_len(horizon)) {
          vals[u_idx, i,, h] <- if (!is.null(t_val))
            (1 - w[h]) * B[, h] + w[h] * t_val
          else
            B[, h] * (1 + (factor - 1) * w[h])
        }
      }
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
#' Re-bakes the FE block's `values` with a deterministic constant baseline.
#' The array is shaped `c(nsim, inner_sims, horizon)` for time FE and
#' `c(n_units, horizon)` for unit FE, preserving the uniform consumer contract.
#'
#' @param block An `endogenr_fe_param` from [setup_param()].
#' @param value Numeric scalar baseline (default `0`).
#' @inheritParams fe_resample
#'
#' @return The modified `endogenr_fe_param` block.
#' @family simulation
#' @export
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
#' Sets each simulation unit's fixed-effect offset to its historically
#' estimated value, optionally converging toward a target. With all defaults
#' (`to = NULL`, `factor = 1`, `path = "constant"`), `active = FALSE` and
#' `factor(unit)` stays in the design matrix **unchanged** — identical to
#' today's behaviour with no scenario modification.
#'
#' @param block An `endogenr_fe_param` with `kind == "unit"`.
#' @inheritParams fe_resample
#'
#' @return The modified `endogenr_fe_param` block (`active = FALSE` under all
#'   defaults; `active = TRUE` when any non-default argument is supplied).
#' @family simulation
#' @export
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
#' Sugar over the path mechanics: applies a path-shaped convergence toward
#' `to` from the block's natural baseline. For unit FE the baseline per unit
#' is its historically estimated effect. For time FE a single resample draw per
#' trajectory is held constant across steps (so trajectories vary but each one
#' is internally consistent).
#'
#' @param block An `endogenr_fe_param` from [setup_param()].
#' @param to Convergence target (required): numeric scalar, `"zero"`,
#'   `"mean"`, or a character vector of unit ids whose mean effect is the
#'   target.
#' @param path `"linear"` (default), `"constant"`, or `"sigmoid"`.
#' @param mid Sigmoid inflection step (default: midpoint of the horizon).
#' @param steep Sigmoid steepness (default `1`).
#'
#' @return The modified `endogenr_fe_param` block with `active = TRUE`.
#' @family simulation
#' @export
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
      pool      <- block$effects_by_draw[[i]]
      base_draw <- sample(pool, inner_sims, replace = TRUE)  # one per inner_sim
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
#' Bakes overridden coefficient trajectories into the coef block, replacing
#' the per-draw fitted estimate for the named term at every forecast step.
#' `value` may be:
#' \itemize{
#'   \item A **scalar** — constant across all draws and forecast steps.
#'   \item A **length-`horizon` numeric vector** — per-step, same across draws.
#'   \item A **`function(h, beta_hat)`** — returns the override β* at step `h`
#'     given the per-draw fitted estimate `beta_hat`.
#' }
#' Pass `value = NULL` to clear a previously applied override.
#'
#' Overrides are only supported for `linear` models; calling this on a
#' non-linear block errors immediately.
#'
#' @param block An `endogenr_coef_param` from [setup_param()].
#' @param term Character. Coefficient name (must be in
#'   `colnames(block$beta_by_draw)`).
#' @param value Override specification or `NULL` to clear.
#'
#' @return The modified `endogenr_coef_param` block.
#' @family simulation
#' @export
coef_override <- function(block, term, value) {
  if (!block$linear) {
    stop("coefficient overrides are only supported for `linear` models",
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
      list(
        var     = entry$time_fe$var,
        ref     = entry$time_fe$ref,
        offsets = entry$time_fe$values[i, , ]  # c(inner_sims, horizon)
      )
    } else NULL

    # Unit FE
    unit_fe_s <- if (!is.null(entry$unit_fe) && isTRUE(entry$unit_fe$active)) {
      pt <- isTRUE(attr(entry$unit_fe$values, "per_trajectory"))
      eff <- if (pt) {
        entry$unit_fe$values[, i, , ]   # c(n_units, inner_sims, horizon)
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
#' models in the fitted system, including independent models). All stochastic
#' parameters (time-FE draw trajectories) are **baked at build time** using
#' the current RNG state, so results are inspectable before
#' [simulate_system()].
#'
#' Default strategies applied at build time:
#' \itemize{
#'   \item **Time FE** (`factor(timevar)`): [fe_resample()] — one independent
#'     draw per `(outer draw, inner sim, forecast step)`.
#'   \item **Unit FE** (`factor(unitvar)`): [fe_persist()] with `active =
#'     FALSE` — `factor(unit)` stays in the design matrix unchanged;
#'     equivalent to the current default behaviour.
#' }
#'
#' Each entry contains:
#' \describe{
#'   \item{`type`}{Model type string.}
#'   \item{`outcome`}{Outcome variable name.}
#'   \item{`independent`}{Logical.}
#'   \item{`adjustable`}{Logical: any non-NULL adjustable block present.}
#'   \item{`time_fe`}{`endogenr_fe_param` or `NULL`.}
#'   \item{`unit_fe`}{`endogenr_fe_param` or `NULL`.}
#'   \item{`coef`}{`endogenr_coef_param` or `NULL`.}
#' }
#'
#' The returned object carries a `dims` attribute:
#' `list(nsim, inner_sims, horizon, test_start, units, timevar, unitvar)`.
#'
#' @param fitted_system An `endogenr_fitted_system` from [fit_system()].
#'
#' @return An `endogenr_scenario_params` object.
#' @seealso [fe_resample()], [fe_distribution()], [fe_fixed()], [fe_persist()],
#'   [fe_converge()], [coef_override()], [simulate_system()]
#' @family simulation
#' @export
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
      blk <- structure(
        list(
          kind            = "time",
          term            = time_fe_rep$term,
          var             = time_fe_rep$timevar,
          ref             = time_fe_rep$ref_value,
          levels          = names(time_fe_rep$effects),
          effects         = time_fe_rep$effects,
          effects_by_draw = ebd,
          dims            = dims,
          active          = TRUE,
          heuristic       = NULL,
          values          = NULL
        ),
        class = "endogenr_fe_param"
      )
      fe_resample(blk)  # default: resample independently per (sim, step)
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
          dims         = dims,
          overrides    = list(),
          effective    = list()
        ),
        class = "endogenr_coef_param"
      )
    } else NULL

    entries[[oc]] <- list(
      type        = type,
      outcome     = oc,
      independent = isTRUE(model$independent),
      adjustable  = adjustable,
      time_fe     = time_fe_blk,
      unit_fe     = unit_fe_blk,
      coef        = coef_blk
    )
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
