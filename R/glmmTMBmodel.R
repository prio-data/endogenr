# glmmTMB model support for endogenr
#
# Integrates glmmTMB mixed-effects models as step-by-step dependent models
# inside the endogenr dynamic simulation loop. The model is predicted one time
# step at a time via predict.glmmTMB(newdata=…), exactly like glm/heterolm.
#
# Two-stage approach mirroring heterolmmodel.R:
#   Stage 1: shared AST-level materialization of lag/roll/cum sub-expressions
#            across the main formula, dispformula, and ziformula into .pt#
#            columns via panel_materialize.
#   Stage 2: pooled glmmTMB::glmmTMB fit on the rewritten aliased formula.
#
# IMPORTANT: fitted class is "glmmTMB_endogenr" (not "glmmTMB") to avoid
# clobbering glmmTMB's own predict.glmmTMB S3 method.

# --------------------------------------------------------------------------
# Part 1: Formula decomposition helpers (pure base-R; no glmmTMB needed)
# --------------------------------------------------------------------------

# Test whether an expression is a random-effects bar (| or ||).
.glmmTMB_is_bar <- function(e) {
  is.call(e) && (identical(e[[1L]], as.name("|")) ||
                   identical(e[[1L]], as.name("||")))
}

# Classify one additive formula term.
# Returns list(effect, group, covstruct) for RE terms, NULL for fixed terms.
#
#  - Plain RE:    (1 + lag(dem) | region)  -> effect=`1 + lag(dem)`, group=`region`, covstruct=FALSE
#  - Cov-struct:  ar1(times + 0 | g)       -> effect=`times + 0`, group=`g`,      covstruct=TRUE
#  - Fixed:       lag(dem)                 -> NULL
.glmmTMB_as_re <- function(term) {
  # Peel a leading `(` wrapping call
  inner <- if (is.call(term) && identical(term[[1L]], as.name("("))) term[[2L]] else term

  if (.glmmTMB_is_bar(inner)) {
    return(list(effect = inner[[2L]], group = inner[[3L]], covstruct = FALSE))
  }

  # Cov-struct wrapper: a one-argument call whose argument (after peeling `(`)
  # is a bar — e.g. ar1(times + 0 | g), us(…|g), exp(pos+0|g).
  if (is.call(inner) && length(inner) == 2L) {
    arg <- inner[[2L]]
    arg <- if (is.call(arg) && identical(arg[[1L]], as.name("("))) arg[[2L]] else arg
    if (.glmmTMB_is_bar(arg)) {
      return(list(effect = arg[[2L]], group = arg[[3L]], covstruct = TRUE))
    }
  }

  NULL  # fixed term
}

# Walk an additive (+/-) tree, returning a list of individual terms.
# Subtracted terms are dropped (they carry no predictor variable).
# Shared by .glmmTMB_decompose and .gamlss_decompose.
.flatten_additive_terms <- function(expr) {
  if (!is.call(expr)) return(list(expr))
  head <- as.character(expr[[1L]])
  if (head == "+") {
    return(c(.flatten_additive_terms(expr[[2L]]),
             .flatten_additive_terms(expr[[3L]])))
  }
  if (head == "-" && length(expr) == 3L) {
    return(.flatten_additive_terms(expr[[2L]]))  # drop subtracted part
  }
  list(expr)
}

# Helper: is this a trivial (intercept-only or empty) RHS?
# Used to skip ~1 / ~0 dispformula/ziformula in materialization.
.is_trivial_rhs <- function(f) {
  is.null(f) ||
    identical(rlang::f_rhs(f), 1L) ||
    identical(rlang::f_rhs(f), 0L) ||
    identical(rlang::f_rhs(f), 1) ||
    identical(rlang::f_rhs(f), 0)
}

#' Decompose a glmmTMB formula into a bar-free graph formula
#'
#' Extracts every predictor variable (for dependency-graph edges) from the main
#' formula plus dispformula/ziformula. Random-effects grouping factors and
#' covariance-structure coordinates are treated as ordinary predictors: they are
#' read at every forecast step, so — like any predictor — they must be produced
#' by some model (e.g. an `exogen` that carries the grouping column forward).
#' Returns a bar-free graph_formula that `stats::terms()` and
#' `.edges_from_formula()` can consume safely.
#'
#' @param formula Main two-sided glmmTMB formula (may contain RE bars,
#'   cov-struct wrappers).
#' @param dispformula Optional one-sided dispersion formula.
#' @param ziformula Optional one-sided zero-inflation formula.
#'
#' @return A list with:
#'   \item{outcome}{Character: the LHS variable name.}
#'   \item{graph_formula}{A bar-free formula for dependency analysis.}
#' @keywords internal
.glmmTMB_decompose <- function(formula, dispformula = NULL, ziformula = NULL) {
  outcome <- all.vars(rlang::f_lhs(formula))

  pred_terms <- list()   # language objects to join with +

  parts <- Filter(function(f) inherits(f, "formula"),
                  list(formula, dispformula, ziformula))

  add_pred <- function(lst, expr) {
    if (length(all.vars(expr)) > 0L) c(lst, list(expr)) else lst
  }

  for (part in parts) {
    terms_list <- .flatten_additive_terms(rlang::f_rhs(part))
    for (term in terms_list) {
      re <- .glmmTMB_as_re(term)
      if (is.null(re)) {
        # Fixed term — a predictor when it carries a variable (not the intercept)
        pred_terms <- add_pred(pred_terms, term)
      } else {
        # RE bar (effect | group) or cov-struct (ar1(coords | group)): the effect
        # / coordinates AND the grouping factor are all predictors that must be
        # produced (carried forward), so they become dependency-graph edges.
        for (et in .flatten_additive_terms(re$effect)) pred_terms <- add_pred(pred_terms, et)
        pred_terms <- add_pred(pred_terms, re$group)
      }
    }
  }

  # Deduplicate pred_terms by deparse, preserving first occurrence
  seen <- character(0)
  unique_terms <- list()
  for (t in pred_terms) {
    key <- deparse(t, width.cutoff = 500L)
    if (!key %in% seen) { seen <- c(seen, key); unique_terms <- c(unique_terms, list(t)) }
  }
  pred_terms <- unique_terms

  # Build the bar-free graph RHS
  graph_rhs <- if (length(pred_terms) == 0L) {
    as.symbol("1")
  } else {
    Reduce(function(a, b) call("+", a, b), pred_terms)
  }

  graph_formula <- rlang::new_formula(
    rlang::f_lhs(formula),
    graph_rhs,
    env = rlang::f_env(formula)
  )

  list(outcome = outcome, graph_formula = graph_formula)
}


# --------------------------------------------------------------------------
# Part 2: Family response-draw helpers (no glmmTMB needed at definition time)
# --------------------------------------------------------------------------

# Per-family predictive draw from response scale.
# `mu`         — conditional mean vector (after linkinv).
# `family_name`— string from family(fitted)$family.
# `disp`       — per-row dispersion from predict(type="disp"),
#                as defined by sigma.glmmTMB for each family.
.glmmTMB_response_draw <- function(mu, family_name, disp, fam_pars = NULL) {
  n   <- length(mu)
  eps <- .Machine$double.eps

  switch(
    family_name,

    "gaussian" = {
      stats::rnorm(n, mean = mu, sd = pmax(disp, 0))
    },

    "poisson" = {
      stats::rpois(n, lambda = pmax(mu, 0))
    },

    "binomial" = {
      # No trial count in the grid — proportion / Bernoulli draw, matching GLM assumption
      stats::rbinom(n, size = 1L, prob = pmin(pmax(mu, 0), 1))
    },

    "Gamma" = {
      # sigma.glmmTMB: disp = 1/sqrt(shape)  =>  shape = 1/disp^2
      shape <- 1 / disp^2
      scale <- pmax(mu, eps) * disp^2
      stats::rgamma(n, shape = shape, scale = scale)
    },

    "nbinom2" = {
      # var = mu(1 + mu/theta), disp = theta (size in rnbinom)
      stats::rnbinom(n, mu = pmax(mu, 0), size = pmax(disp, eps))
    },

    "nbinom1" = {
      # var = mu(1 + alpha), disp = alpha => rnbinom size = mu/alpha
      mu_safe   <- pmax(mu, 0)
      disp_safe <- pmax(disp, eps)
      stats::rnbinom(n, mu = mu_safe, size = mu_safe / disp_safe)
    },

    "beta" = {
      # var = mu(1-mu)/(1+phi), disp = phi => Beta(a = mu*phi, b = (1-mu)*phi)
      mu2 <- pmin(pmax(mu, eps), 1 - eps)
      stats::rbeta(n,
                   shape1 = pmax(mu2 * disp, eps),
                   shape2 = pmax((1 - mu2) * disp, eps))
    },

    "betabinomial" = {
      # Same Beta parameterisation as "beta"
      mu2 <- pmin(pmax(mu, eps), 1 - eps)
      stats::rbeta(n,
                   shape1 = pmax(mu2 * disp, eps),
                   shape2 = pmax((1 - mu2) * disp, eps))
    },

    "t" = {
      # t_family: location = mu, scale = disp (sigma), df = family_params.
      # link is identity, so mu (=linkinv(eta)) is the location/mean (df>1).
      df <- as.numeric(fam_pars)[1L]
      if (!is.finite(df) || df <= 0) { .glmmTMB_warn_unsupported("t"); return(mu) }
      mu + pmax(disp, 0) * stats::rt(n, df = df)
    },

    "lognormal" = {
      # mu (=exp(eta)) is the response MEAN; disp is the response SD. Recover the
      # log-scale params so E[Y]=mu and SD[Y]=disp, then draw.
      mu_s  <- pmax(mu, eps)
      sdlog <- sqrt(log1p((disp / mu_s)^2))
      stats::rlnorm(n, meanlog = log(mu_s) - sdlog^2 / 2, sdlog = sdlog)
    },

    "skewnormal" = {
      # dskewnorm reparam (glmmTMB src/distrib.h): mu = mean, disp = sigma,
      # alpha = shape (family_params). Draw via the standard CSN construction.
      alpha <- as.numeric(fam_pars)[1L]
      if (!is.finite(alpha)) alpha <- 0
      delta <- alpha / sqrt(1 + alpha^2)
      omega <- pmax(disp, eps) / sqrt(1 - (2 / pi) * delta^2)
      xi    <- mu - omega * delta * sqrt(2 / pi)
      xi + omega * (delta * abs(stats::rnorm(n)) + sqrt(1 - delta^2) * stats::rnorm(n))
    },

    "truncated_poisson" = {
      # mu (=exp(eta)) is the UNtruncated rate; zero-truncate by inverse CDF.
      lambda <- pmax(mu, eps)
      stats::qpois(stats::runif(n, stats::dpois(0L, lambda), 1), lambda)
    },

    "truncated_nbinom2" = {
      # mu = untruncated mean, disp = theta (size); zero-truncate by inverse CDF.
      mu_s <- pmax(mu, eps)
      size <- pmax(disp, eps)
      stats::qnbinom(stats::runif(n, stats::dnbinom(0L, size = size, mu = mu_s), 1),
                     size = size, mu = mu_s)
    },

    "truncated_nbinom1" = {
      # nbinom1: disp = alpha, size = mu/alpha; zero-truncate by inverse CDF.
      mu_s <- pmax(mu, eps)
      size <- mu_s / pmax(disp, eps)
      stats::qnbinom(stats::runif(n, stats::dnbinom(0L, size = size, mu = mu_s), 1),
                     size = size, mu = mu_s)
    },

    "tweedie" = {
      # var = phi*mu^p: phi = disp, p = family_params ("Tweedie power").
      # Needs the tweedie package; fall back to the mean (warn once) if absent.
      p <- as.numeric(fam_pars)[1L]
      if (is.finite(p) && requireNamespace("tweedie", quietly = TRUE)) {
        tweedie::rtweedie(n, mu = pmax(mu, eps), phi = pmax(disp, eps), power = p)
      } else {
        .glmmTMB_warn_unsupported("tweedie")
        mu
      }
    },

    {
      # Unsupported family: warn once, fall back to point estimate (mean-only draw)
      .glmmTMB_warn_unsupported(family_name)
      mu
    }
  )
}

# One-time per-session warning for unsupported glmmTMB families.
.glmmTMB_warn_unsupported <- function(family_name) {
  key <- paste0(".endogenr_glmmTMB_warned_", family_name)
  if (!isTRUE(getOption(key))) {
    warning(
      "glmmTMB family '", family_name, "' has no response-scale predictive draw ",
      "in endogenr; falling back to the conditional mean (parameter-uncertainty-only draw), ",
      "which under-disperses. Supported families: gaussian, poisson, binomial, Gamma, ",
      "nbinom1, nbinom2, beta, betabinomial, t, lognormal, skewnormal, ",
      "truncated_poisson, truncated_nbinom1, truncated_nbinom2, and tweedie ",
      "(requires the 'tweedie' package).",
      call. = FALSE
    )
    do.call(options, stats::setNames(list(TRUE), key))
  }
  invisible(NULL)
}

# Extract the link-inverse function robustly from a fitted glmmTMB object.
.glmmTMB_linkinv <- function(fitted) {
  fam <- family(fitted)
  if (!is.null(fam$linkinv)) return(fam$linkinv)
  stats::make.link(fam$link)$linkinv
}


# --------------------------------------------------------------------------
# Part 3: fit_model dispatch + glmmTMBmodel constructor
# --------------------------------------------------------------------------

#' Fit a glmmTMB model
#'
#' S3 dispatch for `glmmTMB_spec` objects. Delegates to [glmmTMBmodel()].
#'
#' @param spec A `glmmTMB_spec` object from [build_model()].
#' @param data A data.table or data.frame with training data.
#' @param ctx A [panel_context()] object.
#' @param subset Optional list with `start`/`end` for the training window.
#' @param ... Ignored; accepted for generic consistency.
#'
#' @return A fitted `glmmTMB_endogenr` endogenmodel.
#' @family simulation
#' @export
#' @exportS3Method
fit_model.glmmTMB_spec <- function(spec, data = NULL, ctx = NULL, subset = NULL, ...) {
  glmmTMBmodel(
    formula     = spec$formula,
    family      = if (!is.null(spec$args$family))      spec$args$family      else stats::gaussian(),
    dispformula = if (!is.null(spec$args$dispformula)) spec$args$dispformula else ~1,
    ziformula   = if (!is.null(spec$args$ziformula))   spec$args$ziformula   else ~0,
    control     = spec$args$control,
    data = data, ctx = ctx, subset = subset
  )
}

#' glmmTMB model constructor
#'
#' @param formula Two-sided R formula (lme4-style RE bars allowed).
#' @param family A family object (default `stats::gaussian()`).
#' @param dispformula One-sided dispersion formula (default `~1`).
#' @param ziformula One-sided zero-inflation formula (default `~0`).
#' @param control Optional `glmmTMB::glmmTMBControl()` list.
#' @param data A data.table or data.frame.
#' @param ctx A [panel_context()] object.
#' @param subset Optional list with `start`/`end` for the training window.
#' @param ... Additional arguments (currently unused).
#'
#' @return An endogenmodel of class `glmmTMB_endogenr`.
#' @keywords internal
glmmTMBmodel <- function(formula = NULL,
                         family = stats::gaussian(),
                         dispformula = ~1,
                         ziformula = ~0,
                         control = NULL,
                         data = NULL,
                         ctx = NULL,
                         subset = NULL,
                         ...) {
  if (!requireNamespace("glmmTMB", quietly = TRUE)) {
    stop(
      "Package 'glmmTMB' is required for glmmTMB models. ",
      "Install it with install.packages('glmmTMB').",
      call. = FALSE
    )
  }

  model <- new_endogenmodel(formula)
  model$independent <- FALSE

  grp_keys <- ctx_keys(ctx)
  timevar  <- ctx_time(ctx)
  model$timevar  <- timevar
  model$subset   <- subset
  model$family   <- family

  # Store formula components for graph/closure analysis and predict-time re-use
  model$dispformula <- dispformula
  model$ziformula   <- ziformula

  # Dependency-graph formula (bar-free; built by .glmmTMB_decompose)
  dec <- .glmmTMB_decompose(formula, dispformula, ziformula)
  model$graph_formula <- dec$graph_formula

  # Covariance-structure detection + forecast origin (consumed by the predict
  # method). TRUE if any formula part carries a cov-struct RE term (ar1/ou/exp/…),
  # whose multi-step forecast needs the whole forecast-so-far block in one
  # prediction so glmmTMB applies the correct phi^k correlation decay.
  model$has_covstruct <- {
    cs_parts <- Filter(function(f) inherits(f, "formula"),
                       list(formula, dispformula, ziformula))
    any(unlist(lapply(cs_parts, function(p)
      vapply(.flatten_additive_terms(rlang::f_rhs(p)),
             function(tm) { re <- .glmmTMB_as_re(tm); !is.null(re) && isTRUE(re$covstruct) },
             logical(1)))))
  }
  # Forecast origin: the last time index the fit conditions on. The baseline fit
  # uses the whole training frame (max = test_start - 1); a sliding/min_window
  # refit ends at subset$end.
  model$last_train_time <- if (!is.null(subset)) subset$end
                           else max(data[[timevar]], na.rm = TRUE)

  # ── Stage 1: shared materialization across all formula parts ──────────────
  # Build a combined RHS at the AST level (NOT via terms()/reformulate() which
  # break on RE bars).  The combined formula is never fit; it drives panel_materialize
  # so all lag/roll/cum atoms across main + disp + zi formulas share one ts_map.
  parts <- list(rlang::f_rhs(formula))
  if (!.is_trivial_rhs(dispformula)) parts <- c(parts, list(rlang::f_rhs(dispformula)))
  if (!.is_trivial_rhs(ziformula))   parts <- c(parts, list(rlang::f_rhs(ziformula)))
  combined_rhs     <- Reduce(function(a, b) call("+", a, b), parts)
  combined_formula <- rlang::new_formula(rlang::f_lhs(formula), combined_rhs,
                                         env = rlang::f_env(formula))

  pm        <- panel_materialize(combined_formula, data,
                                 groupvar = grp_keys, timevar = timevar)
  alias_map <- .pt_make_aliases(pm$map)
  old <- intersect(names(alias_map), names(pm$data))
  if (length(old) > 0L) data.table::setnames(pm$data, old, alias_map[old])

  model$ts_map       <- pm$map
  model$pt_alias_map <- alias_map
  model$mat_formula  <- combined_formula  # env used by .apply_ts_map at predict time
  model$groupvar     <- grp_keys

  # Rewrite each formula part through the shared, already-complete map.
  # .rewrite_from_map() uses the existing map so no atom is extracted twice.
  fit_formula <- .pt_alias_formula(.rewrite_from_map(formula,      pm$map), alias_map)
  disp_fit    <- if (.is_trivial_rhs(dispformula)) dispformula
                 else .pt_alias_formula(.rewrite_from_map(dispformula, pm$map), alias_map)
  zi_fit      <- if (.is_trivial_rhs(ziformula)) ziformula
                 else .pt_alias_formula(.rewrite_from_map(ziformula,   pm$map), alias_map)

  model$fit_formula <- fit_formula
  model$disp_fit    <- disp_fit
  model$zi_fit      <- zi_fit

  class(model) <- c("glmmTMB_endogenr", class(model))

  # ── Stage 2: fit ─────────────────────────────────────────────────────────
  # Closure captures pm$data; called once here and optionally again per-draw
  # in refit scenarios (min_window mode).
  model$fit <- function(fit_formula, disp_fit, zi_fit, data, family, control, subset, timevar) {
    # Restrict to columns the model actually reads (keeps grouping vars and
    # cov-struct coordinates because all.vars() walks through bars).
    fit_cols <- unique(c(
      intersect(all.vars(fit_formula), names(data)),
      intersect(all.vars(disp_fit),   names(data)),
      intersect(all.vars(zi_fit),     names(data)),
      timevar
    ))
    data <- data[, ..fit_cols]
    if (!is.null(subset)) {
      data <- data[data[[timevar]] >= subset$start & data[[timevar]] <= subset$end]
    }
    data <- stats::na.omit(as.data.frame(data))

    args <- list(
      formula     = fit_formula,
      data        = data,
      family      = family,
      dispformula = disp_fit,
      ziformula   = zi_fit
    )
    if (!is.null(control)) args$control <- control
    do.call(glmmTMB::glmmTMB, args)
  }

  model$fitted <- model$fit(fit_formula, disp_fit, zi_fit, pm$data,
                            family, control, subset, timevar)

  # Extra family parameters (t df, skewnormal shape, tweedie power, …) for the
  # response draw; numeric(0) for families without any. Survives .strip_fit_data
  # (which only nulls model$data and trims model$fitted sub-fields).
  model$family_params <- tryCatch(glmmTMB::family_params(model$fitted),
                                  error = function(e) numeric(0))

  # Predict-invariant quantities cached once at fit time so the hot predict loop
  # (called nsim * horizon times) avoids redundant work. The fitted model is fixed
  # across all draws/steps, so its family, link-inverse, and — when dispformula is
  # trivial — its dispersion are constant.
  model$fam_name <- family(model$fitted)$family
  model$linkinv  <- .glmmTMB_linkinv(model$fitted)
  # Constant-dispersion fast path: when dispformula is trivial (~1 / ~0),
  # predict(type = "disp") returns sigma(fitted) for every row (verified across
  # gaussian/Gamma/nbinom2/poisson/beta), so cache the scalar and skip a full
  # predict.glmmTMB rebuild at every forecast step. NULL signals the predict
  # method to fall back to the per-row predict path (non-trivial dispformula).
  model$disp_const <- if (.is_trivial_rhs(dispformula)) {
    tryCatch(as.numeric(stats::sigma(model$fitted)), error = function(e) NULL)
  } else {
    NULL
  }

  # ── Coefficient + GOF fields for plot_coefficients() ─────────────────────
  co <- tryCatch(summary(model$fitted)$coefficients$cond, error = function(e) NULL)
  model$coefs <- if (!is.null(co) && nrow(co) > 0L) {
    data.frame(
      term      = rownames(co),
      estimate  = co[, 1L],
      std.error = co[, 2L],
      statistic = co[, 3L],
      p.value   = co[, 4L],
      row.names = NULL,
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(term = character(0), estimate = numeric(0), std.error = numeric(0),
               statistic = numeric(0), p.value = numeric(0), stringsAsFactors = FALSE)
  }

  model$gof <- tryCatch(
    data.frame(
      nobs        = stats::nobs(model$fitted),
      logLik      = as.numeric(stats::logLik(model$fitted)),
      AIC         = stats::AIC(model$fitted),
      BIC         = stats::BIC(model$fitted),
      df.residual = stats::df.residual(model$fitted)
    ),
    error = function(e) data.frame(nobs = NA_integer_, logLik = NA_real_,
                                   AIC = NA_real_, BIC = NA_real_, df.residual = NA_real_)
  )

  model$outcome <- dec$outcome

  # Required history: max over all formula parts (.required_history descends
  # through | and ( as generic calls, so it correctly finds lag() inside RE bars).
  model$required_history <- max(
    .required_history(formula),
    .required_history(dispformula),
    .required_history(ziformula)
  )

  model$data <- pm$data  # stripped by .strip_fit_data() after fit_system()

  return(model)
}


# --------------------------------------------------------------------------
# Part 4: predict.glmmTMB_endogenr
# --------------------------------------------------------------------------

#' Predict method for glmmTMB models
#'
#' Called by the endogenr dynamic simulation loop at each forecast time step.
#' Re-materialises time-series columns, predicts on the link scale with
#' parameter uncertainty, adds family-specific response noise, and optionally
#' applies a structural-zero Bernoulli mask for zero-inflated models.
#'
#' @param model A `glmmTMB_endogenr` endogenmodel.
#' @param data A data.table (the full simulation grid with history rows).
#' @param t Integer. Time step to predict.
#' @param ctx A [panel_context()] object.
#' @param what Either `"pi"` (predictive interval draw) or `"expectation"`.
#' @param ... Ignored; accepted for S3 generic consistency.
#'
#' @return A data.table with columns `c(ctx_keys, ctx_time, outcome)`,
#'   one row per `(unit, sim)` at time `t`.
#' @family simulation
#' @export
predict.glmmTMB_endogenr <- function(model, data, t, ctx, what = "pi", ...) {
  idx      <- ctx_time(ctx)
  all_keys <- ctx_keys(ctx)

  # Covariance-structure terms (ar1, ou, exp, …) need the WHOLE forecast-so-far
  # block in one prediction so glmmTMB applies the correct multi-step correlation
  # decay (phi^k) and growing variance. Pass rows [last_train_time + 1 .. t]
  # (plus required_history rows for lag materialization) and keep only the rows
  # at t. Plain-RE / fixed-effect models predict only the rows at t.
  if (isTRUE(model$has_covstruct)) {
    lo   <- model$last_train_time + 1L - model$required_history
    data <- data[data[[idx]] >= lo & data[[idx]] <= t]
  } else {
    data <- .history_subset(data, idx, t, model$required_history)
  }

  # Re-materialise ts columns per (unit, sim) group
  env <- rlang::f_env(model$mat_formula)
  mat <- .apply_ts_map(model$ts_map, data, all_keys, idx, env = env, copy = FALSE)
  old <- intersect(names(model$pt_alias_map), names(mat))
  if (length(old) > 0L) data.table::setnames(mat, old, model$pt_alias_map[old])

  # Block of rows to predict over; only the rows at t are written back.
  if (isTRUE(model$has_covstruct)) {
    mat <- mat[mat[[idx]] >= model$last_train_time + 1L & mat[[idx]] <= t]
  } else {
    mat <- mat[mat[[idx]] == t]
  }
  at_t <- mat[[idx]] == t
  df   <- as.data.frame(mat)

  result_cols <- c(all_keys, idx, model$outcome)
  result      <- mat[at_t, ..result_cols]
  if (nrow(mat) == 0L) return(result)

  if (what == "expectation") {
    val <- stats::predict(model$fitted, newdata = df, type = "response",
                          re.form = NULL, allow.new.levels = TRUE)
    data.table::set(result, j = model$outcome, value = as.numeric(val[at_t]))

  } else if (what == "pi") {
    # Link-scale prediction with SE for parameter uncertainty
    pl  <- stats::predict(model$fitted, newdata = df, type = "link", se.fit = TRUE,
                          re.form = NULL, allow.new.levels = TRUE)
    # Asymptotic normal draw on link scale (standard for ML mixed models)
    eta <- pl$fit + pl$se.fit * stats::rnorm(length(pl$fit))
    mu  <- model$linkinv(eta)

    # Per-row dispersion. Trivial dispformula -> constant sigma cached at fit time
    # (skips a full predict.glmmTMB rebuild each step). Non-trivial dispformula ->
    # per-row predict as before.
    disp <- if (!is.null(model$disp_const)) {
      rep(model$disp_const, length(mu))
    } else {
      tryCatch(
        as.numeric(stats::predict(model$fitted, newdata = df, type = "disp",
                                  allow.new.levels = TRUE)),
        error = function(e) rep(1, length(mu))
      )
    }

    fam_name <- model$fam_name
    draw <- .glmmTMB_response_draw(mu, fam_name, disp, model$family_params)

    # Structural-zero mask for zero-inflated models
    if (!.is_trivial_rhs(model$ziformula)) {
      p_zero <- tryCatch(
        as.numeric(stats::predict(model$fitted, newdata = df, type = "zprob",
                                  allow.new.levels = TRUE)),
        error = function(e) rep(0, length(draw))
      )
      p_zero   <- pmin(pmax(p_zero, 0), 1)
      not_zero <- stats::rbinom(length(draw), size = 1L, prob = 1 - p_zero)
      draw     <- draw * not_zero
    }

    data.table::set(result, j = model$outcome, value = as.numeric(draw[at_t]))

  } else {
    stop("`what` must be either `pi` or `expectation`", call. = FALSE)
  }

  return(result)
}
