# gamlss model support for endogenr
#
# Integrates GAMLSS (Generalized Additive Models for Location, Scale and Shape)
# as step-by-step dependent models inside the endogenr dynamic simulation loop.
# The model is predicted one time step at a time via gamlss::predictAll(newdata=…,
# data=…), exactly like glm/glmmTMB.
#
# Two-stage approach mirroring glmmTMBmodel.R:
#   Stage 1: shared AST-level materialization of lag/roll/cum sub-expressions
#            across all four parameter formulas (mu, sigma, nu, tau) into .pt#
#            columns via panel_materialize.
#   Stage 2: pooled gamlss::gamlss fit on the rewritten aliased formulas.
#
# IMPORTANT:
#   - fitted class is "endogenr_gamlss" (NOT "gamlss_endogenr" or "gamlss") to
#     avoid clobbering gamlss's own S3 methods AND to avoid colliding with the
#     regex `regexpr("predict.gamlss", sys.calls())` that gamlss.pb uses to find
#     the predict.gamlss frame (a method named predict.gamlss* would hijack it).
#   - predictAll MUST receive data = model$gamlss_data (the exact training frame
#     used during the fit). Without it, predictAll errors with "object not found"
#     when called from inside a function scope; a different frame gives silently
#     wrong predictions.

# --------------------------------------------------------------------------
# Part 1: Formula decomposition helpers (pure base-R; no gamlss needed)
# --------------------------------------------------------------------------

# Extract predictor sub-terms from one additive gamlss formula term.
#
# Returns a list of language objects (each becomes a dependency-graph edge).
# Grouping factors from random()/ra()/re() are ordinary predictors: they are
# read at every forecast step, so — like any predictor — they must be produced
# by some model (e.g. an `exogen` carrying the column forward).
#   - random(g) / ra(g): the grouping expression g is a predictor.
#   - re(...): scan formula arguments; both bar sides (effect | group) and any
#     non-bar formula RHS contribute predictors; non-formula args are ignored.
#   - everything else (fixed, pb/cs/lo smoother wrapping a variable): the term
#     itself is a predictor when it carries a variable.
.gamlss_classify_term <- function(term) {
  preds <- list()
  add <- function(expr) if (length(all.vars(expr)) > 0L) preds[[length(preds) + 1L]] <<- expr

  if (!is.call(term)) {
    add(term)
    return(preds)
  }

  fn <- as.character(term[[1L]])

  # random(g) or ra(g) — first argument is the grouping factor (a predictor)
  if (fn == "random" || fn == "ra") {
    add(term[[2L]])
    return(preds)
  }

  # re(...) — scan all arguments for formula-looking calls (stored as `call`
  # objects in the AST, not as S3 class "formula"; detect via [[1L]] == "~")
  if (fn == "re") {
    args <- as.list(term)[-1L]  # drop the function name
    for (a in args) {
      is_tilde_call <- is.language(a) && is.call(a) &&
                       identical(a[[1L]], as.name("~"))
      if (is_tilde_call) {
        rhs <- a[[length(a)]]  # last element: RHS for both ~rhs and lhs~rhs
        if (.glmmTMB_is_bar(rhs)) {
          # bar LHS = effect, bar RHS = grouping — both are predictors
          for (et in .flatten_additive_terms(rhs[[2L]])) add(et)
          add(rhs[[3L]])
        } else {
          for (et in .flatten_additive_terms(rhs)) add(et)
        }
      }
      # Non-tilde args (method=, control=, …) are ignored
    }
    return(preds)
  }

  # Default: fixed term or smoother (pb, cs, lo, I, poly, …)
  add(term)
  preds
}


#' Decompose gamlss formulas into a smoother/grouping-free graph formula
#'
#' Extracts every predictor variable (for dependency-graph edges) from the four
#' parameter formulas (mu/sigma/nu/tau). Grouping factors from random()/ra()/re()
#' are treated as ordinary predictors: they are read at every forecast step, so
#' they must be produced by some model (e.g. an `exogen` carrying the grouping
#' column forward). Returns a graph_formula that `stats::terms()` and
#' `.edges_from_formula()` can consume safely (no bars, no random() calls).
#'
#' @param formula Main two-sided gamlss formula (mu; may contain pb/cs/lo,
#'   random/ra/re).
#' @param sigma.formula Optional one-sided sigma formula.
#' @param nu.formula Optional one-sided nu formula.
#' @param tau.formula Optional one-sided tau formula.
#'
#' @return A list with:
#'   \item{outcome}{Character: the LHS variable name.}
#'   \item{graph_formula}{A formula suitable for dependency analysis.}
#' @keywords internal
.gamlss_decompose <- function(formula, sigma.formula = NULL,
                              nu.formula = NULL, tau.formula = NULL) {
  outcome    <- all.vars(rlang::f_lhs(formula))
  pred_terms <- list()

  parts <- Filter(function(f) inherits(f, "formula"),
                  list(formula, sigma.formula, nu.formula, tau.formula))

  for (part in parts) {
    for (term in .flatten_additive_terms(rlang::f_rhs(part))) {
      pred_terms <- c(pred_terms, .gamlss_classify_term(term))
    }
  }

  # Deduplicate pred_terms by deparse, preserving first occurrence
  seen         <- character(0)
  unique_terms <- list()
  for (t in pred_terms) {
    key <- deparse(t, width.cutoff = 500L)
    if (!key %in% seen) {
      seen         <- c(seen, key)
      unique_terms <- c(unique_terms, list(t))
    }
  }
  pred_terms <- unique_terms

  # Build the smoother/grouping-free graph RHS
  graph_rhs <- if (length(pred_terms) == 0L) {
    1
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
# Part 2: Family response-draw helper (no gamlss needed at definition time)
# --------------------------------------------------------------------------

# One-time per-session warning when the family has no r<FAM> generator.
.gamlss_warn_no_rfun <- function(family_name) {
  key <- paste0(".endogenr_gamlss_warned_rfun_", family_name)
  if (!isTRUE(getOption(key))) {
    warning(
      "gamlss family '", family_name, "' has no r", family_name, "() generator ",
      "in the gamlss.dist namespace; falling back to the conditional mean (mu). ",
      "Install the relevant gamlss.dist extension or use a standard gamlss.dist family.",
      call. = FALSE
    )
    do.call(options, stats::setNames(list(TRUE), key))
  }
  invisible(NULL)
}


# Predictive draw from the gamlss family's random generator.
#
# `fitted` — a gamlss fitted object.
# `pa`     — named list from predictAll(..., type = "response").
# `n`      — number of draws (= nrow(newdata)).
#
# Looks up r<FAM> from gamlss.dist (or the calling environment) and calls it
# with named arguments built from fitted$parameters. Falls back to pa$mu with
# a one-time warning for families without an r function.
.gamlss_response_draw <- function(fitted, pa, n) {
  fam_short <- fitted$family[1L]

  rfun <- tryCatch(
    getExportedValue("gamlss.dist", paste0("r", fam_short)),
    error = function(e) get0(paste0("r", fam_short), mode = "function")
  )

  if (is.null(rfun)) {
    .gamlss_warn_no_rfun(fam_short)
    return(as.numeric(pa$mu))
  }

  # Build args by name from the active parameter set.
  # bd (trial denominator for binomial-type families) defaults to 1 if omitted.
  args <- list(n = n)
  for (p in fitted$parameters) args[[p]] <- as.numeric(pa[[p]])

  as.numeric(do.call(rfun, args))
}


# --------------------------------------------------------------------------
# Part 3: fit_model dispatch + gamlssmodel constructor
# --------------------------------------------------------------------------

#' Stage-2 pooled gamlss fit on materialized data
#'
#' File-level analogue of `.lm_stage2_fit()` — see that helper for why the
#' fit step must not be a closure stored on the model object.
#'
#' @param mu_fit,sigma_fit,nu_fit,tau_fit Rewritten, aliased parameter formulas.
#' @param data The materialized data.table.
#' @param family A `gamlss.family` object.
#' @param ctrl A `gamlss::gamlss.control()` list.
#' @param subset Optional list with `start`/`end` training window.
#' @param timevar Character. Time column name.
#' @return A list with `fitted` (the gamlss fit) and `frame` (the EXACT
#'   data.frame the fit saw — required by `predictAll(data = ...)`).
#' @keywords internal
.gamlss_stage2_fit <- function(mu_fit, sigma_fit, nu_fit, tau_fit, data,
                               family, ctrl, subset, timevar) {
  fit_cols <- unique(c(
    intersect(all.vars(mu_fit),    names(data)),
    intersect(all.vars(sigma_fit), names(data)),
    intersect(all.vars(nu_fit),    names(data)),
    intersect(all.vars(tau_fit),   names(data)),
    timevar
  ))
  d <- data[, ..fit_cols]
  if (!is.null(subset)) {
    d <- .dt_rows(d, d[[timevar]] >= subset$start & d[[timevar]] <= subset$end)
  }
  d <- stats::na.omit(as.data.frame(d))
  # Ensure gamlss functions (pb, cs, lo, random, ra, re, etc.) are findable
  # when gamlss evaluates formula terms in model.frame. Required because
  # requireNamespace() loads but does not attach the package.
  gamlss_ns <- asNamespace("gamlss")
  .set_gamlss_env <- function(f) {
    environment(f) <- new.env(parent = gamlss_ns)
    f
  }
  mu_fit    <- .set_gamlss_env(mu_fit)
  sigma_fit <- .set_gamlss_env(sigma_fit)
  nu_fit    <- .set_gamlss_env(nu_fit)
  tau_fit   <- .set_gamlss_env(tau_fit)
  fitted <- gamlss::gamlss(
    formula       = mu_fit,
    sigma.formula = sigma_fit,
    nu.formula    = nu_fit,
    tau.formula   = tau_fit,
    family        = family,
    data          = d,
    control       = ctrl
  )
  list(fitted = fitted, frame = d)
}

#' Fit a gamlss model
#'
#' S3 dispatch for `gamlss_spec` objects. Delegates to [gamlssmodel()].
#'
#' @param spec A `gamlss_spec` object from [build_model()].
#' @param data A data.table or data.frame with training data.
#' @param ctx A [panel_context()] object.
#' @param subset Optional list with `start`/`end` for the training window.
#' @param ... Ignored; accepted for generic consistency.
#'
#' @return A fitted `endogenr_gamlss` endogenmodel.
#' @family simulation
#' @export
#' @exportS3Method
fit_model.gamlss_spec <- function(spec, data = NULL, ctx = NULL, subset = NULL, ...) {
  gamlssmodel(
    formula       = spec$formula,
    sigma.formula = if (!is.null(spec$args[["sigma.formula"]])) spec$args[["sigma.formula"]] else ~1,
    nu.formula    = if (!is.null(spec$args[["nu.formula"]]))    spec$args[["nu.formula"]]    else ~1,
    tau.formula   = if (!is.null(spec$args[["tau.formula"]]))   spec$args[["tau.formula"]]   else ~1,
    family        = if (!is.null(spec$args[["family"]]))        spec$args[["family"]]        else gamlss.dist::NO(),
    control       = spec$args[["control"]],
    data = data, ctx = ctx, subset = subset
  )
}


#' gamlss model constructor
#'
#' @param formula Two-sided R formula for the mu (location) parameter.
#'   May include gamlss smoothers (`pb()`, `cs()`, `lo()`) and grouping terms
#'   (`random()`, `ra()`, `re()`).
#' @param sigma.formula One-sided formula for the sigma (scale) parameter
#'   (default `~1`).
#' @param nu.formula One-sided formula for the nu (shape 1) parameter
#'   (default `~1`).
#' @param tau.formula One-sided formula for the tau (shape 2) parameter
#'   (default `~1`).
#' @param family A `gamlss.family` object (default `gamlss.dist::NO()`).
#' @param control Optional `gamlss::gamlss.control()` list.
#' @param data A data.table or data.frame.
#' @param ctx A [panel_context()] object.
#' @param subset Optional list with `start`/`end` for the training window.
#' @param ... Additional arguments (currently unused).
#'
#' @return An endogenmodel of class `endogenr_gamlss`.
#' @keywords internal
gamlssmodel <- function(formula       = NULL,
                        sigma.formula = ~1,
                        nu.formula    = ~1,
                        tau.formula   = ~1,
                        family        = gamlss.dist::NO(),
                        control       = NULL,
                        data          = NULL,
                        ctx           = NULL,
                        subset        = NULL,
                        ...) {
  if (!requireNamespace("gamlss", quietly = TRUE)) {
    stop(
      "Package 'gamlss' is required for gamlss models. ",
      "Install it with install.packages('gamlss').",
      call. = FALSE
    )
  }

  model <- new_endogenmodel(formula)
  model$independent <- FALSE

  grp_keys       <- ctx_keys(ctx)
  timevar        <- ctx_time(ctx)
  model$timevar  <- timevar
  model$subset   <- subset
  model$family   <- family

  # Store parameter formulas under underscore names to avoid $ partial matching
  # against "sigma", "nu", "tau" when the model list is accessed.
  model$sigma_formula <- sigma.formula
  model$nu_formula    <- nu.formula
  model$tau_formula   <- tau.formula

  # Dependency-graph formula (smoother/grouping-free; built by .gamlss_decompose)
  dec                 <- .gamlss_decompose(formula, sigma.formula, nu.formula, tau.formula)
  model$graph_formula <- dec$graph_formula
  model$outcome       <- dec$outcome

  # ── Stage 1: shared materialization across all parameter formula parts ──────
  # Build a combined RHS at the AST level (NOT via terms()/reformulate() which
  # break on random()/pb() calls).  The combined formula is never fit; it drives
  # panel_materialize so all lag/roll/cum atoms share one ts_map.
  parts <- list(rlang::f_rhs(formula))
  if (!.is_trivial_rhs(sigma.formula)) parts <- c(parts, list(rlang::f_rhs(sigma.formula)))
  if (!.is_trivial_rhs(nu.formula))    parts <- c(parts, list(rlang::f_rhs(nu.formula)))
  if (!.is_trivial_rhs(tau.formula))   parts <- c(parts, list(rlang::f_rhs(tau.formula)))

  combined_rhs     <- Reduce(function(a, b) call("+", a, b), parts)
  combined_formula <- rlang::new_formula(rlang::f_lhs(formula), combined_rhs,
                                         env = rlang::f_env(formula))

  pm        <- panel_materialize(combined_formula, data,
                                 groupvar = grp_keys, timevar = timevar)
  alias_map <- .pt_make_aliases(pm$map)
  .pt_apply_aliases(pm$data, alias_map)

  model$ts_map       <- pm$map
  model$pt_alias_map <- alias_map
  model$mat_formula  <- combined_formula  # env used by .apply_ts_map at predict time
  model$groupvar     <- grp_keys

  # Rewrite each parameter formula through the shared, already-complete map.
  mu_fit    <- .pt_alias_formula(.rewrite_from_map(formula,       pm$map), alias_map)
  sigma_fit <- if (.is_trivial_rhs(sigma.formula)) sigma.formula else
                 .pt_alias_formula(.rewrite_from_map(sigma.formula, pm$map), alias_map)
  nu_fit    <- if (.is_trivial_rhs(nu.formula)) nu.formula else
                 .pt_alias_formula(.rewrite_from_map(nu.formula,    pm$map), alias_map)
  tau_fit   <- if (.is_trivial_rhs(tau.formula)) tau.formula else
                 .pt_alias_formula(.rewrite_from_map(tau.formula,   pm$map), alias_map)

  class(model) <- c("endogenr_gamlss", class(model))

  # ── Stage 2: fit ─────────────────────────────────────────────────────────
  # File-level helper — see .lm_stage2_fit for why this must not be a closure
  # stored on the model object (serialization payload).
  #
  # CRITICAL: store the exact fit frame so predictAll() can always receive
  # data = model$gamlss_data. Without this, predictAll errors or gives wrong
  # predictions when called from inside a function scope (name lookup fails).
  ctrl <- if (!is.null(control)) control else gamlss::gamlss.control(trace = FALSE)

  fitres            <- .gamlss_stage2_fit(mu_fit, sigma_fit, nu_fit, tau_fit,
                                          pm$data, family, ctrl, subset, timevar)
  model$fitted      <- fitres$fitted
  model$gamlss_data <- fitres$frame  # EXACT frame; never replaced; passed as data= to predictAll

  # ── mu fixed-effect detection (TWFE / TW-Mundlak) ─────────────────────────
  # Detect factor(timevar) and factor(unitkey) in the mu formula only.
  # mu.xlevels is the gamlss analogue of lm$xlevels.
  mu_coefs       <- stats::coef(model$fitted)   # mu coefficients (named)
  model$time_fe  <- tryCatch(
    .detect_time_fe(mu_fit, mu_coefs, model$fitted$mu.xlevels, timevar),
    error = function(e) NULL
  )
  model$unit_fe  <- tryCatch(
    .detect_unit_fe(mu_fit, mu_coefs, model$fitted$mu.xlevels, grp_keys),
    error = function(e) NULL
  )

  # ── Coefficient + GOF fields for plot_coefficients() ─────────────────────
  # Best-effort only; NULL-safe because get_coefficients() skips NULL $coefs.
  model$coefs <- tryCatch({
    co <- stats::coef(model$fitted)   # mu (location) coefficients
    data.frame(
      term      = names(co),
      estimate  = as.numeric(co),
      std.error = NA_real_,
      statistic = NA_real_,
      p.value   = NA_real_,
      row.names = NULL,
      stringsAsFactors = FALSE
    )
  }, error = function(e) NULL)

  model$gof <- tryCatch(
    data.frame(
      nobs        = model$fitted$N,
      df.fit      = model$fitted$df.fit,
      df.residual = model$fitted$df.residual,
      aic         = model$fitted$aic,
      sbc         = model$fitted$sbc,
      deviance    = model$fitted$G.deviance
    ),
    error = function(e) data.frame(nobs = NA_integer_)
  )

  # Required history: max over all four parameter formulas.
  # .required_history descends through pb()/random()/ra() as generic calls and
  # correctly finds lag() inside them.
  model$required_history <- max(
    .required_history(formula),
    .required_history(sigma.formula),
    .required_history(nu.formula),
    .required_history(tau.formula)
  )

  model$data <- pm$data  # stripped by .strip_fit_data() after fit_system()

  return(model)
}


# --------------------------------------------------------------------------
# Part 4: predict.endogenr_gamlss
# --------------------------------------------------------------------------

#' Predict method for gamlss models
#'
#' Called by the endogenr dynamic simulation loop at each forecast time step.
#' Re-materialises time-series columns, calls [gamlss::predictAll()] to get
#' all distribution parameters, then either returns the location mean
#' (`what = "expectation"`) or draws from the family's random generator
#' (`what = "pi"`).
#'
#' @param model An `endogenr_gamlss` endogenmodel.
#' @param data A data.table (the full simulation grid with history rows).
#' @param t Integer. Time step to predict.
#' @param ctx A [panel_context()] object.
#' @param what Either `"pi"` (predictive interval draw) or `"expectation"`
#'   (returns the location parameter mu).
#' @param ... Ignored; accepted for S3 generic consistency.
#'
#' @details
#'   Prediction rows whose predictors are `NA` at step `t` (e.g. a lagged
#'   exogenous predictor reaching a period the exogenous series does not cover)
#'   yield `NA` for that row rather than erroring, matching
#'   `predict.lm`/`predict.glm`. Supply every predictor across the full
#'   forecast window `[test_start, test_start + horizon - 1]` to avoid `NA`
#'   (and, under `bounds`, midpoint-filled) forecasts.
#'
#' @return A data.table with columns `c(ctx_keys, ctx_time, outcome)`,
#'   one row per `(unit, sim)` at time `t`.
#' @family simulation
#' @export
predict.endogenr_gamlss <- function(model, data, t, ctx, what = "pi",
                                    scenario = NULL, test_start = NULL, ...) {
  idx      <- ctx_time(ctx)
  all_keys <- ctx_keys(ctx)

  # Subset to the per-unit history window the RHS actually needs
  data <- .history_subset(data, idx, t, model$required_history)

  # Re-materialise ts columns per (unit, sim) group
  env <- rlang::f_env(model$mat_formula)
  mat <- .apply_ts_map(model$ts_map, data, all_keys, idx, env = env, copy = FALSE)
  .pt_apply_aliases(mat, model$pt_alias_map)

  # Filter to the prediction time step
  mat <- .dt_rows(mat, mat[[idx]] == t)

  result_cols <- c(all_keys, idx, model$outcome)
  result      <- mat[, ..result_cols]

  if (nrow(mat) == 0L) return(result)

  # --- Scenario: extract per-draw slices ------------------------------------
  oc_scen       <- if (!is.null(scenario)) scenario[[model$outcome]] else NULL
  time_fe_slice <- if (!is.null(oc_scen)) oc_scen$time_fe else NULL
  unit_fe_slice <- if (!is.null(oc_scen)) oc_scen$unit_fe else NULL
  eff_coef      <- if (!is.null(oc_scen) && !is.null(oc_scen$coef))
                     oc_scen$coef$effective
                   else list()
  h             <- if (!is.null(test_start)) t - test_start + 1L else NULL

  # --- Neutralise FE columns on a copy of mat --------------------------------
  # When a time or unit FE slice is active, replace the factor column with its
  # reference level so predictAll() contributes 0 for that contrast; the baked
  # eta offset is added below (on the link scale).
  pred_frame <- mat

  if (!is.null(time_fe_slice)) {
    pred_frame <- data.table::copy(mat)
    data.table::set(pred_frame, j = time_fe_slice$var, value = time_fe_slice$ref)
  }

  if (!is.null(unit_fe_slice)) {
    if (identical(pred_frame, mat)) pred_frame <- data.table::copy(mat)
    ref_val   <- unit_fe_slice$ref
    col_class <- class(mat[[unit_fe_slice$var]])[1L]
    if (col_class == "integer")       ref_val <- as.integer(ref_val)
    else if (col_class == "numeric")  ref_val <- as.numeric(ref_val)
    data.table::set(pred_frame, j = unit_fe_slice$var, value = ref_val)
  }

  # --- Compute per-row eta offset (link scale) --------------------------------
  # Non-NULL only when h is known and at least one scenario channel is active.
  has_active_scenario <- !is.null(h) &&
    (!is.null(time_fe_slice) || !is.null(unit_fe_slice) || length(eff_coef) > 0L)

  eta_off <- if (has_active_scenario) {
    off <- numeric(nrow(pred_frame))

    # sim column — needed to index per-trajectory FE values
    sim_col <- ctx_sim(ctx)
    sim_ids <- if (!is.null(sim_col) && sim_col %in% names(pred_frame))
                 pred_frame[[sim_col]]
               else rep(1L, nrow(pred_frame))

    # Time FE offset (baked on eta / link scale)
    if (!is.null(time_fe_slice)) {
      if (what == "pi") {
        off <- off + time_fe_slice$offsets[sim_ids, h]
      } else {
        off <- off + mean(time_fe_slice$offsets[, h])
      }
    }

    # Unit FE offset
    if (!is.null(unit_fe_slice)) {
      units_here <- as.character(mat[[unit_fe_slice$var]])
      unit_idx   <- match(units_here, rownames(unit_fe_slice$eff))
      if (what == "pi" && isTRUE(unit_fe_slice$per_trajectory)) {
        off <- off + unit_fe_slice$eff[cbind(unit_idx, sim_ids,
                                             rep(h, length(unit_idx)))]
      } else {
        off <- off + unit_fe_slice$eff[unit_idx, h]
      }
    }

    # Coef overrides: delta on the mu linear predictor
    if (length(eff_coef) > 0L) {
      X  <- .gamlss_mu_designmatrix(model$fitted, pred_frame)
      bh <- stats::coef(model$fitted)
      for (nm in names(eff_coef)) {
        if (!nm %in% colnames(X))
          stop("coef override term '", nm, "' has no mu design column", call. = FALSE)
        off <- off + (eff_coef[[nm]][[h]] - bh[[nm]]) * X[, nm]
      }
    }

    off
  } else NULL

  # Distribution-parameter prediction and draws happen inside draw_predictive;
  # the family carries no parameter-uncertainty draw, so n_param = 0.
  if (what == "expectation") {
    data.table::set(result, j = model$outcome,
                    value = as.vector(
                      draw_predictive(model, pred_frame, n_param = 0L, n_innov = 0L,
                                      eta_offset = eta_off)))
  } else if (what == "pi") {
    data.table::set(result, j = model$outcome,
                    value = as.vector(
                      draw_predictive(model, pred_frame, n_param = 0L, n_innov = 1L,
                                      eta_offset = eta_off)))
  } else {
    stop("`what` must be either `pi` or `expectation`", call. = FALSE)
  }

  return(result)
}

#' @rdname draw_predictive
#' @export
draw_predictive.endogenr_gamlss <- function(model, newdata, n_param = 1L, n_innov = 1L,
                                            eta_offset = NULL, ...) {
  .check_draw_counts(n_param, n_innov)
  if (n_param > 0L) .warn_no_param_draw("gamlss")

  # newdata for predictAll: predictor/grouping columns of the fit frame, WITHOUT
  # the outcome (predictAll only needs the predictors; outcome at t is NA anyway).
  newdata <- data.table::as.data.table(newdata)
  keep <- intersect(setdiff(names(model$gamlss_data), model$outcome), names(newdata))
  df   <- as.data.frame(newdata[, ..keep])

  n   <- nrow(df)
  K   <- max(n_innov, 1L)
  out <- matrix(NA_real_, n, K)

  # gamlss::predictAll() row-binds the training frame with `newdata`, then indexes
  # the pooled design matrix with a logical vector built BEFORE model.frame() drops
  # NA rows — so any NA predictor row makes that logical longer than the matrix and
  # crashes with "(subscript) logical subscript too long". Predict only the
  # complete-case rows (mirroring predict.lm's na.pass) and leave NA for the rest.
  pred_cols <- setdiff(names(df), model$timevar)
  ok <- if (length(pred_cols) == 0L) rep(TRUE, n)
        else stats::complete.cases(df[, pred_cols, drop = FALSE])
  if (!any(ok)) return(out)

  df_ok <- df[ok, , drop = FALSE]
  n_ok  <- sum(ok)

  # predictAll requires data = the EXACT training frame used during the fit.
  # suppressMessages() silences chatty random()-related output.
  pa <- suppressMessages(
    gamlss::predictAll(model$fitted, newdata = df_ok,
                       data = model$gamlss_data, type = "response")
  )

  # Apply eta offset on the mu link scale when a scenario is active.
  # eta_offset is length nrow(newdata); eta_offset[ok] aligns with df_ok rows.
  if (!is.null(eta_offset)) {
    lf    <- model$family$mu.linkfun
    li    <- model$family$mu.linkinv
    pa$mu <- li(lf(as.numeric(pa$mu)) + eta_offset[ok])
  }

  for (k in seq_len(K)) {
    out[ok, k] <- if (n_innov == 0L) as.numeric(pa$mu)
                  else as.numeric(.gamlss_response_draw(model$fitted, pa, n_ok))
  }
  out
}

# Helper: reconstruct the mu design matrix for newdata using the fitted mu terms.
# Used by predict.endogenr_gamlss for coefficient-override delta computation.
.gamlss_mu_designmatrix <- function(fitted, newdata) {
  tt <- stats::delete.response(fitted$mu.terms)
  mf <- stats::model.frame(tt, data = as.data.frame(newdata),
                            xlev = fitted$mu.xlevels)
  stats::model.matrix(tt, mf)
}
