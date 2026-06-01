# Helper: substitute a symbol name in an R expression (recursive AST walk)
.sub_sym <- function(expr, from, to) {
  if (is.symbol(expr)) {
    if (identical(as.character(expr), from)) return(as.symbol(to))
    return(expr)
  }
  if (is.call(expr)) {
    return(as.call(lapply(as.list(expr), .sub_sym, from = from, to = to)))
  }
  expr
}

# Helper: clean a single variable name using janitor conventions
.clean_name <- function(x) {
  names(janitor::clean_names(stats::setNames(as.data.frame(list(1L)), x)))
}


#' Lead a vector by a forecast horizon (positional within-group shift)
#'
#' Returns the value `h` steps ahead within an ordered vector, padding the tail
#' with `NA`. This is the long-horizon analogue of the positional `lag()` used in
#' simulation formulas, and the marker recognised on the left-hand side of
#' long-horizon formulas (see [setup_long_horizon()]).
#'
#' Inside a long-horizon formula, write the outcome on the LHS wrapped in
#' `lead_horizon()`, e.g. `lead_horizon(gdppc_grwt) ~ gdppc_grwt + log(gdppc)`.
#' The horizon `h` is supplied internally per horizon by [setup_long_horizon()];
#' you do not pass it in the formula. The right-hand side is evaluated at the
#' forecast origin (the last observed period), so write `lag()` only when you
#' deliberately want older history.
#'
#' @param x A vector, assumed sorted in time order within a group.
#' @param h Integer. Number of steps to lead.
#'
#' @return A vector the same length and type as `x`, led by `h` with `NA` padding.
#' @family long_horizon
#' @export
#'
#' @examples
#' lead_horizon(1:5, 2)
lead_horizon <- function(x, h = 1L) {
  n <- length(x)
  if (h <= 0L) return(x)
  if (h >= n) return(x[rep(NA_integer_, n)])
  c(x[(h + 1L):n], x[rep(NA_integer_, h)])
}

# Helper: extract the inner expression EXPR from a `lead_horizon(EXPR)` LHS.
# Errors if the formula LHS is not wrapped in lead_horizon().
.lh_lhs_inner <- function(formula) {
  lhs <- rlang::f_lhs(formula)
  if (is.null(lhs) || !is.call(lhs) ||
      !identical(lhs[[1L]], as.name("lead_horizon"))) {
    stop("Long-horizon formula LHS must be wrapped in lead_horizon(), e.g. ",
         "lead_horizon(gdppc_grwt) ~ gdppc_grwt + log(gdppc). Got: ",
         deparse(lhs), call. = FALSE)
  }
  lhs[[2L]]
}

# Helper: the raw outcome variable name of a long-horizon formula.
.lh_outcome <- function(formula) {
  all.vars(.lh_lhs_inner(formula))[1L]
}

# Helper: normalise `formulas` into a uniquely named list of formulas. Accepts a
# bare formula or a list, and fills blank/missing names with `model1`,
# `model2`, ... so every variant has a stable identifier for the `variant`
# column (an unnamed list would otherwise be skipped silently).
.name_formulas <- function(formulas) {
  if (inherits(formulas, "formula")) formulas <- list(formulas)
  if (!is.list(formulas) || length(formulas) == 0L) {
    stop("`formulas` must be a formula or a non-empty list of formulas.",
         call. = FALSE)
  }
  nm <- names(formulas)
  if (is.null(nm)) nm <- rep("", length(formulas))
  blank <- is.na(nm) | !nzchar(nm)
  nm[blank] <- paste0("model", seq_along(formulas))[blank]
  names(formulas) <- make.unique(nm)
  formulas
}


#' Create aligned (t, t+h) training data for long-horizon models
#'
#' For each group, computes the h-step-ahead value of the outcome (via
#' [lead_horizon()]) and stores it in a new `.target` column. The original
#' outcome column is left untouched so that formula RHS terms are evaluated at
#' the baseline period t, not the target period t+h.
#'
#' @param data A data.frame or data.table containing all variables.
#' @param outcome Character. The raw outcome variable name (without transformation).
#' @param h Integer. The forecast horizon.
#' @param groupvar Character. The group variable name.
#' @param timevar Character. The time variable name.
#' @param test_start Numeric. First forecast year; only rows whose target year
#'   `t + h < test_start` are included (strict inequality prevents leakage).
#'
#' @return A data.table with all original columns plus a `.target` column holding
#'   the h-step-ahead value of the outcome.
#' @family long_horizon
#' @export
create_lh_data <- function(data, outcome, h, groupvar, timevar, test_start) {
  if (!data.table::is.data.table(data)) {
    data <- data.table::as.data.table(as.data.frame(data))
  }
  dt <- data.table::copy(data)
  data.table::setkeyv(dt, c(groupvar, timevar))

  # Add .target = lead(outcome, h) per group (positional, sorted by time)
  dt[, .target := lead_horizon(get(outcome), h), by = c(groupvar)]

  # Filter to valid training rows
  dt <- dt[dt[[timevar]] < (test_start - h) & !is.na(.target)]
  return(dt)
}


#' Fit a long-horizon OLS model
#'
#' Fits an ordinary least-squares model (with optional bootstrap) that predicts
#' the h-step-ahead outcome from baseline covariates at the forecast origin t.
#'
#' The formula LHS must be wrapped in [lead_horizon()] and names the outcome,
#' optionally transformed, e.g. `lead_horizon(gdppc_grwt) ~ ...` or
#' `lead_horizon(asinh(gdppc)) ~ ...`. Internally the `lead_horizon()` wrapper is
#' stripped and the outcome variable is replaced with `.target` (the led column
#' added by [create_lh_data()]) so that:
#' * LHS = transformation applied to the h-step-ahead value
#' * RHS = covariates at the baseline year t (write `lag()` only for older history)
#'
#' @param formula Formula. Model specification; LHS is `lead_horizon(outcome)`
#'   (optionally transformed), RHS uses baseline covariates at the origin.
#' @param aligned_data A data.table from [create_lh_data()].
#' @param h Integer. Forecast horizon (stored for bookkeeping).
#' @param groupvar Character. Group variable name.
#' @param timevar Character. Time variable name.
#' @param boot Character or NULL. Bootstrap type: `"resid"`, `"wild"`, or `NULL`
#'   for standard OLS.
#'
#' @return A list with model components for use in [forecast_long_horizon()]:
#'   `formula`, `aligned_formula`, `naive_formula`, `h`, `fitted`, `sigma`,
#'   `df`, `fit_data`, `boot`.
#' @family long_horizon
#' @export
fit_lh_model <- function(formula, aligned_data, h, groupvar, timevar, boot = NULL) {
  if (!data.table::is.data.table(aligned_data)) {
    aligned_data <- data.table::as.data.table(as.data.frame(aligned_data))
  }

  # Build aligned formula: strip the lead_horizon() marker from the LHS and
  # replace the outcome variable with .target (the led column from create_lh_data()).
  inner_lhs   <- .lh_lhs_inner(formula)
  outcome_var <- all.vars(inner_lhs)[1L]
  new_lhs     <- .sub_sym(inner_lhs, outcome_var, ".target")
  aligned_formula <- rlang::new_formula(new_lhs, rlang::f_rhs(formula),
                                        env = rlang::f_env(formula))

  # Extend formula with timevar for row tracking
  aligned_formula_ext <- stats::update(aligned_formula,
                                       paste(". ~ . +", timevar))

  # Inject positional lag for per-group evaluation
  aligned_formula_ext <- inject_positional_lag(aligned_formula_ext)

  # Evaluate formula terms per group using data.table
  model_data <- aligned_data[, {
    mf <- tryCatch(
      stats::model.frame(aligned_formula_ext, data = .SD, na.action = NULL),
      error = function(e) NULL
    )
    if (is.null(mf)) return(data.table::data.table())
    mm <- stats::model.matrix(stats::terms(mf), data = mf,
                              na.action = stats::na.pass)
    mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
    resp <- stats::setNames(data.frame(mf[[1L]], stringsAsFactors = FALSE),
                            names(mf)[1L])
    data.table::as.data.table(cbind(resp, as.data.frame(mm)))
  }, by = c(groupvar)]

  # Name the response `.target` (it sits right after the group key) and clean
  # only the remaining columns. Keeping the dotted, reserved-style name and
  # excluding it from clean_names() guarantees the response can never collide
  # with a user covariate, and that predictor names match those produced at
  # predict time (where the response is absent).
  resp_pos <- length(groupvar) + 1L
  data.table::setnames(model_data, resp_pos, ".target")
  other_cols <- setdiff(names(model_data), ".target")
  data.table::setnames(model_data, other_cols,
                       names(janitor::clean_names(model_data[, ..other_cols])))

  # Identify columns to use for fitting (drop groupvar and timevar)
  grp_clean  <- .clean_name(groupvar)
  time_clean <- .clean_name(timevar)
  all_vars   <- names(model_data)
  fit_vars   <- all_vars[!all_vars %in% c(grp_clean, time_clean)]

  fit_data <- stats::na.omit(model_data[, ..fit_vars])

  # Naive formula over cleaned column names (response is .target)
  naive_formula <- stats::as.formula(
    paste(fit_vars[1L], "~", paste(fit_vars[-1L], collapse = " + "))
  )

  if (!is.null(boot)) {
    fitted <- bootstraplm(naive_formula, fit_data, type = boot)
  } else {
    fitted <- stats::lm(naive_formula, fit_data)
  }

  list(
    formula         = formula,
    aligned_formula = aligned_formula,
    naive_formula   = naive_formula,
    h               = h,
    fitted          = fitted,
    sigma           = summary(fitted)$sigma,
    df              = fitted$df.residual,
    fit_data        = fit_data,
    boot            = boot
  )
}


#' Set up long-horizon models for all horizons and formula variants
#'
#' For each combination of horizon and named formula, calls [create_lh_data()]
#' and [fit_lh_model()] to produce a complete set of fitted models. No dynamic
#' loop is needed — models are independent across time and estimated once.
#'
#' Each formula must have its LHS wrapped in [lead_horizon()]; for example
#' `list(lh_linear = lead_horizon(gdppc_grwt) ~ gdppc_grwt + log(gdppc) + best)`.
#' The right-hand side is evaluated at the forecast origin (`test_start - 1`),
#' matching the dynamic simulator's information set. Use `lag()` on the RHS only
#' when you deliberately want history older than the origin.
#'
#' @param data A data.frame or data.table.
#' @param formulas A formula or a list of formulas, one per model variant. Each
#'   LHS must be wrapped in [lead_horizon()]. A bare formula is accepted as a
#'   single variant; unnamed (or partially named) list entries are auto-named
#'   `model1`, `model2`, ... and that name becomes the `variant` label.
#' @param horizons Integer vector of forecast horizons (e.g. `1:12`).
#' @param groupvar Character. Group variable name.
#' @param timevar Character. Time variable name.
#' @param test_start Numeric. First forecast year. Training uses only rows whose
#'   h-step target falls strictly before `test_start`, so `test_start` itself is
#'   the first out-of-sample target. Same meaning as `test_start` in
#'   [setup_simulator()].
#' @param boot Character or NULL. Bootstrap type: `"resid"`, `"wild"`, or `NULL`.
#'
#' @return A `lh_setup` list suitable for [forecast_long_horizon()].
#' @family long_horizon
#' @export
setup_long_horizon <- function(data, formulas, horizons, groupvar, timevar,
                               test_start, boot = NULL) {
  formulas <- .name_formulas(formulas)
  models <- list()

  for (variant in names(formulas)) {
    f           <- formulas[[variant]]
    outcome_var <- .lh_outcome(f)
    models[[variant]] <- list()

    for (h in horizons) {
      aligned <- create_lh_data(data, outcome_var, h, groupvar, timevar, test_start)
      models[[variant]][[as.character(h)]] <-
        fit_lh_model(f, aligned, h, groupvar, timevar, boot)
    }
  }

  list(
    models    = models,
    horizons  = horizons,
    formulas  = formulas,
    groupvar  = groupvar,
    timevar   = timevar,
    test_start = test_start,
    boot      = boot
  )
}


# Internal helper: draw predictive samples for one (lh_model, test_start) pair.
# Nested Monte Carlo: `nsim` bootstrap refits (parameter uncertainty) x
# `inner_sims` predictive draws per refit (residual uncertainty) = nsim *
# inner_sims draws per unit. With boot = NULL there is one fit and only the
# product matters. Returns a data.table with columns: groupvar, .draws.
.predict_lh <- function(lh_model, data, groupvar, timevar, test_start,
                        nsim, inner_sims) {
  if (!data.table::is.data.table(data)) {
    data <- data.table::as.data.table(as.data.frame(data))
  }

  formula <- lh_model$formula
  boot    <- lh_model$boot

  # Build one-sided RHS formula extended with timevar for row tracking
  rhs_str         <- deparse(rlang::f_rhs(formula), width.cutoff = 500L)
  rhs_formula_ext <- stats::as.formula(paste("~", rhs_str, "+", timevar))

  # Inject positional lag
  rhs_formula_ext <- inject_positional_lag(rhs_formula_ext)

  # Covariates are observed at the forecast origin (test_start - 1), the last
  # observed period; include enough history for any lag() on the RHS.
  origin <- test_start - 1L
  recent <- data[data[[timevar]] <= origin]
  data.table::setkeyv(recent, c(groupvar, timevar))
  recent <- recent[, utils::tail(.SD, 20L), by = c(groupvar)]

  # Evaluate RHS formula per group
  rhs_data <- recent[, {
    mf <- tryCatch(
      stats::model.frame(rhs_formula_ext, data = .SD, na.action = NULL),
      error = function(e) NULL
    )
    if (is.null(mf)) return(data.table::data.table())
    mm <- stats::model.matrix(stats::terms(mf), data = mf,
                              na.action = stats::na.pass)
    mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
    data.table::as.data.table(as.data.frame(mm))
  }, by = c(groupvar)]

  rhs_data <- janitor::clean_names(rhs_data)

  grp_clean  <- .clean_name(groupvar)
  time_clean <- .clean_name(timevar)

  rhs_at_t <- rhs_data[rhs_data[[time_clean]] == origin]
  # Remove rows with any NA
  complete <- stats::complete.cases(rhs_at_t)
  rhs_at_t <- rhs_at_t[complete]

  if (nrow(rhs_at_t) == 0L) {
    return(data.table::data.table(
      placeholder_group = character(0L),
      .draws = list()
    ))
  }

  groups    <- rhs_at_t[[grp_clean]]
  pred_cols <- setdiff(names(rhs_at_t), c(grp_clean, time_clean))
  pred_data <- rhs_at_t[, ..pred_cols]
  n_total   <- nsim * inner_sims

  if (is.null(boot)) {
    # Standard: single fit, draw n_total samples
    pred   <- predict(lh_model$fitted, newdata = pred_data, se.fit = TRUE)
    se_pi  <- sqrt(pred$se.fit^2 + lh_model$sigma^2)
    draws_mat <- pred$fit + outer(se_pi, stats::rt(n_total, df = lh_model$df))
  } else {
    # Bootstrap: refit nsim times, draw inner_sims per refit
    draws_mat <- matrix(NA_real_, nrow = nrow(pred_data), ncol = n_total)
    for (i in seq_len(nsim)) {
      boot_fit <- bootstraplm(lh_model$naive_formula, lh_model$fit_data,
                              type = boot)
      pred_i   <- predict(boot_fit, newdata = pred_data, se.fit = TRUE)
      sigma_i  <- summary(boot_fit)$sigma
      df_i     <- boot_fit$df.residual
      se_pi_i  <- sqrt(pred_i$se.fit^2 + sigma_i^2)
      col_idx  <- seq_len(inner_sims) + (i - 1L) * inner_sims
      draws_mat[, col_idx] <-
        pred_i$fit + outer(se_pi_i, stats::rt(inner_sims, df = df_i))
    }
  }

  draws_list <- lapply(seq_len(nrow(draws_mat)), function(i) draws_mat[i, ])

  result <- data.table::data.table(
    placeholder_group = groups,
    .draws = draws_list
  )
  data.table::setnames(result, "placeholder_group", groupvar)
  result
}


#' Generate probabilistic long-horizon forecasts
#'
#' For each (variant, horizon) combination in a `lh_setup` object, draws
#' `nsim * inner_sims` samples from the predictive distribution conditioned on
#' covariates observed at the forecast origin (`test_start - 1`).
#'
#' @details
#' The predictive sample is drawn with a nested Monte Carlo, identical in spirit
#' to the dynamic simulator. `nsim` is the number of bootstrap *refits*
#' (parameter / model uncertainty) and `inner_sims` the number of predictive
#' draws *per refit* (residual / sampling uncertainty), for a total of
#' `nsim * inner_sims` draws per unit and horizon. Refitting is the expensive
#' step, so taking several cheap predictive draws per refit (`inner_sims > 1`)
#' is an efficiency lever, not an extra source of randomness. With `boot = NULL`
#' (plain OLS) there is a single fit and only the product `nsim * inner_sims`
#' matters. These arguments deliberately mirror `nsim` in [simulate_endogenr()]
#' and `inner_sims` in [setup_simulator()]; pass the same values to both
#' approaches to keep a benchmark apples-to-apples.
#'
#' @param lh_setup A setup object from [setup_long_horizon()].
#' @param data The full dataset (original, not aligned). Used to extract
#'   baseline covariates at the origin.
#' @param test_start Numeric or NULL. The first forecast year (same convention as
#'   [simulate_endogenr()]). Covariates are observed at `test_start - 1`, horizon
#'   h=1 targets `test_start`, h=2 targets `test_start + 1`, etc. Defaults to the
#'   `test_start` stored in `lh_setup`; a value below it triggers a leakage
#'   warning.
#' @param nsim Integer. Number of bootstrap refits (parameter uncertainty) when
#'   `boot` is set; with plain OLS it only scales the total draw count. Mirrors
#'   `nsim` in [simulate_endogenr()].
#' @param inner_sims Integer. Number of predictive draws drawn per refit
#'   (residual uncertainty). Total draws per unit and horizon are
#'   `nsim * inner_sims`. Mirrors `inner_sims` in [setup_simulator()].
#'
#' @return A data.table with columns: `groupvar`, `test_start`, `horizon`,
#'   `year_target`, `variant`, `.draws` (list-column of numeric vectors). The
#'   result carries `lh_formulas`, `panel_unit`, and `panel_time` attributes so
#'   [get_lh_accuracy()] can recover the per-variant outcomes and panel context.
#' @family long_horizon
#' @export
forecast_long_horizon <- function(lh_setup, data, test_start = NULL,
                                  nsim = 500L, inner_sims = 10L) {
  groupvar <- lh_setup$groupvar
  timevar  <- lh_setup$timevar

  if (is.null(test_start)) test_start <- lh_setup$test_start
  if (!is.null(lh_setup$test_start) && test_start < lh_setup$test_start) {
    warning("test_start (", test_start, ") is before the setup's test_start (",
            lh_setup$test_start, "); forecasts may overlap the training ",
            "window (leakage).", call. = FALSE)
  }

  results <- list()

  for (variant in names(lh_setup$models)) {
    for (h_char in names(lh_setup$models[[variant]])) {
      lh_model <- lh_setup$models[[variant]][[h_char]]
      h        <- lh_model$h

      pred_dt <- .predict_lh(lh_model, data, groupvar, timevar, test_start,
                             nsim, inner_sims)

      pred_dt[, `:=`(test_start = test_start,
                     horizon = h,
                     year_target = test_start + h - 1L,
                     variant = variant)]
      results[[paste(variant, h_char, sep = "_")]] <- pred_dt
    }
  }

  out <- data.table::rbindlist(results)
  # Stamp metadata so get_lh_accuracy() can recover the per-variant formulas
  # (hence outcomes) and the panel context without re-passing them.
  data.table::setattr(out, "lh_formulas", lh_setup$formulas)
  data.table::setattr(out, "panel_unit", groupvar)
  data.table::setattr(out, "panel_time", timevar)
  out
}


#' Compute CRPS, MAE, and Winkler accuracy scores for long-horizon forecasts
#'
#' Scores forecast draws against observed truth per `(variant, unit, horizon)`
#' using the same NA-safe kernel as [get_accuracy()]. The output shares the
#' `(unit, horizon)` key with `get_accuracy(..., test_start = , by = )`, so the
#' two approaches can be stacked directly with [compare_approaches()].
#'
#' The outcome is taken *per variant* from each formula's LHS (via `lh_setup` or
#' the `lh_formulas` attribute stamped by [forecast_long_horizon()]), so a setup
#' mixing outcomes — e.g.
#' `list(m1 = lead_horizon(gdppc) ~ ., m2 = lead_horizon(population) ~ .)` —
#' is scored correctly, each variant against its own outcome. `truth` must
#' therefore contain every outcome column referenced by the formulas.
#'
#' Group and time columns come from `lh_setup` when supplied, otherwise inferred
#' from the `panel_unit` / `panel_time` attributes on `lh_forecasts` (then
#' `truth`). Pass `lh_setup` for robustness when those attributes may have been
#' dropped (e.g. after manually rebinding forecasts, as [cv_long_horizon()] does).
#'
#' Forecast draws are on the *modelled* scale (the formula LHS). `scale` controls
#' the scoring scale:
#' * `"native"` (default): score raw draws against the raw outcome. For a
#'   non-identity LHS (e.g. `asinh(gdppc)`) supply `inverse` (e.g. `sinh`) to
#'   back-transform the draws; with several transformed variants pass a named
#'   list of inverse functions keyed by variant.
#' * `"model"`: apply each formula's LHS transform to the truth and score on the
#'   modelled scale (no inverse needed).
#'
#' @param lh_forecasts Output from [forecast_long_horizon()] or [cv_long_horizon()].
#' @param truth A data.frame or data.table containing the outcome column(s)
#'   referenced by the formulas, plus the group and time columns.
#' @param lh_setup Optional setup object from [setup_long_horizon()]. Supplies
#'   the group/time columns and the per-variant formulas. When `NULL`, these are
#'   inferred from attributes on `lh_forecasts`.
#' @param scale One of `"native"` or `"model"`. See Details.
#' @param inverse Function, named list of functions (keyed by variant), or
#'   `NULL`. Inverse of the LHS transform applied to the draws under
#'   `scale = "native"` for variants with a non-identity LHS.
#' @param level Numeric. Coverage level for the Winkler score (default 50).
#'
#' @return A data.table with columns: `variant`, the unit column, `horizon`,
#'   `crps`, `mae`, `winkler`.
#' @family long_horizon
#' @export
get_lh_accuracy <- function(lh_forecasts, truth, lh_setup = NULL,
                            scale = c("native", "model"),
                            inverse = NULL, level = 50) {
  scale <- match.arg(scale)

  if (!data.table::is.data.table(lh_forecasts)) {
    lh_forecasts <- data.table::as.data.table(lh_forecasts)
  }
  if (!data.table::is.data.table(truth)) {
    truth <- data.table::as.data.table(as.data.frame(truth))
  }

  # Resolve group/time columns and the per-variant formulas.
  if (!is.null(lh_setup)) {
    groupvar <- lh_setup$groupvar
    timevar  <- lh_setup$timevar
    formulas <- lh_setup$formulas
  } else {
    formulas <- attr(lh_forecasts, "lh_formulas", exact = TRUE)
    ctx <- .infer_ctx(lh_forecasts, truth)
    if (is.null(formulas) || is.null(ctx)) {
      stop("Cannot recover formulas / panel context. Pass `lh_setup` from ",
           "setup_long_horizon(), or score a forecast object that still ",
           "carries its `lh_formulas` / `panel_unit` / `panel_time` attributes.",
           call. = FALSE)
    }
    groupvar <- ctx_unit(ctx)
    timevar  <- ctx_time(ctx)
  }
  formulas <- .name_formulas(formulas)

  variants <- unique(as.character(lh_forecasts$variant))

  # Per variant: inner LHS expression and the raw outcome variable.
  inner_by_v   <- lapply(stats::setNames(variants, variants),
                         function(v) .lh_lhs_inner(formulas[[v]]))
  outcome_by_v <- vapply(inner_by_v, function(e) all.vars(e)[1L], character(1))

  # Truth columns needed: all LHS vars (model) or the outcomes (native), + keys.
  needed <- if (scale == "model") {
    unique(unlist(lapply(inner_by_v, all.vars)))
  } else {
    unique(unname(outcome_by_v))
  }
  truth_cols <- unique(c(groupvar, timevar, needed))
  truth_sub  <- truth[, ..truth_cols]
  data.table::setnames(truth_sub, timevar, "year_target")

  scored <- merge(lh_forecasts, truth_sub, by = c(groupvar, "year_target"),
                  all.x = TRUE)

  # Per-variant truth on the scoring scale (and back-transformed draws for a
  # transformed LHS scored on the native scale).
  scored[, .truth_score := NA_real_]
  for (v in variants) {
    idx <- which(scored$variant == v)
    if (length(idx) == 0L) next
    inner_v   <- inner_by_v[[v]]
    outcome_v <- outcome_by_v[[v]]

    if (scale == "model") {
      data.table::set(scored, i = idx, j = ".truth_score",
                      value = as.numeric(eval(inner_v, scored[idx])))
    } else {
      data.table::set(scored, i = idx, j = ".truth_score",
                      value = as.numeric(scored[[outcome_v]][idx]))
      is_identity <- is.name(inner_v) &&
        identical(as.character(inner_v), outcome_v)
      if (!is_identity) {
        inv_v <- if (is.list(inverse)) inverse[[v]] else inverse
        if (is.null(inv_v)) {
          stop("scale = 'native' with a transformed LHS (", deparse(inner_v),
               ") for variant '", v, "' requires an `inverse` function.",
               call. = FALSE)
        }
        for (k in idx) scored$.draws[[k]] <- inv_v(scored$.draws[[k]])
      }
    }
  }

  sc <- .score_draws(scored$.draws, scored$.truth_score, level)
  scored[, `:=`(crps = sc$crps, mae = sc$mae, winkler = sc$winkler)]

  result <- scored[, .(crps = mean(crps, na.rm = TRUE),
                       mae = mean(mae, na.rm = TRUE),
                       winkler = mean(winkler, na.rm = TRUE)),
                   by = c("variant", groupvar, "horizon")]
  return(result)
}


#' Cross-validate long-horizon prediction models
#'
#' Runs the full pipeline — [setup_long_horizon()] + [forecast_long_horizon()] —
#' for each training window defined by `test_starts`, stacks all forecasts, and
#' returns aggregated accuracy scores per `(variant, unit, horizon)`.
#'
#' @param data A data.frame or data.table.
#' @param formulas A formula or list of formulas (see [setup_long_horizon()]);
#'   bare formulas and unnamed lists are accepted and auto-named.
#' @param horizons Integer vector of forecast horizons.
#' @param groupvar Character. Group variable name.
#' @param timevar Character. Time variable name.
#' @param test_starts Numeric vector of first forecast years. Each value is used
#'   as the `test_start` for its fold.
#' @param boot Character or NULL. Bootstrap type.
#' @param nsim Integer. Number of bootstrap refits per fold (see
#'   [forecast_long_horizon()]); total draws per cell are `nsim * inner_sims`.
#' @param inner_sims Integer. Number of predictive draws per refit (see
#'   [forecast_long_horizon()]).
#' @param scale One of `"native"` or `"model"`; see [get_lh_accuracy()].
#' @param inverse Function or NULL; see [get_lh_accuracy()].
#' @param level Numeric. Coverage level for the Winkler score (default 50).
#'
#' @return A data.table with columns: `variant`, the unit column, `horizon`,
#'   `crps`, `mae`, `winkler`, averaged across all folds.
#' @family long_horizon
#' @export
cv_long_horizon <- function(data, formulas, horizons, groupvar, timevar,
                            test_starts, boot = NULL, nsim = 500L,
                            inner_sims = 10L, scale = c("native", "model"),
                            inverse = NULL, level = 50) {
  scale <- match.arg(scale)
  formulas <- .name_formulas(formulas)
  all_forecasts <- list()

  for (ts in test_starts) {
    lh_setup  <- setup_long_horizon(data, formulas, horizons, groupvar, timevar,
                                    test_start = ts, boot = boot)
    forecasts <- forecast_long_horizon(lh_setup, data, test_start = ts,
                                       nsim = nsim, inner_sims = inner_sims)
    all_forecasts[[as.character(ts)]] <- forecasts
  }

  combined <- data.table::rbindlist(all_forecasts)

  if (!data.table::is.data.table(data)) {
    data <- data.table::as.data.table(as.data.frame(data))
  }

  # Truth needs every column referenced by a formula LHS (covers both scales).
  needed <- unique(unlist(lapply(formulas, function(f) all.vars(.lh_lhs_inner(f)))))
  truth_cols <- unique(c(groupvar, timevar, needed))
  truth <- data[, ..truth_cols]

  # All folds share group/time/formulas; the last fold's setup carries that
  # metadata (rbindlist drops the stamped attributes).
  get_lh_accuracy(combined, truth, lh_setup = lh_setup,
                  scale = scale, inverse = inverse, level = level)
}


#' Stack long-horizon and simulation accuracy results for comparison
#'
#' Combines the output of [get_lh_accuracy()] / [cv_long_horizon()] with
#' simulation accuracy from [get_accuracy()], adding an `approach` column. For a
#' like-for-like comparison both inputs should share the `(unit, horizon)` key —
#' score the simulation with `get_accuracy(..., test_start = , by = c(unit, "horizon"))`.
#'
#' @param lh_accuracy A data.table from [get_lh_accuracy()] or [cv_long_horizon()].
#' @param sim_accuracy A data.table from [get_accuracy()]. For per-horizon
#'   comparison it must carry a `horizon` column (and the unit column).
#'
#' @return A data.table with an `approach` column plus the shared score columns;
#'   `variant` is `NA` for simulation rows.
#' @family long_horizon
#' @export
compare_approaches <- function(lh_accuracy, sim_accuracy) {
  if (!data.table::is.data.table(lh_accuracy)) {
    lh_accuracy <- data.table::as.data.table(lh_accuracy)
  }
  if (!data.table::is.data.table(sim_accuracy)) {
    sim_accuracy <- data.table::as.data.table(sim_accuracy)
  }

  lh_tbl <- data.table::copy(lh_accuracy)
  lh_tbl[, approach := "long_horizon"]

  sim_tbl <- data.table::copy(sim_accuracy)
  if (!"variant" %in% names(sim_tbl)) sim_tbl[, variant := NA_character_]
  sim_tbl[, approach := "simulation"]

  data.table::rbindlist(list(lh_tbl, sim_tbl), fill = TRUE)
}
