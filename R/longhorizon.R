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


#' Create aligned (t, t+h) training data for long-horizon models
#'
#' For each group, computes the h-step-ahead value of the outcome and stores it
#' in a new `.target` column. The original outcome column is left untouched so
#' that formula RHS terms (e.g. `lag(outcome)`) are evaluated at the baseline
#' period t, not the target period t+h.
#'
#' @param data A data.frame or tsibble containing all variables.
#' @param outcome Character. The raw outcome variable name (without transformation).
#' @param h Integer. The forecast horizon.
#' @param groupvar Character. The group variable name.
#' @param timevar Character. The time variable name.
#' @param train_end Numeric. Training cutoff; only rows where `t + h < train_end`
#'   are included (strict inequality prevents outcome leakage).
#'
#' @return A tibble with all original columns plus a `.target` column holding
#'   the h-step-ahead value of the outcome.
#' @export
create_lh_data <- function(data, outcome, h, groupvar, timevar, train_end) {
  outcome_sym <- rlang::sym(outcome)
  group_sym   <- rlang::sym(groupvar)
  time_sym    <- rlang::sym(timevar)

  data |>
    dplyr::as_tibble() |>
    dplyr::group_by(!!group_sym) |>
    dplyr::arrange(!!time_sym, .by_group = TRUE) |>
    dplyr::mutate(.target = dplyr::lead(!!outcome_sym, n = h)) |>
    dplyr::ungroup() |>
    dplyr::filter(!!time_sym < (train_end - h), !is.na(.target))
}


#' Fit a long-horizon OLS model
#'
#' Fits an ordinary least-squares model (with optional bootstrap) that predicts
#' the h-step-ahead outcome from baseline covariates at time t.
#'
#' The formula LHS specifies the outcome variable and its transformation, e.g.
#' `asinh(gdppc) ~ lag(asinh(gdppc))`. Internally the LHS variable is replaced
#' with `.target` (the column added by [create_lh_data()]) so that:
#' * LHS = transformation applied to the h-step-ahead value
#' * RHS = covariates at the baseline year t
#'
#' @param formula Formula. Model specification; LHS is the transformed outcome,
#'   RHS uses baseline covariates (may include `lag()`, raw variables, etc.).
#' @param aligned_data A tibble from [create_lh_data()].
#' @param h Integer. Forecast horizon (stored for bookkeeping).
#' @param groupvar Character. Group variable name.
#' @param timevar Character. Time variable name.
#' @param boot Character or NULL. Bootstrap type: `"resid"`, `"wild"`, or `NULL`
#'   for standard OLS.
#'
#' @return A list with model components for use in [forecast_long_horizon()]:
#'   `formula`, `aligned_formula`, `naive_formula`, `h`, `fitted`, `sigma`,
#'   `df`, `fit_data`, `boot`.
#' @export
fit_lh_model <- function(formula, aligned_data, h, groupvar, timevar, boot = NULL) {
  # Build aligned formula: replace outcome var in LHS with .target
  lhs_expr    <- rlang::f_lhs(formula)
  outcome_var <- all.vars(lhs_expr)[1L]
  new_lhs     <- .sub_sym(lhs_expr, outcome_var, ".target")
  aligned_formula <- rlang::new_formula(new_lhs, rlang::f_rhs(formula),
                                        env = rlang::f_env(formula))

  # Extend formula with timevar for row tracking inside group_modify
  aligned_formula_ext <- stats::update(aligned_formula,
                                       paste(". ~ . +", timevar))

  # Evaluate formula terms per group (handles lag() correctly)
  model_data <- aligned_data |>
    dplyr::group_by(!!rlang::sym(groupvar)) |>
    dplyr::group_modify(~ {
      mf <- tryCatch(
        stats::model.frame(aligned_formula_ext, data = ., na.action = NULL),
        error = function(e) NULL
      )
      if (is.null(mf)) return(tibble::tibble())
      tibble::as_tibble(mf)
    }) |>
    dplyr::ungroup() |>
    janitor::clean_names()

  # Identify columns to use for fitting (drop groupvar and timevar)
  grp_clean  <- .clean_name(groupvar)
  time_clean <- .clean_name(timevar)
  all_vars   <- names(model_data)
  fit_vars   <- all_vars[!all_vars %in% c(grp_clean, time_clean)]

  fit_data <- stats::na.omit(model_data[, fit_vars, drop = FALSE])

  # Naive formula over cleaned column names
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
#' @param data A data.frame or tsibble.
#' @param formulas Named list of formulas, one per model variant. Example:
#'   ```r
#'   list(
#'     baseline = asinh(gdppc) ~ lag(asinh(gdppc)),
#'     full     = asinh(gdppc) ~ lag(asinh(gdppc)) + lag(yjbest) + lag(dem)
#'   )
#'   ```
#' @param horizons Integer vector of forecast horizons (e.g. `1:12`).
#' @param groupvar Character. Group variable name.
#' @param timevar Character. Time variable name.
#' @param train_end Numeric. Training cutoff; only observations before this year
#'   (shifted by h) are used for training.
#' @param boot Character or NULL. Bootstrap type: `"resid"`, `"wild"`, or `NULL`.
#'
#' @return A `lh_setup` list suitable for [forecast_long_horizon()].
#' @export
setup_long_horizon <- function(data, formulas, horizons, groupvar, timevar,
                               train_end, boot = NULL) {
  models <- list()

  for (variant in names(formulas)) {
    f           <- formulas[[variant]]
    outcome_var <- all.vars(rlang::f_lhs(f))[1L]
    models[[variant]] <- list()

    for (h in horizons) {
      aligned <- create_lh_data(data, outcome_var, h, groupvar, timevar, train_end)
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
    train_end = train_end,
    boot      = boot
  )
}


# Internal helper: draw predictive samples for one (lh_model, test_start) pair.
# Returns a tibble with columns: groupvar, .draws (list of numeric vectors).
.predict_lh <- function(lh_model, data, groupvar, timevar, test_start,
                        nsim, inner_sims) {
  formula <- lh_model$formula
  boot    <- lh_model$boot

  # Build one-sided RHS formula extended with timevar for row tracking
  rhs_str         <- deparse(rlang::f_rhs(formula), width.cutoff = 500L)
  rhs_formula_ext <- stats::as.formula(paste("~", rhs_str, "+", timevar))

  # Include enough history for lag() computation
  recent <- data |>
    dplyr::as_tibble() |>
    dplyr::group_by(!!rlang::sym(groupvar)) |>
    dplyr::filter(!!rlang::sym(timevar) <= test_start) |>
    dplyr::arrange(!!rlang::sym(timevar)) |>
    dplyr::slice_tail(n = 20L) |>
    dplyr::ungroup()

  # Evaluate RHS formula per group
  rhs_data <- recent |>
    dplyr::group_by(!!rlang::sym(groupvar)) |>
    dplyr::group_modify(~ {
      mf <- tryCatch(
        stats::model.frame(rhs_formula_ext, data = ., na.action = NULL),
        error = function(e) NULL
      )
      if (is.null(mf)) return(tibble::tibble())
      tibble::as_tibble(mf)
    }) |>
    dplyr::ungroup() |>
    janitor::clean_names()

  grp_clean  <- .clean_name(groupvar)
  time_clean <- .clean_name(timevar)

  rhs_at_t <- rhs_data |>
    dplyr::filter(!!rlang::sym(time_clean) == test_start) |>
    dplyr::filter(dplyr::if_all(dplyr::everything(), ~ !is.na(.)))

  if (nrow(rhs_at_t) == 0L) {
    return(tibble::tibble(!!rlang::sym(groupvar) := character(0L),
                          .draws = list()))
  }

  groups    <- rhs_at_t[[grp_clean]]
  pred_data <- rhs_at_t[, !names(rhs_at_t) %in% c(grp_clean, time_clean),
                        drop = FALSE]
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

  tibble::tibble(
    !!rlang::sym(groupvar) := groups,
    .draws = draws_list
  )
}


#' Generate probabilistic long-horizon forecasts
#'
#' For each (variant, horizon) combination in a `lh_setup` object, draws
#' `nsim * inner_sims` samples from the predictive distribution conditioned on
#' covariates observed at `test_start`.
#'
#' @param lh_setup A setup object from [setup_long_horizon()].
#' @param data The full dataset (original, not aligned). Used to extract
#'   baseline covariates at `test_start`.
#' @param test_start Numeric. The forecast origin (baseline year).
#' @param nsim Integer. Number of outer simulations (bootstrap draws if boot is
#'   set; otherwise multiplied with `inner_sims` for total draws).
#' @param inner_sims Integer. Number of inner draws per simulation.
#'
#' @return A tibble with columns: `groupvar`, `test_start`, `horizon`,
#'   `year_target`, `variant`, `.draws` (list-column of numeric vectors).
#' @export
forecast_long_horizon <- function(lh_setup, data, test_start,
                                  nsim = 500L, inner_sims = 10L) {
  groupvar <- lh_setup$groupvar
  timevar  <- lh_setup$timevar

  results <- list()

  for (variant in names(lh_setup$models)) {
    for (h_char in names(lh_setup$models[[variant]])) {
      lh_model <- lh_setup$models[[variant]][[h_char]]
      h        <- lh_model$h

      pred_df <- .predict_lh(lh_model, data, groupvar, timevar, test_start,
                             nsim, inner_sims)

      results[[paste(variant, h_char, sep = "_")]] <- pred_df |>
        dplyr::mutate(
          test_start  = test_start,
          horizon     = h,
          year_target = test_start + h,
          variant     = variant
        )
    }
  }

  dplyr::bind_rows(results)
}


#' Compute CRPS and MAE accuracy scores for long-horizon forecasts
#'
#' Wraps draw samples in [distributional::dist_sample()], creates a fable, and
#' uses [fabletools::accuracy()] — mirroring the pattern in [get_accuracy()].
#' Output format is compatible with [get_accuracy()] results (after aggregating
#' by horizon) for cross-approach comparison via [compare_approaches()].
#'
#' @param lh_forecasts Output from [forecast_long_horizon()] or [cv_long_horizon()].
#' @param truth A tsibble with `groupvar` as key, `timevar` as index, and the
#'   outcome column in the same scale as the forecast draws.
#' @param outcome Character. Name of the outcome column in `truth`; must match
#'   the scale of the forecast draws (e.g. if draws are `asinh(gdppc)`, pass a
#'   truth tsibble that also has an `asinh(gdppc)` column under this name).
#' @param groupvar Character. Group variable name.
#' @param timevar Character. Time variable name.
#'
#' @return A tibble with columns: `variant`, `horizon`, `crps`, `mae`.
#' @export
get_lh_accuracy <- function(lh_forecasts, truth, outcome, groupvar, timevar) {
  results <- list()

  for (v in unique(lh_forecasts$variant)) {
    for (h in unique(lh_forecasts$horizon)) {
      sub <- lh_forecasts |>
        dplyr::filter(variant == v, horizon == h)

      dist_tsbl <- sub |>
        dplyr::select(dplyr::all_of(c(groupvar, "year_target", ".draws"))) |>
        dplyr::mutate(
          !!rlang::sym(outcome) := distributional::dist_sample(.draws),
          !!rlang::sym(timevar) := year_target
        ) |>
        dplyr::select(-`.draws`, -year_target) |>
        tsibble::as_tsibble(key = groupvar, index = timevar)

      fbl <- fabletools::fable(dist_tsbl,
                               key    = groupvar,
                               index  = timevar,
                               response     = outcome,
                               distribution = outcome)

      metrics <- list(crps = fabletools::CRPS, mae = fabletools::MAE)
      acc     <- fabletools::accuracy(fbl, truth, metrics)

      results[[paste(v, h, sep = "_")]] <- dplyr::mutate(acc,
                                                          variant = v,
                                                          horizon = h)
    }
  }

  dplyr::bind_rows(results) |>
    dplyr::group_by(variant, horizon) |>
    dplyr::summarise(
      crps = mean(crps, na.rm = TRUE),
      mae  = mean(mae,  na.rm = TRUE),
      .groups = "drop"
    )
}


#' Cross-validate long-horizon prediction models
#'
#' Runs the full pipeline — [setup_long_horizon()] + [forecast_long_horizon()] —
#' for each training window defined by `test_starts`, stacks all forecasts, and
#' returns aggregated accuracy scores.
#'
#' @param data A data.frame or tsibble.
#' @param formulas Named list of formulas (see [setup_long_horizon()]).
#' @param horizons Integer vector of forecast horizons.
#' @param groupvar Character. Group variable name.
#' @param timevar Character. Time variable name.
#' @param test_starts Numeric vector of forecast origins, e.g. `c(2000, 2005, 2010)`.
#'   Each value is also used as `train_end` for that fold.
#' @param boot Character or NULL. Bootstrap type.
#' @param nsim Integer. Number of outer simulations.
#' @param inner_sims Integer. Number of inner draws per simulation.
#'
#' @return A tibble with columns: `variant`, `horizon`, `crps`, `mae`, averaged
#'   across all folds, units, and horizons.
#' @export
cv_long_horizon <- function(data, formulas, horizons, groupvar, timevar,
                            test_starts, boot = NULL, nsim = 500L,
                            inner_sims = 10L) {
  all_forecasts <- list()

  for (ts in test_starts) {
    lh_setup  <- setup_long_horizon(data, formulas, horizons, groupvar, timevar,
                                    train_end = ts, boot = boot)
    forecasts <- forecast_long_horizon(lh_setup, data, test_start = ts,
                                       nsim = nsim, inner_sims = inner_sims)
    all_forecasts[[as.character(ts)]] <- forecasts
  }

  combined    <- dplyr::bind_rows(all_forecasts)
  outcome_var <- all.vars(rlang::f_lhs(formulas[[1L]]))[1L]

  truth <- data |>
    dplyr::as_tibble() |>
    dplyr::select(dplyr::all_of(c(groupvar, timevar, outcome_var))) |>
    tsibble::as_tsibble(key = groupvar, index = timevar)

  get_lh_accuracy(combined, truth, outcome_var, groupvar, timevar)
}


#' Stack long-horizon and simulation accuracy results for comparison
#'
#' Combines the output of [get_lh_accuracy()] / [cv_long_horizon()] with
#' simulation accuracy aggregated by horizon, adding an `approach` column.
#'
#' @param lh_accuracy A tibble from [get_lh_accuracy()] or [cv_long_horizon()].
#' @param sim_accuracy A tibble from [get_accuracy()] after summarising by horizon.
#'   Must contain at least `horizon`, `crps`, and `mae` columns.
#'
#' @return A tidy tibble with columns: `approach`, `variant` (NA for simulation),
#'   `horizon`, `crps`, `mae`.
#' @export
compare_approaches <- function(lh_accuracy, sim_accuracy) {
  lh_tbl <- lh_accuracy |>
    dplyr::mutate(approach = "long_horizon")

  sim_tbl <- sim_accuracy |>
    dplyr::mutate(approach = "simulation",
                  variant  = NA_character_)

  dplyr::bind_rows(lh_tbl, sim_tbl)
}
