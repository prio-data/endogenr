
#' Prepares simulation data (subset data, expand simulations, include future horizon)
#'
#' @param data A data.table with training data.
#' @param ctx A panel_context object.
#' @param train_start Integer.
#' @param test_start Integer.
#' @param horizon Integer.
#' @param inner_sims Integer.
#'
#' @return A data.table with CJ grid of (unit, time, sim) filled with training data.
#' @keywords internal
prepare_simulation_data <- function(data, ctx, train_start, test_start, horizon, inner_sims) {
  unit_col <- ctx_unit(ctx)
  time_col <- ctx_time(ctx)

  # Training data only (data is already filtered to train_start..test_start+horizon-1)
  train <- data[data[[time_col]] < test_start]

  # Create CJ grid
  all_units <- unique(train[[unit_col]])
  all_times <- seq(train_start, test_start + horizon - 1)

  # Memory warning for large grids
  grid_rows <- length(all_units) * length(all_times) * inner_sims
  ncols <- ncol(train)
  est_mb <- grid_rows * ncols * 8 / 1e6  # rough estimate: 8 bytes per cell
  if (est_mb > 1000) {
    warning("Simulation grid will have ~", format(grid_rows, big.mark = ","),
            " rows (estimated ", round(est_mb), " MB). ",
            "Consider reducing inner_sims, horizon, or number of units.",
            call. = FALSE)
  }

  grid <- data.table::CJ(
    unit_placeholder = all_units,
    time_placeholder = all_times,
    sim = seq_len(inner_sims),
    sorted = FALSE
  )
  data.table::setnames(grid, c("unit_placeholder", "time_placeholder"), c(unit_col, time_col))

  # Merge training data onto grid
  result <- merge(grid, train, by = c(unit_col, time_col), all.x = TRUE, allow.cartesian = TRUE)

  # Sort by unit, sim, time for positional lag correctness
  data.table::setkeyv(result, c(unit_col, "sim", time_col))

  return(result)
}

#' Imputes the time-independent forecasts
#'
#' @param simulation_data A data.table.
#' @param models List of fitted models.
#' @param ctx A panel_context object.
#' @param test_start Integer.
#' @param horizon Integer.
#' @param inner_sims Integer.
#'
#' @return The simulation_data data.table, updated by reference.
#' @keywords internal
process_independent_models <- function(simulation_data, models, ctx, test_start, horizon, inner_sims) {
  if (!is.numeric(test_start)) {
    stop("`test_start` must be numeric")
  }

  join_keys <- c(ctx_keys(ctx), ctx_time(ctx))

  independent_models <- models[vapply(models, function(x) x$independent, logical(1))]

  for (model in independent_models) {
    outcome <- tryCatch(
      parse_formula(model)$outcome,
      error = function(e) "unknown"
    )

    pred <- tryCatch(
      predict(model, data = simulation_data, ctx = ctx,
              test_start = test_start, horizon = horizon, inner_sims = inner_sims),
      error = function(e) {
        stop("Prediction failed for independent model '", outcome, "': ",
             conditionMessage(e), call. = FALSE)
      }
    )

    cols <- setdiff(names(pred), join_keys)
    simulation_data[pred, (cols) := mget(paste0("i.", cols)), on = join_keys]
  }

  return(simulation_data)
}

#' Dynamic simulation of the time-dependent models
#'
#' @param simulation_data A data.table.
#' @param models List of fitted models.
#' @param ctx A panel_context object.
#' @param test_start Integer.
#' @param horizon Integer.
#' @param execution_order Character vector of model outcomes in topological order.
#'
#' @return The simulation_data data.table, updated by reference.
#' @keywords internal
process_dependent_models <- function(simulation_data, models, ctx, test_start, horizon, execution_order) {
  dependent_models <- models[vapply(models, function(x) !x$independent, logical(1))]
  outcomes <- vapply(dependent_models, function(x) parse_formula(x)$outcome, character(1))
  names(dependent_models) <- outcomes

  execution_order <- execution_order[execution_order %in% outcomes]
  dependent_models <- dependent_models[execution_order]

  join_keys <- c(ctx_keys(ctx), ctx_time(ctx))

  for (t in test_start:(test_start + horizon - 1)) {
    for (model_name in names(dependent_models)) {
      model <- dependent_models[[model_name]]
      pred <- tryCatch(
        predict(model, t = t, data = simulation_data, ctx = ctx),
        error = function(e) {
          stop("Prediction failed for model '", model_name,
               "' at time step ", t, ": ",
               conditionMessage(e), call. = FALSE)
        }
      )
      cols <- setdiff(names(pred), join_keys)
      simulation_data[pred, (cols) := mget(paste0("i.", cols)), on = join_keys]
    }
  }

  return(simulation_data)
}


#' Set up an endogenr simulator
#'
#' Validates the input panel, validates that the model system is closed, fits
#' all model specs once on the training window, builds the dependency graph,
#' and prepares the simulation grid. The returned setup object is the input
#' to [simulate_endogenr()].
#'
#' @details
#' During setup the following happens, in order:
#' \enumerate{
#'   \item The panel is coerced to a `data.table`, sorted by
#'     `(groupvar, timevar)`, and filtered to
#'     `[train_start, test_start + horizon - 1]`.
#'   \item [validate_panel()] checks contiguous integer time, balanced units,
#'     correct sort order, and complete initial state at `test_start - 1`.
#'   \item [validate_system_closure()] checks that every RHS variable is
#'     produced either by another model or by a data column, that all
#'     formula-referenced columns exist in `data`, and that no two models
#'     share an outcome.
#'   \item Each spec is fit via [fit_model()] once (or with a randomly drawn
#'     window from [get_train_window()] for `linear`/`glm`/`heterolm` when
#'     `min_window` is set).
#'   \item The dependency graph is built and a topological execution order is
#'     computed via [get_execution_order()].
#'   \item The cross-join simulation grid is created via
#'     `prepare_simulation_data()`.
#' }
#'
#' The returned list carries two contexts: `ctx` (with `sim = "sim"`, used at
#' predict time) and `fit_ctx` (without `sim`, used by inner refits). Both
#' are produced by [panel_context()].
#'
#' `min_window` triggers a random training window per outer simulation for
#' the parametric regression families. Leave it `NULL` for a fixed window
#' `[train_start, test_start - 1]`.
#'
#' `globals` is forwarded as-is to `future.apply::future_lapply()` so that
#' user-defined functions referenced by formulas (e.g. closure helpers) are
#' visible on parallel workers.
#'
#' @param models A list of model specs created by [build_model()].
#' @param data A data.frame, data.table, or tsibble. Must contain all
#'   variables referenced by model formulas, plus `groupvar`, `timevar`, and
#'   future exogenous data through `test_start + horizon - 1`.
#' @param train_start Integer. Earliest training time to include.
#' @param test_start Integer. First time step in the forecast period.
#' @param horizon Integer. Number of forecast time steps.
#' @param groupvar Character. Column name identifying panel units.
#' @param timevar Character. Column name identifying time periods.
#' @param inner_sims Integer. Number of inner simulations per outer sim.
#' @param min_window Integer or `NULL`. If set, `linear`/`glm`/`heterolm`
#'   models are refit on a random training window of at least this length
#'   per outer simulation.
#' @param globals Named list of user functions to export to parallel workers.
#'
#' @return A list of setup parameters for [simulate_endogenr()].
#' @seealso [build_model()], [simulate_endogenr()], [validate_panel()],
#'   [validate_system_closure()]
#' @family simulation
#' @export
#'
#' @examples
#' \dontrun{
#' df <- endogenr::example_data
#' system <- list(
#'   build_model("deterministic",
#'               formula = gdppc ~ I(abs(lag(gdppc) * (1 + gdppc_grwt)))),
#'   build_model("parametric_distribution",
#'               formula = ~gdppc_grwt, distribution = "norm"),
#'   build_model("exogen", formula = ~population)
#' )
#' setup <- setup_simulator(system, df,
#'   train_start = 1970, test_start = 2010, horizon = 5,
#'   groupvar = "gwcode", timevar = "year", inner_sims = 2)
#' }
setup_simulator <- function(models, data, train_start, test_start, horizon, groupvar, timevar,
                            inner_sims, min_window = NULL, globals = NULL) {
  # Create panel context
  ctx <- panel_context(unit = groupvar, time = timevar)

  # Convert to data.table
  if (!data.table::is.data.table(data)) {
    data <- data.table::as.data.table(as.data.frame(data))
  }

  # Filter to relevant time range
  data <- data[data[[timevar]] >= train_start & data[[timevar]] <= (test_start + horizon - 1)]
  data.table::setkeyv(data, c(groupvar, timevar))
  train <- data[data[[timevar]] < test_start]

  # Extract model outcomes for validation
  specs <- models
  independent_types <- get_independent_models()
  model_outcomes <- vapply(specs, function(spec) {
    f <- spec$formula
    if (spec$type %in% independent_types) {
      all.vars(f)
    } else {
      all.vars(rlang::f_lhs(f))
    }
  }, character(1))

  # Validate panel integrity
  validate_panel(data, ctx, test_start, model_outcomes)

  # Validate system closure (before fitting)
  validate_system_closure(specs, names(data))

  # Fit all models from specs
  fit_one <- function(spec) {
    type <- spec$type
    subset <- if (!is.null(min_window) && type %in% c("linear", "glm", "heterolm"))
      get_train_window(train_start, test_start, min_window) else NULL

    switch(type,
      "deterministic" = fit_model(spec, ctx = ctx),
      "linear" =, "glm" =, "heterolm" =
        fit_model(spec, data = train, ctx = ctx, subset = subset),
      "exogen" = fit_model(spec, newdata = data, ctx = ctx,
                           test_start = test_start, inner_sims = inner_sims),
      "parametric_distribution" =, "univariate_fable" =
        fit_model(spec, data = train, ctx = ctx),
      "spatial_lag" = fit_model(spec, ctx = ctx),
      stop("Unknown spec type: ", type)
    )
  }
  fitted_models <- lapply(specs, fit_one)

  # Build dependency graph
  dependency_graph <- igraph::make_empty_graph(directed = TRUE)
  for (model in fitted_models) {
    dependency_graph <- update_dependency_graph(model, dependency_graph)
  }
  execution_order <- get_execution_order(dependency_graph)

  # Create a simulation context WITH sim dimension
  sim_ctx <- panel_context(unit = groupvar, time = timevar, sim = "sim")

  simulation_data <- prepare_simulation_data(data, sim_ctx, train_start, test_start, horizon, inner_sims)

  return(list(
    "simulation_data" = simulation_data,
    "specs" = specs,
    "fitted_models" = fitted_models,
    "train_data" = train,
    "full_data" = data,
    "train_start" = train_start,
    "test_start" = test_start,
    "horizon" = horizon,
    "execution_order" = execution_order,
    "ctx" = sim_ctx,
    "fit_ctx" = ctx,
    "groupvar" = groupvar,
    "timevar" = timevar,
    "inner_sims" = inner_sims,
    "min_window" = min_window,
    "globals" = globals
  ))
}


#' The inner simulation loop
#'
#' @param i Simulation index.
#' @param specs List of endogenr_spec objects.
#' @param pre_fitted List of pre-fitted model objects.
#' @param train_data A data.table with training data.
#' @param simulation_data A data.table.
#' @param ctx A panel_context object (with sim dimension, for predict).
#' @param fit_ctx A panel_context object (without sim, for fitting).
#' @param test_start Integer.
#' @param horizon Integer.
#' @param execution_order Character vector.
#' @param inner_sims Integer.
#' @param min_window Integer or NULL.
#' @param train_start Integer.
#'
#' @return A data.table with simulation results.
#' @keywords internal
inner_simulation <- function(i, specs, pre_fitted, train_data, simulation_data,
                             ctx, fit_ctx, test_start, horizon, execution_order,
                             inner_sims, min_window, train_start) {
  refit_types <- c("linear", "glm", "heterolm")

  fitted_models <- lapply(seq_along(specs), function(j) {
    spec <- specs[[j]]
    if (spec$type %in% refit_types && !is.null(min_window)) {
      subset <- get_train_window(train_start, test_start, min_window)
      fit_model(spec, data = train_data, ctx = fit_ctx, subset = subset)
    } else {
      pre_fitted[[j]]
    }
  })

  sim <- data.table::copy(simulation_data)
  sim <- process_independent_models(sim, fitted_models, ctx, test_start, horizon, inner_sims)
  sim <- process_dependent_models(sim, fitted_models, ctx, test_start, horizon, execution_order)
  return(sim)
}

#' Run the dynamic simulation of an endogenr system
#'
#' Estimates models, sequences the forecast according to the dependency
#' graph computed in [setup_simulator()], and simulates outcomes for all
#' units across the forecast horizon. Returns one long `data.table` with a
#' `.sim` column identifying each draw.
#'
#' @details
#' Parallel execution is controlled by the user via the `future` package.
#' Set a plan before calling `simulate_endogenr()`:
#'
#' ```r
#' future::plan(future::multisession, workers = 6)
#' progressr::with_progress({
#'   res <- simulate_endogenr(nsim = 50, simulator_setup = setup)
#' })
#' future::plan(future::sequential)
#' ```
#'
#' Progress reporting uses [progressr::progressor()] when the `progressr`
#' package is installed; wrap the call in [progressr::with_progress()] to
#' see a bar.
#'
#' The output carries `panel_unit` / `panel_time` attributes (the
#' `paneltools::as_panel()` convention) so downstream helpers like
#' [get_accuracy()], [plotsim()], and [sim_to_dist()] infer the panel
#' context automatically.
#'
#' `parallel` and `ncores` are deprecated and retained only for backward
#' compatibility; they emit a deprecation message and, when `parallel =
#' TRUE`, register a temporary `future::multisession` plan that is restored
#' on exit. Prefer setting the plan yourself.
#'
#' @param nsim Integer. Number of outer-loop simulations.
#' @param simulator_setup A list from [setup_simulator()].
#' @param parallel Logical. Deprecated. Use [future::plan()] before calling.
#' @param ncores Integer. Deprecated. Use [future::plan()] before calling.
#'
#' @return A `data.table` with all simulation columns plus `.sim` (integer
#'   `1..nsim * inner_sims`). Carries `panel_unit` and `panel_time`
#'   attributes.
#' @seealso [setup_simulator()], [get_accuracy()], [plotsim()],
#'   [sim_to_dist()]
#' @family simulation
#' @export
simulate_endogenr <- function(nsim, simulator_setup, parallel = FALSE, ncores = 6) {
  if (!missing(parallel) || !missing(ncores)) {
    .Deprecated(
      msg = paste0(
        "The `parallel` and `ncores` arguments are deprecated.\n",
        "Set your own plan before calling simulate_endogenr(), e.g.:\n",
        "  future::plan(future::multisession, workers = 6)\n",
        "  result <- simulate_endogenr(nsim, setup)\n",
        "  future::plan(future::sequential)"
      )
    )
    if (parallel) {
      old_plan <- future::plan()
      on.exit(future::plan(old_plan), add = TRUE)
      future::plan(future::multisession, workers = ncores, gc = TRUE)
    }
  }

  # Build globals for future workers
  future_globals <- list(simulator_setup = simulator_setup)
  if (!is.null(simulator_setup$globals)) {
    for (fn_name in names(simulator_setup$globals)) {
      future_globals[[fn_name]] <- simulator_setup$globals[[fn_name]]
    }
  }

  # progressr support: users opt in via progressr::with_progress()
  p <- if (requireNamespace("progressr", quietly = TRUE)) {
    progressr::progressor(steps = nsim)
  } else {
    function(...) invisible(NULL)
  }

  simulation_results <- future.apply::future_lapply(
    seq_len(nsim),
    function(i) {
      result <- inner_simulation(
        i,
        specs = simulator_setup$specs,
        pre_fitted = simulator_setup$fitted_models,
        train_data = simulator_setup$train_data,
        simulation_data = simulator_setup$simulation_data,
        ctx = simulator_setup$ctx,
        fit_ctx = simulator_setup$fit_ctx,
        test_start = simulator_setup$test_start,
        horizon = simulator_setup$horizon,
        execution_order = simulator_setup$execution_order,
        inner_sims = simulator_setup$inner_sims,
        min_window = simulator_setup$min_window,
        train_start = simulator_setup$train_start
      )
      p()
      result
    },
    future.seed = TRUE,
    future.globals = future_globals,
    future.packages = c("endogenr", "data.table"),
    future.chunk.size = max(1L, ceiling(nsim / future::nbrOfWorkers()))
  )

  # Bind results and create .sim ID
  results <- data.table::rbindlist(simulation_results, idcol = ".id")

  # Create unique .sim across (outer_id, inner_sim) pairs
  ids <- unique(results[, .(.id, sim)])
  ids[, .sim := .I]
  results <- merge(results, ids, by = c(".id", "sim"))
  results[, c(".id", "sim") := NULL]

  # Stamp panel attributes (paneltools::as_panel convention) so downstream
  # helpers can infer the panel context without the user passing it.
  data.table::setattr(results, "panel_unit", ctx_unit(simulator_setup$ctx))
  data.table::setattr(results, "panel_time", ctx_time(simulator_setup$ctx))
  results[]

  return(results)
}

#' Nest simulation draws into list-columns
#'
#' Groups simulation results by (unit, time) and collapses each output variable
#' into a list of draws per cell.
#'
#' @param simulation_results A data.table from [simulate_endogenr()].
#' @param outputs Character vector of column names to nest.
#' @param ctx A `panel_context` object. Optional; inferred from
#'   `simulation_results` if it carries `panel_unit` / `panel_time` attributes
#'   (e.g. produced by [simulate_endogenr()] or `paneltools::as_panel()`).
#' @param sim_var Character. Name of the simulation ID column (default ".sim").
#'
#' @return A data.table with one row per (unit, time), output columns as list-columns.
#' @family postprocess
#' @export
sim_to_dist <- function(simulation_results, outputs, ctx = NULL, sim_var = ".sim") {
  if (is.null(ctx)) {
    ctx <- .infer_ctx(simulation_results)
    if (is.null(ctx)) {
      stop(
        "Could not infer panel context. Pass `ctx = panel_context(unit, time)` ",
        "or supply a `paneltools::as_panel()` data object.",
        call. = FALSE
      )
    }
  }
  keys <- c(ctx_unit(ctx), ctx_time(ctx))
  dt <- data.table::as.data.table(simulation_results)
  dt[, lapply(.SD, list), by = keys, .SDcols = outputs]
}

# Internal NA-safe scoring kernel shared by get_accuracy() and
# get_lh_accuracy(). Scores a list of draw vectors against a truth vector,
# returning CRPS, MAE, and Winkler scores. NA truth or an empty (all-NA) draw
# vector yields NA for that row.
.score_draws <- function(draws, y, level = 50) {
  alpha <- 1 - level / 100
  n <- length(draws)
  crps <- rep(NA_real_, n)
  mae  <- rep(NA_real_, n)
  winkler <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    d <- draws[[i]]
    d <- d[!is.na(d)]
    yi <- y[i]
    if (length(d) == 0L || is.na(yi)) next
    crps[i] <- scoringRules::crps_sample(y = yi, dat = d)
    mae[i]  <- abs(stats::median(d) - yi)
    lo <- unname(stats::quantile(d, alpha / 2))
    hi <- unname(stats::quantile(d, 1 - alpha / 2))
    penalty <- ifelse(yi < lo, (2 / alpha) * (lo - yi),
                ifelse(yi > hi, (2 / alpha) * (yi - hi), 0))
    winkler[i] <- (hi - lo) + penalty
  }
  list(crps = crps, mae = mae, winkler = winkler)
}

#' Calculate probabilistic accuracy scores for simulations
#'
#' Computes CRPS, MAE, and Winkler scores by comparing simulation draws against
#' observed truth. Uses [scoringRules::crps_sample()] for CRPS computation.
#'
#' By default scores are aggregated per unit. Pass `test_start` to add a
#' `horizon` column (`time - test_start + 1`) and `by` to control the grouping
#' (e.g. `by = c(groupvar, "horizon")` for a per-(unit, horizon) breakdown that
#' lines up with [get_lh_accuracy()] and [compare_approaches()]). `transform`
#' applies a function to both draws and truth before scoring.
#'
#' @param simulation_results A data.table from [simulate_endogenr()], filtered to
#'   the forecast period.
#' @param outcome Character. Name of the outcome variable to score.
#' @param truth A data.frame/data.table with observed values. Must contain the
#'   groupvar, timevar, and outcome columns.
#' @param ctx A `panel_context` object. Optional; inferred from
#'   `simulation_results` first, then `truth`, if either carries
#'   `panel_unit` / `panel_time` attributes (e.g. produced by
#'   [simulate_endogenr()] or `paneltools::as_panel()`).
#' @param sim_var Character. Name of the simulation ID column (default ".sim").
#' @param level Numeric. Coverage level for the Winkler score (default 50).
#' @param test_start Numeric or NULL. When supplied, a `horizon` column is added
#'   as `time - test_start + 1` before aggregation.
#' @param by Character vector of grouping columns for aggregation. Defaults to
#'   the unit column. Including `"horizon"` requires `test_start`.
#' @param transform Function or NULL. Applied to both draws and truth before
#'   scoring (score on a transformed scale).
#'
#' @return A data.table of scores, one row per group in `by`, with columns
#'   `crps`, `mae`, `winkler`.
#' @family postprocess
#' @export
get_accuracy <- function(simulation_results, outcome, truth, ctx = NULL,
                         sim_var = ".sim", level = 50, test_start = NULL,
                         by = NULL, transform = NULL) {
  if (is.null(ctx)) {
    ctx <- .infer_ctx(simulation_results, truth)
    if (is.null(ctx)) {
      stop(
        "Could not infer panel context. Pass `ctx = panel_context(unit, time)` ",
        "or supply a `paneltools::as_panel()` data object.",
        call. = FALSE
      )
    }
  }
  unit_col <- ctx_unit(ctx)
  time_col <- ctx_time(ctx)

  if (!data.table::is.data.table(simulation_results)) {
    simulation_results <- data.table::as.data.table(simulation_results)
  }
  if (!data.table::is.data.table(truth)) {
    truth <- data.table::as.data.table(as.data.frame(truth))
  }

  group_by <- if (is.null(by)) unit_col else by
  if ("horizon" %in% group_by && is.null(test_start)) {
    stop("`by` includes 'horizon' but `test_start` was not supplied.",
         call. = FALSE)
  }

  # Nest draws per (unit, time)
  draws_dt <- simulation_results[, .(draws = list(get(outcome))), by = c(unit_col, time_col)]

  # Merge with truth
  truth_cols <- c(unit_col, time_col, outcome)
  truth_sub <- truth[, ..truth_cols]
  data.table::setnames(truth_sub, outcome, ".truth")
  scored <- merge(draws_dt, truth_sub, by = c(unit_col, time_col))

  # Optional scoring-scale transform applied to both draws and truth
  if (!is.null(transform)) {
    scored[, draws := lapply(draws, transform)]
    scored[, .truth := transform(.truth)]
  }

  # Score per (unit, time) with the shared NA-safe kernel
  sc <- .score_draws(scored$draws, scored$.truth, level)
  scored[, `:=`(crps = sc$crps, mae = sc$mae, winkler = sc$winkler)]

  # Optional horizon dimension
  if (!is.null(test_start)) {
    scored[, horizon := get(time_col) - test_start + 1L]
  }

  # Aggregate
  result <- scored[, .(crps = mean(crps, na.rm = TRUE),
                       mae = mean(mae, na.rm = TRUE),
                       winkler = mean(winkler, na.rm = TRUE)), by = group_by]
  return(result)
}

#' Plot simulations for selected outcome and units
#'
#' Computes quantile ribbons from simulation draws and plots them alongside
#' observed truth data.
#'
#' @param simulation_results A data.table from [simulate_endogenr()], filtered to
#'   the forecast period.
#' @param outcome Character. Name of the outcome variable to plot.
#' @param units Vector of unit IDs to include.
#' @param true_data A data.frame/data.table with observed historical data.
#' @param ctx A `panel_context` object. Optional; inferred from
#'   `simulation_results` first, then `true_data`, if either carries
#'   `panel_unit` / `panel_time` attributes (e.g. produced by
#'   [simulate_endogenr()] or `paneltools::as_panel()`).
#' @param sim_var Character. Name of the simulation ID column (default ".sim").
#'
#' @return A ggplot object.
#' @family postprocess
#' @export
plotsim <- function(simulation_results, outcome, units, true_data, ctx = NULL, sim_var = ".sim") {
  if (is.null(ctx)) {
    ctx <- .infer_ctx(simulation_results, true_data)
    if (is.null(ctx)) {
      stop(
        "Could not infer panel context. Pass `ctx = panel_context(unit, time)` ",
        "or supply a `paneltools::as_panel()` data object.",
        call. = FALSE
      )
    }
  }
  unit_col <- ctx_unit(ctx)
  time_col <- ctx_time(ctx)

  if (!data.table::is.data.table(simulation_results)) {
    simulation_results <- data.table::as.data.table(simulation_results)
  }
  if (!data.table::is.data.table(true_data)) {
    true_data <- data.table::as.data.table(as.data.frame(true_data))
  }

  # Filter to selected units
  sim_sub <- simulation_results[simulation_results[[unit_col]] %in% units]
  truth_sub <- true_data[true_data[[unit_col]] %in% units]

  # Compute quantile ribbons
  probs <- c(0.025, 0.10, 0.25, 0.50, 0.75, 0.90, 0.975)
  ribbons <- sim_sub[, {
    vals <- get(outcome)
    qs <- as.list(stats::quantile(vals, probs = probs, na.rm = TRUE))
    names(qs) <- paste0("q", gsub("\\.", "", sprintf("%.3f", probs)))
    qs
  }, by = c(unit_col, time_col)]

  p <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(data = ribbons,
      ggplot2::aes(x = .data[[time_col]], ymin = q0025, ymax = q0975,
                   group = .data[[unit_col]]),
      fill = "steelblue", alpha = 0.15) +
    ggplot2::geom_ribbon(data = ribbons,
      ggplot2::aes(x = .data[[time_col]], ymin = q0100, ymax = q0900,
                   group = .data[[unit_col]]),
      fill = "steelblue", alpha = 0.2) +
    ggplot2::geom_ribbon(data = ribbons,
      ggplot2::aes(x = .data[[time_col]], ymin = q0250, ymax = q0750,
                   group = .data[[unit_col]]),
      fill = "steelblue", alpha = 0.3) +
    ggplot2::geom_line(data = ribbons,
      ggplot2::aes(x = .data[[time_col]], y = q0500, group = .data[[unit_col]]),
      color = "steelblue", linewidth = 0.8) +
    ggplot2::geom_line(data = truth_sub,
      ggplot2::aes(x = .data[[time_col]], y = .data[[outcome]],
                   group = .data[[unit_col]])) +
    ggplot2::facet_wrap(stats::as.formula(paste("~", unit_col))) +
    ggplot2::labs(y = outcome) +
    ggplot2::theme_minimal()

  return(p)
}
