
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
#' @export
prepare_simulation_data <- function(data, ctx, train_start, test_start, horizon, inner_sims) {
  unit_col <- ctx_unit(ctx)
  time_col <- ctx_time(ctx)

  # Training data only (data is already filtered to train_start..test_start+horizon-1)
  train <- data[data[[time_col]] < test_start]

  # Create CJ grid
  all_units <- unique(train[[unit_col]])
  all_times <- seq(train_start, test_start + horizon - 1)

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
#' @export
process_independent_models <- function(simulation_data, models, ctx, test_start, horizon, inner_sims) {
  if (!is.numeric(test_start)) {
    stop("`test_start` must be numeric")
  }

  join_keys <- c(ctx_keys(ctx), ctx_time(ctx))

  independent_models <- models[vapply(models, function(x) x$independent, logical(1))]

  for (model in independent_models) {
    pred <- predict(model, data = simulation_data, ctx = ctx,
                    test_start = test_start, horizon = horizon, inner_sims = inner_sims)

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
#' @export
process_dependent_models <- function(simulation_data, models, ctx, test_start, horizon, execution_order) {
  dependent_models <- models[vapply(models, function(x) !x$independent, logical(1))]
  outcomes <- vapply(dependent_models, function(x) parse_formula(x)$outcome, character(1))
  names(dependent_models) <- outcomes

  execution_order <- execution_order[execution_order %in% outcomes]
  dependent_models <- dependent_models[execution_order]

  join_keys <- c(ctx_keys(ctx), ctx_time(ctx))

  for (t in test_start:(test_start + horizon - 1)) {
    for (model in dependent_models) {
      pred <- predict(model, t = t, data = simulation_data, ctx = ctx)
      cols <- setdiff(names(pred), join_keys)
      simulation_data[pred, (cols) := mget(paste0("i.", cols)), on = join_keys]
    }
  }

  return(simulation_data)
}


#' Sets up the endogenr simulator, including fitting the models
#'
#' @param models A list of models created using [build_model()]
#' @param data A data.frame, data.table, or tsibble. Must contain all variables referenced by
#'   model formulas, plus groupvar, timevar, and future exogenous data.
#' @param train_start Integer. The earliest training time to include.
#' @param test_start Integer. When to start the forecast period.
#' @param horizon Integer. How many time steps to forecast.
#' @param groupvar Character. The column name identifying panel units.
#' @param timevar Character. The column name identifying time periods.
#' @param inner_sims Integer. Number of inner simulations per outer sim.
#' @param min_window Integer or NULL. If set, linear/glm/heterolm models get a random
#'   training window of at least this length per outer simulation.
#' @param globals Named list of user functions to export to parallel workers.
#'
#' @return A list of setup parameters for [simulate_endogenr()].
#' @export
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
  train <- data[data[[timevar]] < test_start]

  # Partially apply model constructors with data
  models <- lapply(models, function(x) {
    model_types <- c("deterministic", "parametric_distribution", "linear", "glm",
                     "exogen", "univariate_fable", "heterolm", "spatial_lag")
    type <- model_types[model_types %in% class(x)]
    if (!is.null(min_window) && type == "linear") type <- "linear_subset"
    if (!is.null(min_window) && type == "glm") type <- "glm_subset"
    if (!is.null(min_window) && type == "heterolm") type <- "heterolm_subset"

    f <- switch(type,
      "deterministic" = purrr::partial(x, ctx = ctx),
      "parametric_distribution" = purrr::partial(x, data = train, ctx = ctx),
      "linear" = purrr::partial(x, data = train, ctx = ctx),
      "linear_subset" = purrr::partial(x, data = train, ctx = ctx,
                                       subset = get_train_window(train_start, test_start, min_window)),
      "exogen" = purrr::partial(x, newdata = data, impute_from = test_start,
                                inner_sims = inner_sims, ctx = ctx),
      "univariate_fable" = purrr::partial(x, data = train, ctx = ctx),
      "heterolm" = purrr::partial(x, data = train, ctx = ctx),
      "heterolm_subset" = purrr::partial(x, data = train, ctx = ctx,
                                         subset = get_train_window(train_start, test_start, min_window)),
      "glm" = purrr::partial(x, data = train, ctx = ctx),
      "glm_subset" = purrr::partial(x, data = train, ctx = ctx,
                                    subset = get_train_window(train_start, test_start, min_window)),
      "spatial_lag" = purrr::partial(x, ctx = ctx),
      stop("Unknown model type: ", type)
    )
    class(f) <- c(class(f), type)
    return(f)
  })

  fitted_models <- lapply(models, function(x) x())

  # Build dependency graph
  dependency_graph <- igraph::make_empty_graph(directed = TRUE)
  for (model in fitted_models) {
    dependency_graph <- update_dependency_graph(model, dependency_graph)
  }
  execution_order <- get_execution_order(dependency_graph)

  # Create a simulation context WITH sim dimension
  # IMPORTANT: don't mutate the original ctx — it's shared by model partials
  # which need a ctx without sim for refitting on train data.
  sim_ctx <- panel_context(unit = groupvar, time = timevar, sim = "sim")

  simulation_data <- prepare_simulation_data(data, sim_ctx, train_start, test_start, horizon, inner_sims)

  return(list(
    "simulation_data" = simulation_data,
    "models" = models,
    "fitted_models" = fitted_models,
    "train_data" = train,
    "full_data" = data,
    "train_start" = train_start,
    "test_start" = test_start,
    "horizon" = horizon,
    "execution_order" = execution_order,
    "ctx" = sim_ctx,
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
#' @param simulation_data A data.table.
#' @param models List of (partially applied) model closures.
#' @param ctx A panel_context object.
#' @param test_start Integer.
#' @param horizon Integer.
#' @param execution_order Character vector.
#' @param inner_sims Integer.
#'
#' @return A data.table with simulation results.
#' @export
inner_simulation <- function(i, simulation_data, models, ctx, test_start, horizon,
                             execution_order, inner_sims) {
  # Fit linear/glm/heterolm models every outer sim (random training window)
  fitted_models <- lapply(models, function(x) {
    linear_types <- c("linear", "linear_subset", "glm", "glm_subset",
                      "heterolm", "heterolm_subset")
    if (any(linear_types %in% class(x))) {
      return(x())
    } else {
      return(x)
    }
  })

  sim <- data.table::copy(simulation_data)
  sim <- process_independent_models(sim, fitted_models, ctx, test_start, horizon, inner_sims)
  sim <- process_dependent_models(sim, fitted_models, ctx, test_start, horizon, execution_order)
  return(sim)
}

#' Dynamic simulation of the system
#'
#' Runs the endogenr simulation system. Estimates models, sequences the forecast,
#' and simulates outcomes for all units across the forecast horizon.
#'
#' @param nsim Integer. Number of outer-loop simulations.
#' @param simulator_setup A list from [setup_simulator()].
#' @param parallel Logical. Use parallel computation via [future::multisession()].
#' @param ncores Integer. Number of workers for parallel execution.
#'
#' @return A data.table with all original columns plus `.sim` (integer 1..nsim*inner_sims).
#' @export
simulate_endogenr <- function(nsim, simulator_setup, parallel = FALSE, ncores = 6) {
  old_plan <- future::plan()
  if (parallel) {
    future::plan(future::multisession, workers = ncores, gc = TRUE)
  } else {
    future::plan(future::sequential)
  }

  ctx <- simulator_setup$ctx

  # Pre-fit non-linear models (only need fitting once)
  simulator_setup$models <- lapply(simulator_setup$models, function(x) {
    linear_types <- c("linear", "linear_subset", "glm", "glm_subset",
                      "heterolm", "heterolm_subset")
    if (any(linear_types %in% class(x))) {
      return(x)
    } else {
      return(x())
    }
  })

  # Build globals for future
  future_globals <- list(simulator_setup = simulator_setup)
  if (!is.null(simulator_setup$globals)) {
    for (fn_name in names(simulator_setup$globals)) {
      future_globals[[fn_name]] <- simulator_setup$globals[[fn_name]]
    }
  }

  simulation_results <- vector("list", nsim)
  for (i in seq_len(nsim)) {
    simulation_results[[i]] <- future::future({
      inner_simulation(
        i,
        simulation_data = simulator_setup$simulation_data,
        models = simulator_setup$models,
        ctx = simulator_setup$ctx,
        test_start = simulator_setup$test_start,
        horizon = simulator_setup$horizon,
        execution_order = simulator_setup$execution_order,
        inner_sims = simulator_setup$inner_sims
      )
    }, packages = c("endogenr", "data.table"), seed = TRUE, globals = future_globals)
  }
  simulation_results <- lapply(simulation_results, FUN = future::value)
  future::plan(old_plan)

  # Bind results and create .sim ID
  results <- data.table::rbindlist(simulation_results, idcol = ".id")

  # Create unique .sim across (outer_id, inner_sim) pairs
  ids <- unique(results[, .(.id, sim)])
  ids[, .sim := .I]
  results <- merge(results, ids, by = c(".id", "sim"))
  results[, c(".id", "sim") := NULL]

  return(results)
}

#' Nest simulation draws into list-columns
#'
#' Groups simulation results by (unit, time) and collapses each output variable
#' into a list of draws per cell.
#'
#' @param simulation_results A data.table from [simulate_endogenr()].
#' @param outputs Character vector of column names to nest.
#' @param ctx A panel_context object.
#' @param sim_var Character. Name of the simulation ID column (default ".sim").
#'
#' @return A data.table with one row per (unit, time), output columns as list-columns.
#' @export
sim_to_dist <- function(simulation_results, outputs, ctx, sim_var = ".sim") {
  keys <- c(ctx_unit(ctx), ctx_time(ctx))
  dt <- data.table::as.data.table(simulation_results)
  dt[, lapply(.SD, list), by = keys, .SDcols = outputs]
}

#' Calculate probabilistic accuracy scores for simulations
#'
#' Computes CRPS, MAE, and Winkler scores by comparing simulation draws against
#' observed truth. Uses [scoringRules::crps_sample()] for CRPS computation.
#'
#' @param simulation_results A data.table from [simulate_endogenr()], filtered to
#'   the forecast period.
#' @param outcome Character. Name of the outcome variable to score.
#' @param truth A data.frame/data.table with observed values. Must contain the
#'   groupvar, timevar, and outcome columns.
#' @param ctx A panel_context object.
#' @param sim_var Character. Name of the simulation ID column (default ".sim").
#' @param level Numeric. Coverage level for the Winkler score (default 50).
#'
#' @return A data.table with one row per unit, columns: unit, crps, mae, winkler.
#' @export
get_accuracy <- function(simulation_results, outcome, truth, ctx, sim_var = ".sim", level = 50) {
  unit_col <- ctx_unit(ctx)
  time_col <- ctx_time(ctx)

  if (!data.table::is.data.table(simulation_results)) {
    simulation_results <- data.table::as.data.table(simulation_results)
  }
  if (!data.table::is.data.table(truth)) {
    truth <- data.table::as.data.table(as.data.frame(truth))
  }

  # Nest draws per (unit, time)
  draws_dt <- simulation_results[, .(draws = list(get(outcome))), by = c(unit_col, time_col)]

  # Merge with truth
  truth_cols <- c(unit_col, time_col, outcome)
  truth_sub <- truth[, ..truth_cols]
  data.table::setnames(truth_sub, outcome, ".truth")
  scored <- merge(draws_dt, truth_sub, by = c(unit_col, time_col))

  # Compute scores per (unit, time)
  alpha <- 1 - level / 100
  scored[, `:=`(
    crps = vapply(seq_len(.N), function(i) {
      d <- draws[[i]][!is.na(draws[[i]])]
      if (length(d) == 0 || is.na(.truth[i])) return(NA_real_)
      scoringRules::crps_sample(y = .truth[i], dat = d)
    }, numeric(1)),
    mae = vapply(seq_len(.N), function(i) {
      d <- draws[[i]][!is.na(draws[[i]])]
      if (length(d) == 0 || is.na(.truth[i])) return(NA_real_)
      abs(stats::median(d) - .truth[i])
    }, numeric(1)),
    winkler = vapply(seq_len(.N), function(i) {
      d <- draws[[i]][!is.na(draws[[i]])]
      if (length(d) == 0 || is.na(.truth[i])) return(NA_real_)
      lo <- stats::quantile(d, alpha / 2)
      hi <- stats::quantile(d, 1 - alpha / 2)
      width <- hi - lo
      penalty <- ifelse(.truth[i] < lo, (2 / alpha) * (lo - .truth[i]),
                 ifelse(.truth[i] > hi, (2 / alpha) * (.truth[i] - hi), 0))
      width + penalty
    }, numeric(1))
  )]

  # Aggregate per unit
  result <- scored[, .(crps = mean(crps, na.rm = TRUE),
                       mae = mean(mae, na.rm = TRUE),
                       winkler = mean(winkler, na.rm = TRUE)), by = unit_col]
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
#' @param ctx A panel_context object.
#' @param sim_var Character. Name of the simulation ID column (default ".sim").
#'
#' @return A ggplot object.
#' @export
plotsim <- function(simulation_results, outcome, units, true_data, ctx, sim_var = ".sim") {
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
