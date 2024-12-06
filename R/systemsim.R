
#' Prepares simulation data (subset data, expand simulations, include future horizon)
#'
#' @param data
#' @param groupvar
#' @param timevar
#' @param train_start
#' @param test_start
#' @param horizon
#' @param inner_sims
#'
#' @return
#' @export
#'
#' @examples
prepare_simulation_data <- function(data, groupvar, timevar, train_start, test_start, horizon, inner_sims) {
  # Convert groupvar and timevar to symbols for tidy evaluation
  group_sym <- rlang::sym(groupvar)
  time_sym <- rlang::sym(timevar)

  # Filter training data
  simulation_data <- data %>%
    dplyr::filter(!!time_sym < test_start)

  # Create sequence of all times
  all_times <- base::seq(train_start, (test_start + horizon - 1), by = 1)

  # Get unique groups
  all_groups <- simulation_data %>%
    dplyr::distinct(!!group_sym) %>%
    dplyr::pull(!!group_sym)

  # Create expanded grid using tidyr
  expanded <- tidyr::expand_grid(
    !!time_sym := all_times,
    !!group_sym := all_groups,
    sim = 1:inner_sims
  )

  # Convert to tsibble
  result <- expanded %>%
    tsibble::as_tsibble(
      key = c(groupvar, "sim"),
      index = timevar
    ) %>%
    dplyr::left_join(
      simulation_data,
      by = c(groupvar, timevar)
    )

  return(result)
}

#' Imputes the time-independent forecasts
#'
#' @param simulation_data
#' @param models
#' @param test_start
#'
#' @return
#' @export
#'
#' @examples
process_independent_models <- function(simulation_data, models, test_start, horizon, inner_sims) {
  # Input validation
  if (!is.numeric(test_start)) {
    stop("`test_start´ must be numeric")
  }

  # Get key and index variables from tsibble
  grp <- tsibble::key_vars(simulation_data)
  idx <- tsibble::index_var(simulation_data)

  # Filter independent models
  independent_models <- models[sapply(models, function(x) x$independent)]

  # Process each independent model
  for (model in independent_models) {
    # Get predictions
    pred <- predict(model, data = simulation_data, test_start = test_start, horizon = horizon, inner_sims = inner_sims)

    # Update the simulation data with new predictions
    simulation_data <- simulation_data |>
      dplyr::rows_patch(pred, by = c(grp, idx))
  }

  return(simulation_data)
}

#' Dynamic simulation of the time-dependent models
#'
#' @param simulation_data
#' @param models
#' @param test_start
#' @param horizon
#' @param execution_order
#'
#' @return
#' @export
#'
#' @examples
process_dependent_models <- function(simulation_data, models, test_start, horizon, execution_order) {
  # Get dependent models and their outcomes
  dependent_models <- models[sapply(models, function(x) !x$independent)]
  outcomes <- sapply(dependent_models, function(x) parse_formula(x)$outcome)
  names(dependent_models) <- outcomes

  # Order models according to execution order
  execution_order <- execution_order[execution_order %in% outcomes]
  dependent_models <- dependent_models[execution_order]

  # Get key and index variables from tsibble
  grp <- tsibble::key_vars(simulation_data)
  idx <- tsibble::index_var(simulation_data)

  # Process each time point
  for (t in test_start:(test_start + horizon - 1)) {
    # Process each model in the specified order
    for (model in dependent_models) {
      # Get predictions for current time point
      pred <- predict(model, t = t, data = simulation_data)

      # Update simulation data with new predictions using rows_patch
      simulation_data <- simulation_data |>
        dplyr::rows_patch(pred, by = c(grp, idx))
    }
  }

  return(simulation_data)
}


#' Sets up the endogenr simulator, including fitting the models
#'
#' Used as input to [simulate_endogenr()]. The function will impute further parameters into the models defined in [build_model()]
#' using [purrr::partial()].
#'
#' @param models A list of models created using [build_model()]
#' @param data A data.frame. Must contain all variables in the model system, groupvar, timevar, and future exogenous data.
#' @param train_start Integer. When to (at least) start training. The system will remove all data before this time.
#' @param test_start Integer. When to start the forecast period. Forecast includes the test_start time.
#' @param horizon Integer. How many time steps to forecast into.
#' @param groupvar A string. The variable denoting groups in the panel-data.
#' @param timevar A string. The variable denoting time in the panel-data.
#' @param inner_sims Integer. The number of "inner simulations", i.e., simulations using the same estimated models. While a pure bootstrapping
#'   might only want outer simulations, it is much faster to calculate the inner ones, and exploring the uncertainty of any estimated model system
#'   that way.
#' @param min_window Integer. When using linear models (see [build_model()]), if min_window is not NULL, every outer simulation will provide
#'   a random subset time-window of training data as training data for the linear models. This can be useful to explore uncertainties in data-generating
#'   processes that are not stable across time.
#'
#' @return A list of setup parameters for [simulate_endogenr()].
#' @export
#'
#' @examples
#' df <- endogenr::example_data |> tsibble::as_tsibble(key = "gwcode", index = "year")
#' train <- df |> dplyr::filter(year>= 1970, year < 2010) # used for starting values in parametric distribution
#' c1 <- yjbest ~ lag(zoo::rollsumr(yjbest, k = 5, fill = NA)) + lag(log(gdppc))
#' model_system <- list(
#'   build_model("deterministic",formula = gdppc ~ I(abs(lag(gdppc)*(1+gdppc_grwt)))),
#'   build_model("deterministic", formula = gdp ~ I(abs(gdppc*population))),
#'   build_model("parametric_distribution", formula = ~gdppc_grwt, distribution = "t_ls", start = list(df = 1, mu = mean(train$gdppc_grwt), sigma = sd(train$gdppc_grwt))),
#'   build_model("linear", formula = c1, boot = "resid"),
#'   build_model("univariate_fable", formula = dem ~ error("A") + trend("N") + season("N"), method = "ets"),
#'   build_model("exogen", formula = ~psecprop),
#'   build_model("exogen", formula = ~population)
#' )
#'
#' simulator_setup <- setup_simulator(models = model_system,
#'                                   data = df,
#'                                   train_start = 1970,
#'                                   test_start = 2010,
#'                                   horizon = 12,
#'                                   groupvar = "gwcode",
#'                                   timevar = "year",
#'                                   inner_sims = 50,
#'                                   min_window = 10)
setup_simulator <- function(models, data, train_start, test_start, horizon, groupvar, timevar, inner_sims, min_window = NULL){
  data <- data |> dplyr::filter(!!rlang::sym(timevar) >= train_start, !!rlang::sym(timevar) <= (test_start + horizon - 1))
  train <- data |> dplyr::filter(!!rlang::sym(timevar) < test_start)

  models <- lapply(models, function(x){
    model_types <- c("deterministic", "parametric_distribution", "linear", "exogen", "univariate_fable")
    type <- model_types[model_types %in% class(x)]
    if(!is.null(min_window) & type == "linear"){
      type <- "linear_subset"
    }

    f <- switch(type,
           "deterministic" = x,
           "parametric_distribution" = purrr::partial(x, data = train),
           "linear" = purrr::partial(
             x,
             data = train
           ),
           "linear_subset" = purrr::partial(
             x,
             data = train,
             subset = get_train_window(train_start, test_start, min_window)
           ),
           "exogen" = purrr::partial(
             x,
             newdata = data,
             impute_from = test_start,
             inner_sims = inner_sims
           ),
           "univariate_fable" = purrr::partial(
             x,
             data = train
           ),
           stop("Unknown model type: ", type)
    )
    class(f) <- c(class(f), type)
    return(f)
  })

  fitted_models <- lapply(models, function(x) x())

  dependency_graph <-  igraph::make_empty_graph(directed = TRUE)
  for(model in fitted_models){
    dependency_graph <- update_dependency_graph(model, dependency_graph)
  }
  execution_order <- get_execution_order(dependency_graph)

  simulation_data <- prepare_simulation_data(data, groupvar, timevar, train_start, test_start, horizon, inner_sims)

  return(list("simulation_data" = simulation_data,
              "models" = models,
              "test_start" = test_start,
              "horizon" = horizon,
              "execution_order" = execution_order,
              "groupvar" = groupvar,
              "timevar" = timevar,
              "inner_sims" = inner_sims))
}


#' The inner simulation loop
#'
#' @param i
#' @param simulation_data
#' @param models
#' @param test_start
#' @param horizon
#' @param execution_order
#' @param inner_sims
#'
#' @return
#' @export
#'
#' @examples
inner_simulation <- function(i, simulation_data, models, test_start, horizon, execution_order, inner_sims){
  # Fit the linear models every new simulation
  fitted_models <- lapply(models, function(x){
    if(c("linear", "linear_subset") %in% class(x) |> any()){
      return(x())
    } else{
      return(x)
    }
  })
  #fitted_models <- lapply(models, function(x) x()) # fit new
  sim <- process_independent_models(simulation_data, fitted_models, test_start, horizon, inner_sims)
  sim <- process_dependent_models(sim, fitted_models, test_start, horizon, execution_order)
  return(sim)
}

#' Dynamic simulation of the system
#'
#' This is the main function used to run the endogenr simulation system. It will estimate all models,
#' sequence the forecast simulation, and simulate outcomes for all units across the forecast horizon.
#' It supports parallel computation using [future::multisession()] which (should) also work on Windows systems.
#'
#' @param nsim Integer. The number of outer-loop simulations. The outer-loop includes new estimation of models, possibly based on
#'  a subset of training data, or other stochastic elements in the estimation stage. See [setup_simulator()] for details on inner simulations.
#'  The total number of simulations would be nsim * inner_sims (defined in [setup_simulator()]). It is recommended to use a number that has ncores
#'  as a common divisor to effectively utilize all cores.
#' @param simulator_setup A list of simulator parameters generated by [setup_simulator()].
#' @param parallel Boolean. If false, will run sequential on one core. If true, will run using [future::multisession()].
#' @param ncores The number of cores to utilize if running in parallel. Use [parallelly::availableCores()] to find out how many you have available.
#'
#' @return a data.frame with simulation results, as well as the training data. .sim denotes each simulation.
#' @export
#'
#' @examples
#' df <- endogenr::example_data |> tsibble::as_tsibble(key = "gwcode", index = "year")
#' train <- df |> dplyr::filter(year>= 1970, year < 2010) # used for starting values in parametric distribution
#' c1 <- yjbest ~ lag(zoo::rollsumr(yjbest, k = 5, fill = NA)) + lag(log(gdppc))
#' model_system <- list(
#'   build_model("deterministic",formula = gdppc ~ I(abs(lag(gdppc)*(1+gdppc_grwt)))),
#'   build_model("deterministic", formula = gdp ~ I(abs(gdppc*population))),
#'   build_model("parametric_distribution", formula = ~gdppc_grwt, distribution = "t_ls", start = list(df = 1, mu = mean(train$gdppc_grwt), sigma = sd(train$gdppc_grwt))),
#'   build_model("linear", formula = c1, boot = "resid"),
#'   build_model("univariate_fable", formula = dem ~ error("A") + trend("N") + season("N"), method = "ets"),
#'   build_model("exogen", formula = ~psecprop),
#'   build_model("exogen", formula = ~population)
#' )
#'
#' simulator_setup <- setup_simulator(models = model_system,
#'                                   data = df,
#'                                   train_start = 1970,
#'                                   test_start = 2010,
#'                                   horizon = 12,
#'                                   groupvar = "gwcode",
#'                                   timevar = "year",
#'                                   inner_sims = 2,
#'                                   min_window = 10)
#' set.seed(42)
#' res <- simulate_endogenr(nsim = 2, simulator_setup = simulator_setup, parallel = F)
simulate_endogenr <- function(nsim, simulator_setup, parallel = FALSE, ncores = 6){


  old_plan <- future::plan()
  if(parallel){
    future::plan(future::multisession, workers = ncores, gc = TRUE)
  } else{
    future::plan(future::sequential)
  }

  # simulation_results <- lapply(1:nsim, inner_simulation,
  #                                                   simulation_data = simulator_setup$simulation_data,
  #                                                   models = simulator_setup$models,
  #                                                   test_start = simulator_setup$test_start,
  #                                                   horizon = simulator_setup$horizon,
  #                                                   execution_order = simulator_setup$execution_order,
  #                                                   inner_sims = simulator_setup$inner_sims)
  #
  # simulation_results <- future.apply::future_lapply(1:nsim, inner_simulation,
  #                                                   simulation_data = simulator_setup$simulation_data,
  #                                                   models = simulator_setup$models,
  #                                                   test_start = simulator_setup$test_start,
  #                                                   horizon = simulator_setup$horizon,
  #                                                   execution_order = simulator_setup$execution_order,
  #                                                   inner_sims = simulator_setup$inner_sims,
  #                                                   future.seed = TRUE,
  #                                                   future.packages = c("dplyr", "endogenr"))

  # Fit all that only needs to be fit once
  simulator_setup$models <- lapply(simulator_setup$models, function(x){
    if(c("linear", "linear_subset") %in% class(x) |> any()){
      return(x)
    } else{
      return(x())
    }
  })

  simulation_results <- list()
  for(i in 1:nsim){
    simulation_results[[i]] <- future::future({
      inner_simulation(
        i,
        simulation_data = simulator_setup$simulation_data,
        models = simulator_setup$models,
        test_start = simulator_setup$test_start,
        horizon = simulator_setup$horizon,
        execution_order = simulator_setup$execution_order,
        inner_sims = simulator_setup$inner_sims
      )
    }, packages = c("dplyr", "endogenr"), seed = TRUE, globals = "simulator_setup")
  }
  simulation_results <- lapply(simulation_results, FUN = future::value)
  future::plan(old_plan)

  simulation_results <- lapply(simulation_results, dplyr::as_tibble)
  simulation_results <- dplyr::bind_rows(simulation_results, .id = ".id")

  ids <- simulation_results |> dplyr::select(.id, sim) |> unique() |> dplyr::mutate(newid = 1:dplyr::n())
  simulation_results <- dplyr::left_join(simulation_results, ids, by = c(".id", "sim"))
  simulation_results$.id <- NULL
  simulation_results$sim <- NULL
  simulation_results <- simulation_results |> dplyr::rename(.sim = newid)

  return(simulation_results)
}

#' Nest simulation results into distributional::dist_sample
#'
#' @param simulation_results
#' @param outputs
#' @param sim_var
#'
#' @return
#' @export
#'
#' @examples
sim_to_dist = function(simulation_results, outputs, sim_var = ".sim"){
  grp <- tsibble::key_vars(simulation_results)
  grp <- grp[!grp %in% sim_var]
  idx <- tsibble::index_var(simulation_results)
  # Create the summarise expression dynamically
  summarise_expr <- outputs |>
    purrr::map(~ rlang::expr(list(!!rlang::sym(.x)))) |>
    rlang::set_names(outputs)

  # Apply the transformation
  nested_sim <- simulation_results |> dplyr::as_tibble() |>
    dplyr::group_by(!!rlang::sym(idx), !!rlang::sym(grp)) |>
    dplyr::summarise(!!!summarise_expr, .groups = "drop_last") |>
    dplyr::mutate(across(all_of(outputs), ~ distributional::dist_sample(.))) |>
    tsibble::tsibble(key = grp, index = idx)
  return(nested_sim)
}

#' Calculate probabilistic accuracy scores for simulations
#'
#' If you have simulated out-of-sample, i.e., with a true observed outcome for the forecast horizon,
#' you can test the probabilistic accuracy of any simulated outcome using this function. It is essentially
#' a wrapper around [fabletools::accuracy()] and [fabletools::fable()]. As of now, the simulation results from
#' [simulate_endogenr()] returns a tibble. You will need to convert it to a tsibble using c(groupvar, sim_var) and timevar
#' as key and index. You must also subset away the training period. Then it can be used as simulation_results here.
#'
#' The function calculates three different accuracy scores: [fabletools::CRPS()], [fabletools::MAE], and [fabletools::winkler_score(level = 50)].
#' Access to other choices can be implemented, but is not there currently.
#'
#' @param simulation_results A tsibble with a c(groupvar, sim_var) as keys and timevar as index.
#' @param outcome Any one of the other column names in simulation_results.
#' @param truth A tsibble with groupvar as key and timevar as index that includes the observed outcome variable
#'   in the forecast horizon.
#' @param sim_var Defaults to ".sim". Should be the name of the variable denoting the simulation index.
#'
#' @return A tibble with a row per groupvar. Columns: groupvar, .type, crps, mae, winkler.
#' @export
#'
#' @examples
#' df <- endogenr::example_data |> tsibble::as_tsibble(key = "gwcode", index = "year")
#' train <- df |> dplyr::filter(year>= 1970, year < 2010) # used for starting values in parametric distribution
#' c1 <- yjbest ~ lag(zoo::rollsumr(yjbest, k = 5, fill = NA)) + lag(log(gdppc))
#' model_system <- list(
#'   build_model("deterministic",formula = gdppc ~ I(abs(lag(gdppc)*(1+gdppc_grwt)))),
#'   build_model("deterministic", formula = gdp ~ I(abs(gdppc*population))),
#'   build_model("parametric_distribution", formula = ~gdppc_grwt, distribution = "t_ls", start = list(df = 1, mu = mean(train$gdppc_grwt), sigma = sd(train$gdppc_grwt))),
#'   build_model("linear", formula = c1, boot = "resid"),
#'   build_model("univariate_fable", formula = dem ~ error("A") + trend("N") + season("N"), method = "ets"),
#'   build_model("exogen", formula = ~psecprop),
#'   build_model("exogen", formula = ~population)
#' )
#'
#' simulator_setup <- setup_simulator(models = model_system,
#'                                   data = df,
#'                                   train_start = 1970,
#'                                   test_start = 2010,
#'                                   horizon = 12,
#'                                   groupvar = "gwcode",
#'                                   timevar = "year",
#'                                   inner_sims = 2,
#'                                   min_window = 10)
#' set.seed(42)
#' res <- simulate_endogenr(nsim = 2, simulator_setup = simulator_setup, parallel = F)
#' res <- tsibble::tsibble(res, key = c(simulator_setup$groupvar, ".sim"), index = simulator_setup$timevar) |>
#'   dplyr::filter(year >= simulator_setup$test_start)
#' acc <- get_accuracy(res, "gdppc_grwt", df)
#' acc |> summarize(across(crps:winkler, ~ mean(.x))) |> arrange(crps) |> knitr::kable() # return average accuracy across units
get_accuracy <- function(simulation_results, outcome, truth, sim_var = ".sim"){
  grp <- tsibble::key_vars(simulation_results)
  grp <- grp[!grp %in% sim_var]
  idx <- tsibble::index_var(simulation_results)

  forecast <- sim_to_dist(simulation_results, outcome) |>
    fabletools::fable(key = grp, index = idx,  response = outcome, distribution = outcome)

  metrics <- list(crps = fabletools::CRPS, mae = fabletools::MAE, winkler = function(.dist, .actual, ...) fabletools::winkler_score(.dist, .actual, level = 50))
  acc <- forecast |> fabletools::accuracy(truth, metrics)
  return(acc)
  #acc |> summarize(across(crps:winkler, ~ mean(.x))) |> arrange(crps) |> knitr::kable()
}

#' Plot simulations for selected outcome and units
#'
#' As of now, the simulation results from
#' [simulate_endogenr()] returns a tibble. You will need to convert it to a tsibble using c(groupvar, sim_var) and timevar
#' as key and index. You must also subset away the training period. Then it can be used as simulation_results here.
#'
#' @param simulation_results tsibble
#' @param outcome Character string. The name of the variable to plot.
#' @param units A vector of units in groupvar.
#' @param true_data A tsibble with historical data.
#' @param sim_var Defaults to ".sim". Should be the name of the variable denoting the simulation index.
#'
#' @return
#' @export
#'
#' @examples
#' df <- endogenr::example_data |> tsibble::as_tsibble(key = "gwcode", index = "year")
#' train <- df |> dplyr::filter(year>= 1970, year < 2010) # used for starting values in parametric distribution
#' c1 <- yjbest ~ lag(zoo::rollsumr(yjbest, k = 5, fill = NA)) + lag(log(gdppc))
#' model_system <- list(
#'   build_model("deterministic",formula = gdppc ~ I(abs(lag(gdppc)*(1+gdppc_grwt)))),
#'   build_model("deterministic", formula = gdp ~ I(abs(gdppc*population))),
#'   build_model("parametric_distribution", formula = ~gdppc_grwt, distribution = "t_ls", start = list(df = 1, mu = mean(train$gdppc_grwt), sigma = sd(train$gdppc_grwt))),
#'   build_model("linear", formula = c1, boot = "resid"),
#'   build_model("univariate_fable", formula = dem ~ error("A") + trend("N") + season("N"), method = "ets"),
#'   build_model("exogen", formula = ~psecprop),
#'   build_model("exogen", formula = ~population)
#' )
#'
#' simulator_setup <- setup_simulator(models = model_system,
#'                                   data = df,
#'                                   train_start = 1970,
#'                                   test_start = 2010,
#'                                   horizon = 12,
#'                                   groupvar = "gwcode",
#'                                   timevar = "year",
#'                                   inner_sims = 2,
#'                                   min_window = 10)
#' set.seed(42)
#' res <- simulate_endogenr(nsim = 2, simulator_setup = simulator_setup, parallel = F)
#' res <- tsibble::tsibble(res, key = c(simulator_setup$groupvar, ".sim"), index = simulator_setup$timevar) |>
#'   dplyr::filter(year >= simulator_setup$test_start)
#' plotsim(res, "gdppc", c(2, 20, 530), df) # the numbers are Gleditsch-Ward country-codes (USA, Canada, and Ethiopia).
plotsim <- function(simulation_results, outcome, units, true_data, sim_var = ".sim"){
  # Get key and index variables from tsibble
  grp <- tsibble::key_vars(simulation_results)
  grp <- grp[!grp %in% sim_var]
  idx <- tsibble::index_var(simulation_results)
  group_sym <- rlang::sym(grp)

  simulation_results <- simulation_results |>
    dplyr::filter(!!group_sym %in% units)

  myplot <- function(){
    sim_to_dist(simulation_results, outcome) |>
      fabletools::fable(key = grp, index = idx,  response = outcome, distribution = outcome) |>
      fabletools::autoplot(true_data, level = base::seq(5, 95, 5), point_forecast = list("median" = median), alpha = 0.4) +
      ggplot2::scale_y_continuous(labels = scales::comma)}

  suppressWarnings(myplot()) # Warning in fabletools::fable
}
