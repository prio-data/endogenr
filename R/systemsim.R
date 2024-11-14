
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


#' Sets up the simulator, including fitting the models
#'
#' @param models
#' @param data
#' @param train_start
#' @param test_start
#' @param horizon
#' @param groupvar
#' @param timevar
#' @param inner_sims
#'
#' @return
#' @export
#'
#' @examples
setup_simulator <- function(models, data, train_start, test_start, horizon, groupvar, timevar, inner_sims){
  data <- data |> dplyr::filter(!!rlang::sym(timevar) >= train_start, !!rlang::sym(timevar) <= (test_start + horizon - 1))
  train <- data |> dplyr::filter(!!rlang::sym(timevar) < test_start)

  models <- lapply(models, function(x){
    model_types <- c("deterministic", "parametric_distribution", "linear", "exogen", "univariate_fable")
    type <- model_types[model_types %in% class(x)]

    switch(type,
           "deterministic" = x,
           "parametric_distribution" = purrr::partial(x, data = train),
           "linear" = purrr::partial(
             x,
             data = train
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


#' Dynamic simulation of the system
#'
#' @param nsim
#' @param simulator_setup
#' @param parallel
#' @param ncores
#'
#' @return
#' @export
#'
#' @examples
simulate_endogenr <- function(nsim, simulator_setup, parallel = FALSE, ncores = 6){
  simulate <- function(i, simulation_data, models, test_start, horizon, execution_order, inner_sims){
    fitted_models <- lapply(models, function(x) x()) # fit new
    sim <- process_independent_models(simulation_data, fitted_models, test_start, horizon, inner_sims)
    sim <- process_dependent_models(sim, fitted_models, test_start, horizon, execution_order)
    return(sim)
  }

  if(parallel){
    old_plan <- future::plan()
    future::plan(future::multisession, workers = ncores, gc = TRUE)

    simulation_results <- future.apply::future_lapply(1:nsim, simulate,
                                                      simulation_data = simulator_setup$simulation_data,
                                                      models = simulator_setup$models,
                                                      test_start = simulator_setup$test_start,
                                                      horizon = simulator_setup$horizon,
                                                      execution_order = simulator_setup$execution_order,
                                                      inner_sims = simulator_setup$inner_sims,
                                                      future.seed = TRUE,
                                                      future.packages = c("dplyr", "endogenr"))
    future::plan(old_plan)
  } else{
    simulation_results <- lapply(1:nsim, simulate,
                                 simulation_data = simulator_setup$simulation_data,
                                 models = simulator_setup$models,
                                 test_start = simulator_setup$test_start,
                                 horizon = simulator_setup$horizon,
                                 execution_order = simulator_setup$execution_order,
                                 inner_sims = simulator_setup$inner_sims)
  }

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
#' @param simulation_results
#' @param outcome
#' @param truth
#' @param sim_var
#'
#' @return
#' @export
#'
#' @examples
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
#' @param simulation_results
#' @param outcome
#' @param units
#' @param true_data
#' @param sim_var
#'
#' @return
#' @export
#'
#' @examples
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
