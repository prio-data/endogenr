#' Endogenous System Class
#'
#' @description
#' An R6 Class that manages a system of interdependent models for simulation purposes.
#' The EndogenousSystem coordinates model fitting, simulation execution (both sequential
#' and parallel), and result visualization.
#'
#' @details
#' The EndogenousSystem class provides a framework for handling complex systems where
#' multiple variables interact and depend on each other. It manages the simulation
#' process and provides tools for analyzing and visualizing results.
#'
#' @export
EndogenousSystem <- R6::R6Class(
  "EndogenousSystem",
  public = list(
    #' @field simulator The EndogenousSimulator instance managing the core simulation
    simulator = NULL,

    #' @field simulated_data Data.table containing simulation results
    simulated_data = NULL,

    #' @field simulated_models List of fitted models used in simulations
    simulated_models = NULL,

    #' @description
    #' Create a new endogenous system
    #' @param data A data.frame or data.table containing the input data
    #' @param groupvar Character string specifying the grouping variable column name
    #' @param timevar Character string specifying the time variable column name
    #' @param test_start Numeric value indicating when the test period begins
    #' @param train_start Numeric value indicating when the training period begins
    #' @param min_window Numeric value specifying the minimum window size for bootstrap sampling
    #' @return A new `EndogenousSystem` object
    #' @examples
    #' \dontrun{
    #' system <- EndogenousSystem$new(
    #'   data = my_data,
    #'   groupvar = "group",
    #'   timevar = "time",
    #'   test_start = 100,
    #'   train_start = 1
    #' )
    #' }
    initialize = function(data, groupvar, timevar, test_start, train_start, min_window = NULL) {
      self$simulator <- EndogenousSimulator$new(data, groupvar, timevar, test_start, train_start, min_window)
    },

    #' @description
    #' Add a model to the system
    #' @param model A Model object to add to the system
    #' @details Models must not create circular dependencies when added
    add_model = function(model){
      self$simulator$add_model(model)
    },

    #' @description
    #' Fit specified models or all models in the system
    #' @param models Optional list of specific models to fit
    fit = function(models = NULL){
      self$simulator$fit(models)
    },

    #' @description
    #' Run simulations either sequentially or in parallel
    #' @param horizon Number of time periods to simulate
    #' @param nsim Number of simulations to run
    #' @param result_folder Optional path to save results
    #' @param parallel Whether to use parallel processing
    #' @param ncores Number of cores for parallel processing
    #' @param seed Optional random seed for reproducibility
    #' @examples
    #' \dontrun{
    #' system$simulate(
    #'   horizon = 10,
    #'   nsim = 100,
    #'   parallel = TRUE,
    #'   ncores = 4,
    #'   seed = 123
    #' )
    #' }
    simulate = function(horizon, nsim, result_folder = NULL, parallel = TRUE, ncores = 6, seed = NULL){
      if(!is.null(seed)){
        set.seed(seed)
      }

      results <- if (parallel) {
        private$run_parallel_simulation(horizon, nsim, ncores)
      } else {
        private$run_sequential_simulation(horizon, nsim)
      }

      if (!is.null(result_folder)) {
        private$save_results(results, result_folder)
      }

      self$simulated_data <- results$data
      self$simulated_models <- results$model

    },

    #' @description
    #' Visualize simulation results for specific variables and units
    #' @param outcome Character string specifying the outcome variable to plot
    #' @param units Vector of unit identifiers to include in the plot
    #' @return A ggplot2 object showing simulation results with uncertainty bands
    #' @examples
    #' \dontrun{
    #' system$plot("revenue", units = c("unit1", "unit2"))
    #' }
    plot = function(outcome, units){
      grp <- data.table::key(self$simulator$data)
      idx <- data.table::indices(self$simulator$data)

      true_data <- self$simulator$data[group_var %in% units, env = list(group_var = grp)]
      true_data <- suppressMessages(tsibble::tsibble(true_data, key = grp, index = idx)) # Depreciation warning in tsibble function.

      simulated_data <- self$simulated_data[group_var %in% units, env = list(group_var = grp)]

      data.table::setkeyv(simulated_data, grp)
      data.table::setindexv(simulated_data, idx)

      myplot <- function(){
        private$sim_to_dist(simulated_data, outcome) |>
          fabletools::fable(key = grp, index = idx,  response = outcome, distribution = outcome) |>
          fabletools::autoplot(true_data, level = base::seq(5, 95, 5), point_forecast = list("median" = median), alpha = 0.4) +
          ggplot2::scale_y_continuous(labels = scales::comma)}

      suppressWarnings(myplot()) # Warning in fabletools::fable
    }),

  private = list(
    run_sequential_simulation = function(horizon, nsim) {
      res <- pbapply::pblapply(1:nsim, function(x) {
        self$simulator$simulate(horizon)
        return(list("data" = self$simulator$simulation_data, "models" = self$simulator$models))
      })
      private$process_simulation_results(res)
    },

    run_parallel_simulation = function(horizon, nsim, ncores) {


      #### Multisession ####
      #require(future)
      #require(progressr)
      #require(future.apply)
      old_plan <- future::plan()
      future::plan(future::multisession, workers = ncores, gc = TRUE)

      progressr::handlers(global = TRUE)
      progressr::handlers("debug")

      parallel_simulate <- function(simulators, horizon) {
        progressr::with_progress({
          p <- progressr::progressor(along = simulators)
          future.apply::future_lapply(simulators, function(simulator, horizon){
            simulator$simulate(horizon)
            results <- list("data" = data.table::copy(simulator$simulation_data), "models" = simulator$models)
            p()

            results
          }, horizon = horizon, future.seed = TRUE, future.packages = c("dplyr", "zoo"))
        })
      }

      simulators <- lapply(1:nsim, function(x) simulator$simulator$clone(deep = TRUE))
      results <- parallel_simulate(simulators, horizon = horizon)

      progressr::handlers(global = FALSE)
      future::plan(old_plan)
      #####################

      results <- private$process_simulation_results(results)
      return(results)
    },

    create_simulation_chunks = function(nsim, ncores) {
      chunk_size <- ceiling(nsim/ncores)
      chunks <- split(1:nsim, ceiling(base::seq_along(1:nsim)/chunk_size))
      chunks <- lapply(chunks, function(x) 1:length(x))
    },

    sim_to_dist = function(simulation_results, outputs, sim_var = ".id"){
      idx <- data.table::indices(simulation_results)
      grp <- data.table::key(simulation_results)

      tsdat <- tsibble::tsibble(simulation_results, key = c(".id", grp), index = idx)

      # Create the summarise expression dynamically
      summarise_expr <- outputs |>
        purrr::map(~ rlang::expr(list(!!rlang::sym(.x)))) |>
        rlang::set_names(outputs)

      # Apply the transformation
      nested_sim <- tsdat |> dplyr::as_tibble() |>
        dplyr::group_by(!!rlang::sym(idx), !!rlang::sym(grp)) |>
        dplyr::summarise(!!!summarise_expr, .groups = "drop_last") |>
        dplyr::mutate(across(all_of(outputs), ~ distributional::dist_sample(.))) |>
        tsibble::tsibble(key = grp, index = idx)

      return(nested_sim)
    },

    process_simulation_results = function(simulation_results){
      idx <- self$simulator$timevar
      grp <- self$simulator$groupvar

      # simulation_results is a list(list(simulation_data, modelobj))
      simulated_data <- lapply(simulation_results, function(x) x$data)
      simulated_model <- lapply(simulation_results, function(x) x$models)
      # Flatten the results and collect
      simulated_data <- data.table::rbindlist(simulated_data, idcol = ".id")

      simulated_data <- simulated_data[time >= self$simulator$test_start, env = list(time = idx)]
      data.table::setkeyv(simulated_data, c(".id", grp))
      data.table::setindexv(simulated_data, idx)
      return(list("model" = simulated_model, "data" = simulated_data))
    }
  )

)


#' Endogenous Simulator Class
#'
#' @description
#' An R6 Class that handles the core simulation logic for systems of interdependent
#' models. It manages model dependencies, execution order, and the actual simulation
#' process.
#'
#' @details
#' The EndogenousSimulator is responsible for:
#' \itemize{
#'   \item Validating and managing model dependencies
#'   \item Determining the correct order of model execution
#'   \item Handling both exogenous and endogenous variables
#'   \item Managing the simulation state through time
#' }
#'
#' The simulator ensures that models are executed in the correct order based on their
#' dependencies and handles both exogenous and endogenous variables appropriately.
#'
EndogenousSimulator <- R6::R6Class(
  "EndogenousSimulator",
  public = list(
    #' @field data The full dataset used for simulation
    data = NULL,

    #' @field train Training dataset subset
    train = NULL,

    #' @field groupvar Character string specifying the grouping variable
    groupvar = NULL,

    #' @field timevar Character string specifying the time variable
    timevar = NULL,

    #' @field test_start Numeric indicating start of test period
    test_start = NULL,

    #' @field train_start Numeric indicating start of training period
    train_start = NULL,

    #' @field min_window Minimum window size for bootstrap sampling
    min_window = NULL,

    #' @field time_parameters List of time-related parameters
    time_parameters = NULL,

    #' @field models List of models in the system
    models = NULL,

    #' @field dependency_graph igraph object representing model dependencies
    dependency_graph = NULL,

    #' @field outputs Character vector of output variable names
    outputs = NULL,

    #' @field simulation_data Data.table containing current simulation results
    simulation_data = NULL,

    #' @description
    #' Create a new endogenous simulator instance
    #' @param data A data.frame or data.table containing the input data
    #' @param groupvar Character string specifying the grouping variable column name
    #' @param timevar Character string specifying the time variable column name
    #' @param test_start Numeric value indicating when the test period begins
    #' @param train_start Numeric value indicating when the training period begins
    #' @param min_window Numeric value specifying minimum window size for bootstrap
    #' @return A new `EndogenousSimulator` object
    initialize = function(data, groupvar, timevar, test_start, train_start, min_window = NULL) {
      self$models <- list()
      self$outputs <- character(0)
      self$dependency_graph <-  igraph::make_empty_graph(directed = TRUE)

      private$validate_initialization(data, groupvar, timevar, test_start, train_start)
      private$setup_simulator(data, groupvar, timevar, test_start, train_start, min_window)
    },

    #' @description
    #' Add a model to the simulator and update dependency graph
    #' @param model Model object to add to the system
    #' @return Invisible self (for method chaining)
    add_model = function(model) {
      private$validate_model(model)
      private$register_model(model)
      private$update_dependency_graph(model)
      invisible(self)
    },

    #' @description
    #' Fit all models using the training data
    fit = function() {
      for (model in self$models) {
        model$fit(self$train, self$time_parameters)
      }
    },

    #' @description
    #' Run a single simulation for the specified horizon
    #' @param horizon Number of time periods to simulate
    simulate = function(horizon) {
      self$time_parameters <- private$calculate_time_window(self$train_start, self$test_start, self$min_window)

      # Setup simulation data frame
      self$simulation_data <- NULL # Make sure this is not there from before.
      self$simulation_data <- private$prepare_simulation_data(horizon)

      # Update model fits
      self$fit()

      # Find which models to calculate when.
      execution_order <- private$get_execution_order()

      # Process exogenous variables
      self$simulation_data <- private$process_exogenous_variables(self$simulation_data, self$models, execution_order)
      # Process endogenous variables
      self$simulation_data <- private$process_endogenous_variables(self$simulation_data, self$models, execution_order)
    },

    #' @description
    #' Get current state of the simulation for debugging
    #' @return List containing data summary, model summary, and graph summary
    get_simulation_state = function() {
      list(
        data_summary = private$summarize_data(),
        model_summary = private$summarize_models(),
        graph_summary = private$summarize_dependency_graph()
      )
    }
  ),

  private = list(
    deep_clone = function(name, value){

      if(name == "simulation_data"){
        if (is.null(value)) return(NULL)
        # These objects will be mutated during a simulation.
        if(inherits(value, "data.table")){
          return(data.table::copy(value))
        } else if(is.list(value)){
          y <- list(value) |> unlist(recursive = FALSE)
          if(tracemem(y) == tracemem(value)){
            stop("List was not copied")
          }
          return(y)
        } else{
          stop("Not implemented deep copy for this simulation_data object.")
        }
      }
      if(name == "models"){
        if (is.null(value)) return(NULL)
        if(is.list(value)){
          y <- list(value) |> unlist(recursive = FALSE)
          if(tracemem(y) == tracemem(value)){
            stop("List was not copied")
          }
          return(y)
        } else{
          stop("Not implemented deep copy for this model object.")
        }
      }

      # For all other objects
      return(value)
    },

    validate_initialization = function(data, groupvar, timevar, test_start, train_start) {
      if (!data.table::is.data.table(data)) {
        data <- data.table::as.data.table(data)
      }

      if (train_start < data[, min(time), env = list(time = timevar)]) {
        stop("No data for time of train start. Please change train_start or add data.")
      }

      if (!all(c(groupvar, timevar) %in% names(data))) {
        stop("groupvar and timevar must be columns in the data")
      }
    },

    #### Model management methods ####
    validate_model = function(model) {
      if (!base::inherits(model, "Model")) {
        stop("Model must inherit from the Model class")
      }

      if ("ExogenModel" %in% class(model)) {
        if (!identical(data.table::key(model$exogendata), data.table::key(self$train))) {
          stop("Exogen input data must have the same key as groupvar")
        }
        if (!identical(data.table::indices(model$exogendata), data.table::indices(self$train))) {
          stop("Exogen input data must have timevar as its index")
        }
      }

      if (any(c("ExogenModel", "StochasticStaticModel") %in% class(model))) {
        output <- base::all.vars(model$formula)
        if (length(output) != 1) {
          stop("ExogenModels and StochasticStaticModels should only contain one variable in formula")
        }
      }
    },

    register_model = function(model) {
      self$models[[model$name]] <- model
    },

    update_dependency_graph = function(model) {
      if (!any(c("ExogenModel", "StochasticStaticModel") %in% class(model))) {
        edges <- private$parse_model_formula(model$formula)
        self$outputs <- unique(c(self$outputs, edges$out))
        vertices <- unique(unlist(edges))
        private$add_to_graph(vertices, edges)
      }

      if (any(c("ExogenModel", "StochasticStaticModel") %in% class(model))) {
        vertex <- base::all.vars(model$formula)
        self$outputs <- unique(c(self$outputs, vertex))
        edges <- cbind("in" = vertex, "out" = vertex)
        private$add_to_graph(vertex, edges)
      }
    },

    add_to_graph = function(vertices, edges) {
      # Add new vertices if they don't exist
      existing_vertices <- igraph::V(self$dependency_graph)$name
      new_vertices <- base::setdiff(vertices, existing_vertices)

      if (length(new_vertices) > 0) {
        self$dependency_graph <- igraph::add_vertices(
          self$dependency_graph,
          length(new_vertices),
          name = new_vertices
        )
      }

      # Add edges
      edge_list <- as.matrix(edges)
      self$dependency_graph <- igraph::add_edges(
        self$dependency_graph,
        t(edge_list)
      )
    },

    get_execution_order = function() {
      second_edges <- igraph::E(self$dependency_graph)[.to(igraph::V(self$dependency_graph)[self$outputs])]
      second <- self$dependency_graph |> igraph::delete_edges(second_edges)

      # Use topological sort to determine calculation order
      second_order <- igraph::topo_sort(second, mode = "out")
      execution_order <- igraph::V(self$dependency_graph)$name[second_order]
      return(execution_order[execution_order %in% names(self$models)])
    },

    parse_model_formula = function(formula) {
      terms <- base::all.vars(rlang::f_rhs(formula))
      outcome <- base::all.vars(rlang::f_lhs(formula))
      edges <- data.frame("in" = terms, "out" = outcome)
      return(edges)
    },

    #### Simulation methods ####
    validate_simulation_params = function(horizon, nsim, ncores) {
      if (!is.numeric(horizon) || horizon < 1) {
        stop("horizon must be a positive integer")
      }
      if (!is.numeric(nsim) || nsim < 1) {
        stop("nsim must be a positive integer")
      }
      if (ncores > parallel::detectCores()) {
        warning("Requested cores exceed available cores, using maximum available")
      }
    },

    setup_simulator = function(data, groupvar, timevar, test_start, train_start, min_window) {
      data.table::setkeyv(data, groupvar)
      data.table::setindexv(data, timevar)

      train <- data[time < test_start & time >= train_start, env = list(time = timevar)]
      data.table::setkeyv(train, groupvar)
      data.table::setindexv(train, timevar)

      self$data <- data
      self$train <- train
      self$groupvar <- groupvar
      self$timevar <- timevar
      self$test_start <- test_start
      self$train_start <- train_start
      self$min_window <- min_window
      self$dependency_graph <- igraph::make_empty_graph(directed = TRUE)
      self$outputs <- character(0)
    },

    calculate_time_window = function(train_start, test_start, min_window) {
      if (is.null(min_window)) {
        return(list("start" = train_start, "end" = test_start - 1, "window" = NULL))
      }

      largest_t <- test_start - train_start
      if (min_window > largest_t) {
        stop("Train min_window must be smaller or equal to largest possible train-set")
      }
      if (min_window < 1) {
        stop("Train min_window must be 1 or larger")
      }

      # Use private$ to access other private methods
      window <- private$sample_window_size(min_window, largest_t)
      start <- private$sample_start_time(test_start, window, train_start, largest_t)

      return(list("start" = start, "window" = window, "end" = start + window))
    },

    sample_window_size = function(min_window, largest_t) {
      window <- 0
      while(min_window > window) {
        window <- base::sample.int(largest_t, 1)
      }
      return(window)
    },

    sample_start_time = function(test_start, window, train_start, largest_t) {
      start <- test_start
      while(start + window >= test_start) {
        start <- base::sample.int(largest_t, 1) + train_start - 1
      }
      return(start)
    },

    process_endogenous_variables = function(simulation_data, models, execution_order){
      endogen_outcomes <- names(models)[!sapply(models, function(x) "ExogenModel" %in% class(x))]
      endogen_execution_order <- execution_order[execution_order %in% endogen_outcomes]
      test_stop <- simulation_data[,max(time), env = list(time = data.table::indices(simulation_data))]

      grp <- data.table::key(simulation_data)
      idx <- data.table::indices(simulation_data)

      for (t in self$test_start:test_stop) {
        for (model in models[endogen_execution_order]) {
          pred <- model$predict(simulation_data, t)
          if(!is.null(model$censor)){
            pred <- pred[,outcome := private$censorfunc(outcome, lower = model$censor$lower, upper = model$censor$upper), env = list(outcome = model$name)]
          }
          simulation_data[pred, on = c(grp, idx), outcome := i.outcome, env = list(outcome = model$name, i.outcome = paste0("i.", model$name))]
        }
      }
      return(simulation_data)
    },

    process_exogenous_variables = function(simulation_data, models, execution_order){
      exogen_outcomes <- names(models)[sapply(models, function(x) "ExogenModel" %in% class(x))]
      exogen_execution_order <- execution_order[execution_order %in% exogen_outcomes]

      grp <- data.table::key(self$data)
      idx <- data.table::indices(self$data)

      for (model in models[exogen_execution_order]) {
        exogen_input <- model$predict()
        n <- names(exogen_input)
        simulation_data[exogen_input, (n):=mget(paste0("i.", n)), on = c(grp, idx)]
        data.table::setkeyv(simulation_data, grp)
        data.table::setindexv(simulation_data, idx)
      }
      return(simulation_data)
    },

    # Helper methods
    clone_models = function() {
      lapply(self$models, function(m) m$clone(deep = TRUE))
    },

    censorfunc = function(x, lower, upper){
      ifelse(x < lower, lower, ifelse(x > upper, upper, x))
    },

    prepare_simulation_data = function(horizon) {
      idx <- data.table::indices(self$data)
      grp <- data.table::key(self$data)

      # Create initial subset
      simulation_data <- self$data[time < self$test_start, env = list(time = idx)]
      # Create complete sequence of times
      all_times <- base::seq(self$train_start, self$test_start + horizon, by = 1)
      # Get unique groups
      all_groups <- unique(simulation_data[[grp]])
      # Create expanded grid
      expanded <- data.table::CJ(
        time = all_times,
        group = all_groups,
        sorted = FALSE
      )

      # Rename columns to match original
      data.table::setnames(expanded, c("time", "group"), c(idx, grp))
      # Merge with original data
      result <- merge(expanded, simulation_data, by = c(grp, idx), all.x = TRUE)

      # Set the keys
      data.table::setkeyv(result, grp)
      data.table::setindexv(result, idx)

      return(result)
    },

    # Debug helper methods
    summarize_data = function() {
      list(
        n_observations = nrow(self$data),
        n_groups = length(unique(self$data[[self$groupvar]])),
        time_range = base::range(self$data[[self$timevar]]),
        variables = names(self$data)
      )
    },

    summarize_models = function() {
      lapply(self$models, function(m) {
        list(
          class = class(m),
          formula = base::deparse(m$formula),
          name = m$name
        )
      })
    },

    summarize_dependency_graph = function() {
      list(
        n_vertices = igraph::vcount(self$dependency_graph),
        n_edges = igraph::ecount(self$dependency_graph),
        outputs = self$outputs
      )
    }
  )
)
