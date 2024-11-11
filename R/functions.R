get_independent_models <- function(){
  c("exogen", "parametric_distribution")
}

new_endogenmodel <- function(formula){
  structure(
    list(
      formula = formula
    ),
    class = "endogenmodel"
  )
}

create_panel_frame <- function(formula, data) {
  get_naive_formula <- function(data) {
    # Get key and index variables from tsibble
    grp <- tsibble::key_vars(data)
    idx <- tsibble::index_var(data)

    # Get all variable names excluding key and index
    varnames <- names(data)
    varnames <- varnames[!varnames %in% c(grp, idx)]

    # Create new formula
    new_formula <- paste(varnames[1], "~", paste(varnames[2:length(varnames)], collapse = "+")) |>
      stats::formula()

    return(new_formula)
  }

  # Get key and index variables from tsibble
  grp <- tsibble::key_vars(data)
  idx <- tsibble::index_var(data)

  # Update formula to include index variable
  formula <- stats::update(formula, paste(c(". ~ .", idx), collapse = "+"))

  # Create model frame by group
  result <- data |>
    dplyr::group_by(!!!rlang::syms(grp)) |>
    dplyr::group_modify(~ {
      mf <- stats::model.frame(formula, data = ., na.action = NULL)
      tibble::as_tibble(mf)
    }) |>
    dplyr::ungroup() |>
    # Convert back to tsibble
    tsibble::as_tsibble(
      key = !!rlang::sym(grp),
      index = !!rlang::sym(idx)
    )

  # Clean names
  result <- result |>
    janitor::clean_names() |>
    # Ensure tsibble structure is maintained
    tsibble::as_tsibble(
      key = !!rlang::sym(grp),
      index = !!rlang::sym(idx)
    )

  # Get naive formula for the cleaned result
  naive_formula <- get_naive_formula(result)

  return(list(
    "data" = result,
    "naive_formula" = naive_formula
  ))
}

bootstraplm <- function(formula, data, type){
  data <- na.omit(data)
  fitted <- stats::lm(formula, data)
  resid <- residuals(fitted)

  if (type == "resid") {
    resampled_residuals <- base::sample(resid, size = length(resid), replace = TRUE)
  } else if (type == "wild") {
    resampled_residuals <- resid * stats::rnorm(length(resid))
  } else {
    stop("Unknown bootstrap type")
  }

  outcome <- fitted$terms |> rlang::f_lhs() |> as.character()
  data[[outcome]] <- fitted$fitted.values + resampled_residuals
  stats::lm(formula, data)
}



linearmodel <- function(formula = NULL, boot = NULL, data = NULL, ...){
  model <- new_endogenmodel(formula)
  model$boot <- boot
  model$fit_args <- rlang::list2(...)
  model$independent <- FALSE

  panel_frame <- create_panel_frame(model$formula, data)
  model$naive_formula <- panel_frame$naive_formula
  class(model) <- c("linear", class(model))

  if(!is.null(boot)){
    model$fitted <- bootstraplm(model$naive_formula, panel_frame$data, type = boot)
  } else(
    model$fitted = stats::lm(model$naive_formula, panel_frame$data)
  )

  model$coefs <- broom::tidy(model$fitted)
  model$gof <- broom::glance(model$fitted)

  model$outcome <- parse_formula(model)$outcome

  return(model)
}

get_sepi <- function(lmpred){
  #lmpred <- predict(lmfit, se.fit = T)
  se <- lmpred$se.fit
  scale <- lmpred$residual.scale
  sqrt(se^2 + scale^2)
}

select_col_per_row <- function(mat, column_ids){
  cidx <- cbind(1:nrow(mat), column_ids)
  mat[cidx]
}

getpi <- function(lmpred, single_sample = TRUE){
  sepi <- get_sepi(lmpred)
  pi <- lmpred$fit + outer(sepi, stats::rt(100, lmpred$df))

  if(single_sample){
    mycols <- base::sample.int(100, dim(pi)[1], replace = T)
    pi <- select_col_per_row(pi, mycols)
  }
  return(pi)
}

predict.linear <- function(model, data, t, what = "pi", ...) {
  # Get index and key variables
  idx <- tsibble::index_var(data)
  grp <- tsibble::key_vars(data)

  # Find any numbers in the RHS and use that as the max number of past time periods to include.
  #max_history <- stringr::str_extract_all(m$formula |> as.character(), "[0-9]+")[[3]] |> as.numeric() |> max()

  #data <- data |>
  #  dplyr::filter(!!rlang::sym(idx) <= t, !!rlang::sym(idx) > (t-max_history-1))

  # Create panel frame using the tsibble version
  data <- create_panel_frame(model$formula, data)$data

  # Filter for specific time point
  data <- data |>
    dplyr::filter(!!rlang::sym(idx) == t)

  # Make predictions
  if (!is.null(model$boot)) {
    pred <- predict(model$fitted, newdata = data, se.fit = TRUE)
  } else {
    pred <- predict(model$fitted, newdata = data, se.fit = TRUE)
  }

  # Create result tsibble with only necessary columns
  result <- data |>
    dplyr::select(!!!rlang::syms(c(grp, idx, model$outcome))) |>
    tsibble::as_tsibble(
      key = !!rlang::sym(grp),
      index = !!rlang::sym(idx)
    )

  # Update outcome column based on prediction type
  if (what == "expectation") {
    result <- result |>
      dplyr::mutate(!!rlang::sym(model$outcome) := pred$fit)
  } else if (what == "pi") {
    result <- result |>
      dplyr::mutate(
        !!rlang::sym(model$outcome) := getpi(pred)
      )
  } else{
    stop("`what´ must be either `pi´ or `expectation´")
  }

  return(result)
}


deterministicmodel <- function(formula = NULL){
  model <- new_endogenmodel(formula)
  class(model) <- c("deterministic", class(model))
  model$independent <- FALSE

  model$outcome <- parse_formula(model)$outcome
  return(model)
}

exogenmodel <- function(formula = NULL, impute_from = NULL, newdata = NULL){
  fit_exogenous_model <- function(formula, impute_from, data) {
    # Get key and index variables from tsibble
    grp <- tsibble::key_vars(data)
    idx <- tsibble::index_var(data)

    # Filter data based on impute_from
    filtered_data <- data |>
      dplyr::filter(!!rlang::sym(idx) >= impute_from)

    # Update formula to include group and index variables, removing intercept
    model_formula <- stats::update(
      formula,
      paste("~", grp, idx, ". - 1", sep = "+")
    )

    # Create model matrix
    result <- stats::model.matrix(model_formula, filtered_data) |>
      # Convert to tibble and add back key and index columns
      tibble::as_tibble() |>
      # Convert to tsibble
      tsibble::as_tsibble(
        key = !!rlang::sym(grp),
        index = !!rlang::sym(idx)
      )

    # Return a fresh copy
    return(result)
  }

  model <- new_endogenmodel(formula)
  model$impute_from <- impute_from
  model$fitted <- fit_exogenous_model(formula, impute_from, newdata)
  class(model) <- c("exogen", class(model))
  model$independent <- TRUE

  outcome <- parse_formula(model)$outcome
  return(model)
}

predict.exogen <- function(model, ...){
  return(model$fitted)
}

create_distribution_object = function(fitobj) {
  get_nbinom <- function(size, mu){
    prob <- size / (size + mu)
    distributional::dist_negative_binomial(size = size, prob = prob)
  }

  switch(fitobj$distname,
         "norm" = distributional::dist_normal(
           mean = fitobj$estimate["mean"],
           sd = fitobj$estimate["sd"]
         ),
         "cauchy" = distributional::dist_cauchy(
           location = fitobj$estimate["location"],
           scale = fitobj$estimate["scale"]
         ),
         "gumbel" = distributional::dist_gumbel(
           alpha = fitobj$estimate["alpha"],
           scale = fitobj$estimate["scale"]
         ),
         "gamma" = distributional::dist_gamma(
           shape = fitobj$estimate["shape"],
           rate = fitobj$estimate["rate"]
         ),
         "t_ls" = distributional::dist_student_t(
           df = fitobj$estimate["df"],
           mu = fitobj$estimate["mu"],
           sigma = fitobj$estimate["sigma"]
         ),
         "nbinom" = get_nbinom(fit$estimate["size"], fitobj$estimate["mu"]),
         "pois" = distributional::dist_poisson(
           lambda = fitobj$estimate["lambda"]
         ),
         stop("Unsupported distribution")
  )
}

fit_parametric_distribution_model <- function(model, data){
  if (is.null(model$fit_args)) {
    fitted <- fitdistrplus::fitdist(data[[model$outcome]], model$distribution)
  } else {
    args <- model$fit_args
    args$data <- data[[model$outcome]]
    args$distr <- model$distribution
    fitted <- do.call(fitdistrplus::fitdist, args)
  }
  create_distribution_object(fitted)
}

parametric_distribution_model <- function(formula = NULL, distribution = NULL, data = NULL, ...){
  model <- new_endogenmodel(formula)
  model$distribution <- distribution
  model$fit_args <- rlang::list2(...)

  class(model) <- c("parametric_distribution", class(model))
  model$independent <- TRUE
  model$outcome <- parse_formula(model)$outcome
  model$fitted <- fit_parametric_distribution_model(model, data)

  return(model)
}

predict.deterministic <- function(model, t, data, ...) {
  # Get key and index variables from tsibble
  grp <- tsibble::key_vars(data)
  idx <- tsibble::index_var(data)

  # Get the transformed variable name (y_star) and original outcome variable (y)
  y_star <- attr(stats::terms(model$formula)[1], "term.labels")
  y <- model$outcome

  # Check if the formula is properly wrapped in I()
  if (!stringr::str_detect(y_star, "I\\(")) {
    stop("The formula to apply must be the first term and wrapped in I()")
  }

  # Update formula to include index variable
  frm <- stats::update(model$formula, paste(c(". ~ .", idx), collapse = "+"))

  # Create model frame by group
  result <- data %>%
    dplyr::group_by(!!!rlang::syms(grp)) %>%
    dplyr::group_modify(~ {
      mf <- stats::model.frame(frm, data = ., na.action = NULL)
      tibble::as_tibble(mf)
    }) %>%
    dplyr::ungroup() %>%
    # Remove the original outcome column if it exists
    dplyr::select(-!!rlang::sym(y)) %>%
    # Filter for specific time point
    dplyr::filter(!!rlang::sym(idx) == t) %>%
    # Rename the transformed variable to the original outcome name
    dplyr::rename(!!rlang::sym(y) := !!rlang::sym(y_star)) %>%
    # Convert to tsibble
    tsibble::as_tsibble(
      key = !!rlang::sym(grp),
      index = !!rlang::sym(idx)
    )
  return(result)
}

dt_ls <- function(x, df=1, mu=0, sigma=1) 1/sigma * stats::dt((x - mu)/sigma, df)
pt_ls <- function(q, df=1, mu=0, sigma=1) stats::pt((q - mu)/sigma, df)
qt_ls <- function(p, df=1, mu=0, sigma=1) stats::qt(p, df)*sigma + mu
rt_ls <- function(n, df=1, mu=0, sigma=1) stats::rt(n, df)*sigma + mu
#
# predict.parametric_distribution <- function(model, data, test_start, ...){
#   idx <- data.table::indices(data)
#   grp <- data.table::key(data)
#
#   pred_data <- data[time >= test_start, .(group, time, outcome), env = list(time = idx, group = grp, outcome = model$outcome)]
#   pred_data[,outcome := distributional::generate(model$fitted, nrow(pred_data)), env = list(outcome = model$outcome)]
#   data.table::setindexv(pred_data, idx)
#   return(pred_data |> data.table::copy())
# }

predict.parametric_distribution <- function(model, data, test_start, ...) {
  # Get index and key variables from tsibble
  idx <- tsibble::index_var(data)
  grp <- tsibble::key_vars(data)

  # Create prediction data frame with only necessary columns
  pred_data <- data |>
    dplyr::filter(!!rlang::sym(idx) >= test_start) |>
    dplyr::select(
      !!!rlang::syms(grp),
      !!rlang::sym(idx),
      !!rlang::sym(model$outcome)
    )
    # Generate new outcomes using the fitted distribution
  pred_data <- pred_data |>
    dplyr::mutate(
      !!rlang::sym(model$outcome) := distributional::generate(
        model$fitted,
        dplyr::n()
      )$mean
    ) %>%
    # Ensure tsibble structure
    tsibble::as_tsibble(
      key = !!rlang::sym(grp),
      index = !!rlang::sym(idx)
    )

  # Return a fresh copy
  return(pred_data)
}

get_execution_order = function(dependency_graph) {
  outputs <- igraph::V(dependency_graph)$name
  second_edges <- igraph::E(dependency_graph)[.to(igraph::V(dependency_graph)[outputs])]
  second <- dependency_graph |> igraph::delete_edges(second_edges)

  # Use topological sort to determine calculation order
  second_order <- igraph::topo_sort(second, mode = "out")
  execution_order <- igraph::V(dependency_graph)$name[second_order]
  return(execution_order[execution_order %in% outputs])
}

parse_formula <- function(model){
  independent_models <- get_independent_models()
  formula <- model$formula
  outcome <- model$outcome

  if (!any(independent_models %in% class(model))) {
    terms <- base::all.vars(rlang::f_rhs(formula))
    outcome <- base::all.vars(rlang::f_lhs(formula))
    edges <- data.frame("in" = terms, "out" = outcome)
    vertices <- unique(unlist(edges))
  }

  if (any(independent_models %in% class(model))) {
    outcome <- vertices <- base::all.vars(formula)
    edges <- cbind("in" = vertices, "out" = outcome)
  }

  return(list("edges" = edges, "vertices" = vertices, "outcome" = outcome))
}

update_dependency_graph <- function(model, dependency_graph) {
  edges_vertices_outcome <- parse_formula(model)
  edges <- edges_vertices_outcome$edges
  vertices <- edges_vertices_outcome$vertices

  # Add new vertices if they don't exist
  existing_vertices <- igraph::V(dependency_graph)$name
  new_vertices <- base::setdiff(vertices, existing_vertices)

  if (length(new_vertices) > 0) {
    dependency_graph <- igraph::add_vertices(
      dependency_graph,
      length(new_vertices),
      name = new_vertices
    )
  }

  # Add edges
  edge_list <- as.matrix(edges)
  dependency_graph <- igraph::add_edges(
    dependency_graph,
    t(edge_list)
  )

  return(dependency_graph)
}


prepare_simulation_data <- function(data, groupvar, timevar, train_start, test_start, horizon) {
  # Convert groupvar and timevar to symbols for tidy evaluation
  group_sym <- rlang::sym(groupvar)
  time_sym <- rlang::sym(timevar)

  # Filter training data
  simulation_data <- data %>%
    dplyr::filter(!!time_sym < test_start)

  # Create sequence of all times
  all_times <- base::seq(train_start, test_start + horizon, by = 1)

  # Get unique groups
  all_groups <- simulation_data %>%
    dplyr::distinct(!!group_sym) %>%
    dplyr::pull(!!group_sym)

  # Create expanded grid using tidyr
  expanded <- tidyr::expand_grid(
    !!time_sym := all_times,
    !!group_sym := all_groups
  )

  # Convert to tsibble
  result <- expanded %>%
    tsibble::as_tsibble(
      key = !!group_sym,
      index = !!time_sym
    ) %>%
    dplyr::left_join(
      simulation_data,
      by = c(groupvar, timevar)
    )

  return(result)
}

process_independent_models <- function(simulation_data, models, test_start) {
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
    pred <- predict(model, data = simulation_data, test_start = test_start)

    # Update the simulation data with new predictions
    simulation_data <- simulation_data |>
      dplyr::rows_patch(pred, by = c(grp, idx))
  }

  return(simulation_data)
}

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
  for (t in test_start:(test_start + horizon)) {
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

# process_dependent_models = function(simulation_data, models, test_start, horizon, execution_order){
#   dependent_models <- models[sapply(models, function(x) !x$independent)]
#   outcomes <- sapply(dependent_models, function(x) parse_formula(x)$outcome)
#   names(dependent_models) <- outcomes
#   execution_order <- execution_order[execution_order %in% outcomes]
#   dependent_models <- dependent_models[execution_order]
#
#   grp <- data.table::key(simulation_data)
#   idx <- data.table::indices(simulation_data)
#
#   simulation_data <- data.table::copy(simulation_data)
#   simulation_data[["update_index"]] <- 1:nrow(simulation_data)
#
#   for (t in test_start:(test_start + horizon)) {
#     for (model in dependent_models) {
#       pred <- predict(model, t = t, data = simulation_data)
#       # consider adding possible constraints to the outcome here
#       outcome <- parse_formula(model)$outcome
#
#       data.table::set(simulation_data,
#                       i = simulation_data[pred, update_index, on = c(grp, idx)],
#                       j = outcome,
#                       value = pred[[outcome]])
#     }
#   }
#   simulation_data[,update_index := NULL]
#   data.table::setkeyv(simulation_data, grp)
#   data.table::setindexv(simulation_data, idx)
#   return(data.table::copy(simulation_data))
# }
