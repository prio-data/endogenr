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

create_panel_frame = function(formula, data) {
  get_naive_formula = function(data) {
    grp <- data.table::key(data)
    idx <- data.table::indices(data)

    varnames <- names(data)
    varnames <- varnames[!varnames %in% c(grp, idx)]
    new_formula <- paste(varnames[1], "~", paste(varnames[2:length(varnames)], collapse = "+")) |> stats::formula()
    return(new_formula)
  }

  grp <- data.table::key(data)
  idx <- data.table::indices(data)

  formula <- stats::update(formula, paste(c(". ~ .", idx), collapse = "+"))
  result <- data[, as.list(stats::model.frame(formula, data = .SD, na.action = NULL)), by = grp]
  result <- result |> janitor::clean_names()
  data.table::setkeyv(result, grp)
  data.table::setindexv(result, idx)

  naive_formula <- get_naive_formula(result)
  return(list("data" = result, "naive_formula" = naive_formula))
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

predict.linear <- function(model, data, t, what = "expectation", ...){
  data <- create_panel_frame(model$formula, data)$data
  idx <- data.table::indices(data)
  grp <- data.table::key(data)
  data <- data[time == t, env = list(time = idx)]

  if(!is.null(model$boot)){
    pred <- predict(model$fitted, newdata = data, se.fit = TRUE)
  } else{
    pred <- predict(model$fitted, newdata = data, se.fit = TRUE)
  }

  result <- data[, .(group, time, outcome), env = list(group = grp, time = idx, outcome = model$outcome)]

  if(what == "expecation"){
    result[, outcome := pred$fit, env = list(outcome = model$outcome)]
  } else if(what == "pi"){
    result[, outcome := getpi(pred), env = list(outcome = model$outcome)]
  }

  return(data.table::copy(result))
}


deterministicmodel <- function(formula = NULL){
  model <- new_endogenmodel(formula)
  class(model) <- c("deterministic", class(model))
  model$independent <- FALSE

  model$outcome <- parse_formula(model)$outcome
  return(model)
}

exogenmodel <- function(formula = NULL, impute_from = NULL, newdata = NULL){
  fit_exogenous_model = function(formula, impute_from, data) {
    grp <- data.table::key(data)
    idx <- data.table::indices(data)
    data <- data[time >= impute_from, env = list(time = idx)]
    model_formula <- stats::update(formula, paste(c("~ . -1", grp, idx), collapse = " + "))

    result <- stats::model.matrix(model_formula, data) |> data.table::as.data.table()
    data.table::setkeyv(result, grp)
    data.table::setindexv(result, idx)
    return(result |> data.table::copy())
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
    fitted <- fitdistrplus::fitdist(data[[outcome]], model$distribution)
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

predict.deterministic <- function(model, t, data, ...){
  grp <- data.table::key(data)
  idx <- data.table::indices(data)
  y_star <- attr(stats::terms(model$formula)[1], "term.labels")
  y <- model$outcome

  if (!stringr::str_detect(y_star, "I\\(")) {
    stop("The formula to apply must be the first term and wrapped in I()")
  }

  frm <- stats::update(model$formula, paste(c(". ~ .", idx), collapse = "+"))
  result <- data[, as.list(stats::model.frame(frm, data = .SD, na.action = NULL)), by = grp]
  result <- result[, outcome := NULL, env = list(outcome = y)]
  result <- result[time == t, env = list(time = idx)]
  data.table::setnames(result, y_star, y)
  data.table::setkeyv(result, grp)
  data.table::setindexv(result, idx)
  return(result)
}

dt_ls <- function(x, df=1, mu=0, sigma=1) 1/sigma * stats::dt((x - mu)/sigma, df)
pt_ls <- function(q, df=1, mu=0, sigma=1) stats::pt((q - mu)/sigma, df)
qt_ls <- function(p, df=1, mu=0, sigma=1) stats::qt(p, df)*sigma + mu
rt_ls <- function(n, df=1, mu=0, sigma=1) stats::rt(n, df)*sigma + mu

predict.parametric_distribution <- function(model, data, test_start, ...){
  idx <- data.table::indices(data)
  grp <- data.table::key(data)

  pred_data <- data[time >= test_start, .(group, time, outcome), env = list(time = idx, group = grp, outcome = model$outcome)]
  pred_data[,outcome := distributional::generate(model$fitted, nrow(pred_data)), env = list(outcome = model$outcome)]
  data.table::setindexv(pred_data, idx)
  return(pred_data |> data.table::copy())
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
  if(!is.numeric(train_start)){stop("`train_start´ must be numeric")}

  simulation_data <- data[time < test_start, env = list(time = timevar)]

  all_times <- base::seq(train_start, test_start + horizon, by = 1)

  all_groups <- unique(simulation_data[[groupvar]])

  expanded <- data.table::CJ(
    time = all_times,
    group = all_groups,
    sorted = FALSE
  )

  data.table::setnames(expanded, c("time", "group"), c(timevar, groupvar))
  # Merge with original data
  result <- merge(expanded, simulation_data, by = c(groupvar, timevar), all.x = TRUE)

  # Set the keys
  data.table::setkeyv(result, groupvar)
  data.table::setindexv(result, timevar)

  return(result)
}

process_independent_models <- function(simulation_data, models, test_start){
  if(!is.numeric(test_start)){stop("`test_start´ must be numeric")}

  grp <- data.table::key(simulation_data)
  idx <- data.table::indices(simulation_data)

  simulation_data <- data.table::copy(simulation_data)
  simulation_data[["update_index"]] <- 1:nrow(simulation_data)
  independent_models <- models[sapply(models, function(x) x$independent)]

  for (model in independent_models) {
    result <- predict(model, data = simulation_data, test_start = test_start)
    outcome <- parse_formula(model)$outcome

    data.table::set(simulation_data,
                    i = simulation_data[result, update_index, on = c(grp, idx)],
                    j = outcome,
                    value = result[[outcome]])
  }

  simulation_data[,update_index := NULL]
  data.table::setkeyv(simulation_data, grp)
  data.table::setindexv(simulation_data, idx)
  return(data.table::copy(simulation_data))
}

process_dependent_models = function(simulation_data, models, test_start, horizon, execution_order){
  dependent_models <- models[sapply(models, function(x) !x$independent)]
  outcomes <- sapply(dependent_models, function(x) parse_formula(x)$outcome)
  names(dependent_models) <- outcomes
  execution_order <- execution_order[execution_order %in% outcomes]
  dependent_models <- dependent_models[execution_order]

  grp <- data.table::key(simulation_data)
  idx <- data.table::indices(simulation_data)

  simulation_data <- data.table::copy(simulation_data)
  simulation_data[["update_index"]] <- 1:nrow(simulation_data)

  for (t in test_start:(test_start + horizon)) {
    for (model in dependent_models) {
      pred <- predict(model, t = t, data = simulation_data)
      # consider adding possible constraints to the outcome here
      outcome <- parse_formula(model)$outcome

      data.table::set(simulation_data,
                      i = simulation_data[pred, update_index, on = c(grp, idx)],
                      j = outcome,
                      value = pred[[outcome]])
    }
  }
  simulation_data[,update_index := NULL]
  data.table::setkeyv(simulation_data, grp)
  data.table::setindexv(simulation_data, idx)
  return(data.table::copy(simulation_data))
}
