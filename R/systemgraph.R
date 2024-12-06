#' Get the order of calculation from the dependency graph of a model system
#'
#' @param dependency_graph
#'
#' @return
#' @export
#'
#' @examples
get_execution_order = function(dependency_graph) {
  outputs <- igraph::V(dependency_graph)$name

  lagged_vars <- grepl("lag_", igraph::V(dependency_graph) |> names())

  second <- dependency_graph |> igraph::delete_vertices(lagged_vars)

  isolated_vertices <- which(igraph::degree(second) == 0)
  exogenous_variables <- igraph::V(second)$name[isolated_vertices]

  third <- second |> igraph::delete_vertices(exogenous_variables)

  # Use topological sort to determine calculation order
  endogenous_order <- igraph::topo_sort(third, mode = "out")
  endogenous_order <- igraph::V(third)$name[endogenous_order]
  execution_order <- c(exogenous_variables, endogenous_order)
  return(execution_order)
}

#' Test if a certain function is in each term in a formula
#'
#' @param formula
#' @param func A string
#'
#' @return
#' @export
#'
#' @examples
func_in_term <- function(formula, func = "lag") {
  res <- formula |>
    terms() |>
    attr("variables") |>
    sapply(\(x) inherits(x, "call") && (func %in% (x |> all.vars(functions = TRUE))))
  res[-1] # first term is always just list()
}

#' Parse a model formula to correctly anticipate the input and output variables
#'
#' @param model
#'
#' @return
#' @export
#'
#' @examples
parse_formula <- function(model){
  independent_models <- get_independent_models()
  formula <- model$formula
  outcome <- model$outcome

  if (!any(independent_models %in% class(model))) {
    terms <- base::all.vars(rlang::f_rhs(formula))
    lags <- func_in_term(formula, func = "lag")[-1] # drop lhs
    terms <- ifelse(lags, paste0("lag_", terms), terms)
    outcome <- base::all.vars(rlang::f_lhs(formula))
    edges <- data.frame("in" = terms, "out" = outcome)
    vertices <- unique(unlist(edges))
  }

  if (any(independent_models %in% class(model))) {
    outcome <- base::all.vars(formula)
    input <- paste0("lag_", outcome)
    vertices <- c(input, outcome)
    edges <- cbind("in" = input, "out" = outcome)
  }

  return(list("edges" = edges, "vertices" = vertices, "outcome" = outcome))
}

#' Updates the dependency graph based on a new model
#'
#' @param model
#' @param dependency_graph
#'
#' @return
#' @export
#'
#' @examples
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


#' Gives you the independent model types
#'
#' @return
#' @export
#'
#' @examples
get_independent_models <- function(){
  c("exogen", "parametric_distribution", "univariate_fable")
}
