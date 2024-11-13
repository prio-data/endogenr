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
  second_edges <- igraph::E(dependency_graph)[.to(igraph::V(dependency_graph)[outputs])]
  second <- dependency_graph |> igraph::delete_edges(second_edges)

  # Use topological sort to determine calculation order
  second_order <- igraph::topo_sort(second, mode = "out")
  execution_order <- igraph::V(dependency_graph)$name[second_order]
  return(execution_order[execution_order %in% outputs])
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
