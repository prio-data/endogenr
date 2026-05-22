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
#' @keywords internal
#'
#' @examples
func_in_term <- function(formula, func = "lag") {
  res <- formula |>
    terms() |>
    attr("variables") |>
    sapply(\(x) inherits(x, "call") && (func %in% (x |> all.vars(functions = TRUE))))
  res[-1] # first term is always just list()
}

#' Extract dependency edges from a single formula
#'
#' Iterates over formula terms individually, correctly handling variables that
#' appear both lagged and unlagged (e.g., `y ~ lag(x) + x` produces both
#' `lag_x -> y` and `x -> y` edges).
#'
#' @param formula An R formula.
#' @param outcome Character. The outcome variable name.
#'
#' @return A data.frame with columns `in.` and `out`.
#' @keywords internal
.edges_from_formula <- function(formula, outcome) {
  term_labels <- attr(stats::terms(formula), "term.labels")
  edges <- vector("list", length(term_labels))
  for (i in seq_along(term_labels)) {
    expr <- parse(text = term_labels[[i]])[[1]]
    var_names <- all.vars(expr, functions = FALSE)
    func_names <- setdiff(all.vars(expr, functions = TRUE), var_names)
    has_lag <- "lag" %in% func_names
    prefix <- if (has_lag) "lag_" else ""
    edges[[i]] <- data.frame(in. = paste0(prefix, var_names), out = outcome)
  }
  unique(do.call(rbind, edges))
}

#' Parse a model formula to correctly anticipate the input and output variables
#'
#' @param model
#'
#' @return
#' @keywords internal
#'
#' @examples
parse_formula <- function(model){
  independent_models <- get_independent_models()
  formula <- model$formula

  if (!any(independent_models %in% class(model))) {
    outcome <- base::all.vars(rlang::f_lhs(formula))
    edges <- .edges_from_formula(formula, outcome)

    # Also include variance formula dependencies for heterolm models
    if ("heterolm" %in% class(model) && !is.null(model$variance_formula)) {
      var_edges <- .edges_from_formula(model$variance_formula, outcome)
      if (nrow(var_edges) > 0) {
        edges <- unique(rbind(edges, var_edges))
      }
    }

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
#' @keywords internal
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


#' Compute the maximum lag depth from a formula's AST
#'
#' Walks the formula AST to extract the `n` argument from `lag(expr, n)` calls
#' (defaulting to 1 if omitted) and `k` from rolling functions like
#' `zoo::rollsumr(..., k = k)`. Returns the maximum value found, or 0 if no
#' lag or rolling calls exist.
#'
#' @param formula An R formula.
#' @return Integer. The maximum history depth required.
#' @keywords internal
.max_lag_depth <- function(formula) {
  depths <- integer(0)

  walk <- function(expr) {
    if (!is.call(expr)) return()
    fn <- expr[[1]]
    # Handle namespaced calls like zoo::rollsumr
    if (is.call(fn) && identical(fn[[1]], as.symbol("::"))) {
      fname <- as.character(fn[[3]])
    } else {
      fname <- as.character(fn)
    }

    if (fname == "lag") {
      # lag(expr) or lag(expr, n) or lag(expr, n = 2)
      if (length(expr) >= 3) {
        n_arg <- expr[[3]]
        if (is.numeric(n_arg)) {
          depths[length(depths) + 1L] <<- as.integer(n_arg)
        }
      } else {
        depths[length(depths) + 1L] <<- 1L
      }
    }

    if (fname %in% c("rollsumr", "rollmeanr", "rollapplyr", "rollmaxr",
                      "rollsum", "rollmean", "rollapply", "rollmax")) {
      # zoo rolling functions: rollsumr(x, k, ...) — k is second arg
      if (length(expr) >= 3) {
        k_arg <- expr[[3]]
        if (is.numeric(k_arg)) {
          depths[length(depths) + 1L] <<- as.integer(k_arg)
        }
      }
    }

    # Recurse into sub-expressions
    for (j in seq_along(expr)[-1]) {
      if (is.call(expr[[j]])) walk(expr[[j]])
    }
  }

  rhs <- rlang::f_rhs(formula)
  walk(rhs)
  if (length(depths) == 0L) 0L else max(depths)
}

#' Gives you the independent model types
#'
#' @return
#' @keywords internal
#'
#' @examples
get_independent_models <- function(){
  c("exogen", "parametric_distribution", "univariate_fable")
}
