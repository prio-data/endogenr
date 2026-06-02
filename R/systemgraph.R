#' Determine the order of calculation from a model-system dependency graph
#'
#' Given the directed dependency graph built by [update_dependency_graph()],
#' returns a character vector ordering the model outcomes so that every
#' dependency is computed before any model that requires it. The order is:
#' lagged inputs are stripped, exogenous variables (isolated after the strip)
#' come first, then the remaining endogenous outcomes are topologically sorted
#' along their `in → out` edges.
#'
#' @param dependency_graph A directed `igraph` graph whose vertices are
#'   variable names (with `lag_` prefix for lagged inputs) and whose edges
#'   point from input to outcome.
#'
#' @return Character vector. Variable names in execution order.
#' @family build
#' @export
#'
#' @examples
#' g <- igraph::make_empty_graph(directed = TRUE) |>
#'   igraph::add_vertices(3, name = c("lag_x", "x", "y")) |>
#'   igraph::add_edges(c("lag_x", "y"))
#' get_execution_order(g)
get_execution_order = function(dependency_graph) {
  outputs <- igraph::V(dependency_graph)$name

  lagged_vars <- grepl("lag_", igraph::V(dependency_graph) |> names())

  second <- dependency_graph |> igraph::delete_vertices(lagged_vars)

  isolated_vertices <- which(igraph::degree(second) == 0)
  exogenous_variables <- igraph::V(second)$name[isolated_vertices]

  third <- second |> igraph::delete_vertices(exogenous_variables)

  # A cycle among the same-period (non-lag) edges means a model depends on its
  # own or another model's current-period value, which cannot be sequenced.
  if (igraph::gorder(third) > 0L && !igraph::is_dag(third)) {
    offenders <- .same_period_cycle_vars(third)
    stop("Model system has a same-period dependency cycle among: ",
         paste(offenders, collapse = ", "),
         ". A model cannot depend on its own or another model's current-period ",
         "value; use lag() to reference a prior period.",
         call. = FALSE)
  }

  # Use topological sort to determine calculation order
  endogenous_order <- igraph::topo_sort(third, mode = "out")
  endogenous_order <- igraph::V(third)$name[endogenous_order]
  execution_order <- c(exogenous_variables, endogenous_order)
  return(execution_order)
}

#' Identify variables participating in a same-period dependency cycle
#'
#' @param g A directed `igraph` graph (the non-lag, endogenous subgraph).
#' @return Character vector of vertex names that sit in a non-trivial strongly
#'   connected component or carry a self-loop.
#' @keywords internal
.same_period_cycle_vars <- function(g) {
  comp <- igraph::components(g, mode = "strong")
  multi <- names(which(table(comp$membership) > 1L))
  cyclic <- igraph::V(g)$name[as.character(comp$membership) %in% multi]
  el <- igraph::as_edgelist(g)
  self_loops <- if (nrow(el) > 0L) el[el[, 1L] == el[, 2L], 1L] else character(0)
  unique(c(cyclic, self_loops))
}

#' Test if a certain function is in each term in a formula
#'
#' @param formula An R formula.
#' @param func Character. Function name to look for inside each term.
#'
#' @return Logical vector, one entry per RHS term, indicating whether the
#'   given function appears anywhere inside that term.
#' @keywords internal
func_in_term <- function(formula, func = "lag") {
  res <- formula |>
    terms() |>
    attr("variables") |>
    sapply(\(x) inherits(x, "call") && (func %in% (x |> all.vars(functions = TRUE))))
  res[-1] # first term is always just list()
}

#' Classify the variables in a single formula term by lag context
#'
#' Walks the term's AST and records each variable as lagged (it appears inside
#' a `lag()` call) or same-period (it does not). The classification is
#' per-variable, not per-term: a same-period variable sharing a term with a
#' lagged one (e.g. `I(lag(y) * x)`) is still same-period. A variable that is
#' both (e.g. `I(x - lag(x))`) is reported in both buckets.
#'
#' @param expr A language object (one parsed formula term).
#'
#' @return A list with `lagged` and `plain` character vectors (each unique).
#' @keywords internal
.classify_term_vars <- function(expr) {
  lagged <- character(0)
  plain  <- character(0)
  walk <- function(e, under_lag) {
    if (is.symbol(e)) {
      nm <- as.character(e)
      if (nzchar(nm)) {
        if (under_lag) lagged[[length(lagged) + 1L]] <<- nm
        else           plain[[length(plain) + 1L]]   <<- nm
      }
      return(invisible())
    }
    if (!is.call(e)) return(invisible())
    # A `lag()` call (bare or namespaced) sets the lag context for its args.
    under_lag <- under_lag || identical(.pt_call_name(e), "lag")
    for (j in seq_along(e)[-1L]) walk(e[[j]], under_lag)
  }
  walk(expr, FALSE)
  list(lagged = unique(lagged), plain = unique(plain))
}

#' Extract dependency edges from a single formula
#'
#' Iterates over formula terms and classifies each variable per-term by lag
#' context (see [.classify_term_vars()]). A variable inside a `lag()` call
#' becomes a `lag_`-prefixed input (a prior-period dependency, stripped when
#' computing execution order); a same-period variable keeps its bare name. A
#' variable that is both lagged and unlagged (across terms as in
#' `y ~ lag(x) + x`, or within one term as in `y ~ I(lag(x) - x)`) produces
#' both `lag_x -> y` and `x -> y` edges.
#'
#' @param formula An R formula.
#' @param outcome Character. The outcome variable name.
#'
#' @return A data.frame with columns `in.` and `out` (zero rows when the
#'   formula has no variable terms, e.g. `~ 1`).
#' @keywords internal
.edges_from_formula <- function(formula, outcome) {
  term_labels <- attr(stats::terms(formula), "term.labels")
  ins <- character(0)
  for (lbl in term_labels) {
    vars <- .classify_term_vars(parse(text = lbl)[[1]])
    if (length(vars$lagged)) ins <- c(ins, paste0("lag_", vars$lagged))
    if (length(vars$plain))  ins <- c(ins, vars$plain)
  }
  ins <- unique(ins)
  data.frame(in. = ins, out = rep(outcome, length(ins)), stringsAsFactors = FALSE)
}
#' Parse a model formula to correctly anticipate the input and output variables
#'
#' @param model An endogenmodel with `$formula` and `class()` set.
#'
#' @return A list with `edges` (data.frame of in/out edges), `vertices`
#'   (character vector), and `outcome` (the LHS variable name).
#' @keywords internal
parse_formula <- function(model){
  independent_models <- get_independent_models()
  formula <- model$formula

  if (!any(independent_models %in% class(model))) {
    outcome <- base::all.vars(rlang::f_lhs(formula))
    edges <- .edges_from_formula(formula, outcome)

    # Also include variance formula dependencies for heterolm models
    if ("heterolm" %in% class(model) && inherits(model$variance_formula, "formula")) {
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
#' @param model An endogenmodel produced by `fit_model()`.
#' @param dependency_graph A directed `igraph` graph (possibly empty).
#'
#' @return The updated `igraph` graph with vertices and edges from `model`
#'   merged into it.
#' @keywords internal
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
#' @return Character vector of model-type names that are independent of the
#'   dynamic execution loop (their outputs are populated up front).
#' @keywords internal
get_independent_models <- function(){
  c("exogen", "parametric_distribution", "univariate_fable")
}
