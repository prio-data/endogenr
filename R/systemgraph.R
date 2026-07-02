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
    # Models that carry a prebuilt smoother/bar-free graph_formula (glmmTMB,
    # gamlss) use it directly so that stats::terms() / .edges_from_formula
    # never sees RE bars, random() calls, or other non-standard syntax.
    if (!is.null(model$graph_formula)) {
      gformula <- model$graph_formula
      outcome  <- base::all.vars(rlang::f_lhs(gformula))
      edges    <- .edges_from_formula(gformula, outcome)
      vertices <- unique(unlist(edges))
      return(list(edges = edges, vertices = vertices, outcome = outcome))
    }

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

# Time-series function name sets used by .required_history(). These are the
# canonical registries defined in R/panel_transform.R and imported here —
# a single source of truth shared by both the materialisation walk
# (panel_materialize) and the history-depth / edge-extraction walk below.
# .pt_cum_fns  — cumulative functions requiring the full prior series (Inf depth)
# .pt_roll_fns — rolling/window functions requiring a finite window

# Fetch a call argument by name (preferred) or by positional index among the
# unnamed arguments. `pos = NULL` restricts the lookup to a named match.
.rh_get <- function(call, name, pos) {
  args <- as.list(call)[-1L]
  if (length(args) == 0L) return(NULL)
  nms <- names(args)
  if (!is.null(nms) && any(nzchar(nms)) && name %in% nms[nzchar(nms)]) {
    return(args[[which(nms == name)[1L]]])
  }
  if (is.null(pos)) return(NULL)
  unnamed <- if (is.null(nms)) seq_along(args) else which(!nzchar(nms))
  if (length(unnamed) >= pos) return(args[[unnamed[pos]]])
  NULL
}

# Resolve a count/window argument to a numeric. NULL -> `default`; a numeric
# literal or evaluable constant -> its value; anything undeterminable -> Inf
# (so the caller falls back to the full series, never under-reads).
.rh_count <- function(x, default) {
  if (is.null(x)) return(default)
  if (is.numeric(x) && length(x) == 1L) return(as.numeric(x))
  v <- tryCatch(eval(x, baseenv()), error = function(e) NULL)
  if (is.numeric(v) && length(v) == 1L) return(as.numeric(v))
  Inf
}

# Backward reach of a rolling window given the function name and an optional
# explicit `align` argument (zoo `*r` = right, `*l` = left, bare = center;
# data.table `froll*` = right by default).
.rh_alignment <- function(fname, align_arg, is_f) {
  if (is.character(align_arg) && length(align_arg) == 1L && nzchar(align_arg)) {
    a <- substr(align_arg, 1L, 1L)
    if (a == "r") return("right")
    if (a == "l") return("left")
    if (a == "c") return("center")
  }
  if (is_f) return("right")
  if (endsWith(fname, "r")) return("right")
  if (endsWith(fname, "l")) return("left")
  "center"
}

# Max required history over a non-time-series call's arguments (head dropped).
.rh_max_args <- function(expr) {
  parts <- as.list(expr)[-1L]
  if (length(parts) == 0L) return(0)
  max(vapply(parts, .req_hist_expr, numeric(1)))
}

# Required history (rows strictly before t) for one language object. See
# .required_history() for the composition rules.
.req_hist_expr <- function(expr) {
  if (!is.call(expr)) return(0)                 # symbol or literal
  fname <- .pt_call_name(expr)
  if (is.na(fname)) return(.rh_max_args(expr))

  if (fname %in% .pt_cum_fns) return(Inf)

  if (fname %in% c("lag", "shift")) {
    inner <- expr[[2L]]
    if (fname == "shift") {
      type <- .rh_get(expr, "type", NULL)
      if (is.character(type) && identical(as.character(type), "lead")) {
        return(.req_hist_expr(inner))
      }
    }
    n <- .rh_count(.rh_get(expr, "n", 2L), default = 1)
    return(n + .req_hist_expr(inner))
  }

  if (fname %in% c("lead", "lead_horizon")) {
    return(.req_hist_expr(expr[[2L]]))
  }

  if (fname == "diff") {
    inner <- expr[[2L]]
    lag_n <- .rh_count(.rh_get(expr, "lag", 2L), default = 1)
    dif_n <- .rh_count(.rh_get(expr, "differences", 3L), default = 1)
    return(lag_n * dif_n + .req_hist_expr(inner))
  }

  if (fname %in% .pt_roll_fns) {
    inner <- expr[[2L]]
    is_f  <- startsWith(fname, "froll")
    wname <- if (grepl("apply", fname)) "width" else if (is_f) "n" else "k"
    k <- .rh_count(.rh_get(expr, wname, 2L), default = NA_real_)
    if (!is.finite(k)) return(Inf)
    align <- .rh_alignment(fname, .rh_get(expr, "align", NULL), is_f)
    back <- switch(align, right = k - 1, center = floor(k / 2), left = 0)
    return(back + .req_hist_expr(inner))
  }

  .rh_max_args(expr)
}

#' Compute the history depth a formula's RHS needs before time `t`
#'
#' Walks the right-hand-side AST and composes nested time-series depths so that
#' materialising the RHS over the `need` rows ending at `t` reproduces the value
#' it would take over the full per-unit series. Unlike [.max_lag_depth()] (which
#' takes the maximum of individual depths), nested transforms are summed, so
#' `lag(lag(x))` needs 2 and `lag(rollmeanr(x, 5), 2)` needs 6.
#'
#' Composition rules:
#' \itemize{
#'   \item `lag(expr, n)` / `shift(expr, n)` -> `n + required(expr)` (default
#'     `n = 1`; `shift(type = "lead")` is forward, contributing only
#'     `required(expr)`).
#'   \item `lead(expr, n)` / `lead_horizon(expr, h)` -> `required(expr)`.
#'   \item `diff(expr, lag, differences)` -> `lag * differences + required(expr)`
#'     (defaults `1, 1`).
#'   \item rolling (`roll*`, `froll*`) -> `required(expr)` plus the window's
#'     backward reach (right: `k - 1`, center: `floor(k / 2)`, left: `0`).
#'     Alignment defaults follow the function name and an explicit `align`
#'     argument overrides it.
#'   \item `cumsum`/`cumprod`/`cummax`/`cummin` and the `*_since_event` helpers
#'     depend on the entire prior series -> `Inf`.
#'   \item any other call -> max over its arguments; symbols/literals -> `0`.
#' }
#'
#' An undeterminable lag count or rolling window yields `Inf` (use the full
#' series), which never under-reads.
#'
#' @param formula An R formula (its RHS is walked) or a language object.
#' @return A non-negative numeric scalar, possibly `Inf`.
#' @keywords internal
.required_history <- function(formula) {
  expr <- if (inherits(formula, "formula")) rlang::f_rhs(formula) else formula
  .req_hist_expr(expr)
}

# Subset `data` to the per-unit history window ending at time `t` that
# .required_history() reports the RHS needs. A finite `need` keeps the rows in
# `(t - need, t]` (floored so at least the `(t - 1, t]` window survives); an
# infinite `need` keeps the entire per-unit history up to `t`.
.history_subset <- function(data, idx, t, need) {
  # Build the row index outside the data.table scope so a time column literally
  # named like an argument (e.g. `t`) cannot shadow it.
  tcol <- data[[idx]]
  keep <- if (is.finite(need)) {
    need <- max(need, 1)
    tcol <= t & tcol > (t - need - 1)
  } else {
    tcol <= t
  }
  data[keep]
}

#' Gives you the independent model types
#'
#' @return Character vector of model-type names that are independent of the
#'   dynamic execution loop (their outputs are populated up front).
#' @keywords internal
get_independent_models <- function(){
  c("exogen", "parametric_distribution", "univariate_fable")
}
