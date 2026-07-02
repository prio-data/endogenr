# Tests for R/systemgraph.R ------------------------------------------------

# ── get_independent_models ───────────────────────────────────────────────

test_that("get_independent_models returns expected types", {
  im <- get_independent_models()
  expect_true("exogen" %in% im)
  expect_true("parametric_distribution" %in% im)
  expect_true("univariate_fable" %in% im)
  expect_false("linear" %in% im)
  expect_false("deterministic" %in% im)
})

# ── func_in_term ─────────────────────────────────────────────────────────

test_that("func_in_term detects lag in formula terms", {
  f <- y ~ lag(x) + z
  res <- func_in_term(f, func = "lag")
  # result has one entry per variable on LHS + RHS (minus the list() sentinel)
  # y (LHS), lag(x), z  →  FALSE, TRUE, FALSE
  expect_equal(res, c(FALSE, TRUE, FALSE))
})

test_that("func_in_term returns FALSE when func is absent", {
  f <- y ~ x + z
  res <- func_in_term(f, func = "lag")
  expect_true(all(!res))
})

test_that("func_in_term detects nested function calls", {
  f <- y ~ lag(log(x)) + z
  res <- func_in_term(f, func = "lag")
  # y, lag(log(x)), z
  expect_equal(res, c(FALSE, TRUE, FALSE))
})

test_that("func_in_term works with multiple lagged terms", {
  f <- y ~ lag(a) + b + lag(c)
  res <- func_in_term(f, func = "lag")
  # y, lag(a), b, lag(c)
  expect_equal(res, c(FALSE, TRUE, FALSE, TRUE))
})

# ── parse_formula (dependent models) ─────────────────────────────────────

test_that("parse_formula extracts edges for a dependent model", {
  m <- new_endogenmodel(y ~ lag(x) + z)
  class(m) <- c("linear", class(m))

  result <- parse_formula(m)
  expect_equal(result$outcome, "y")
  expect_true("lag_x" %in% result$edges$`in`)
  expect_true("z" %in% result$edges$`in`)
  expect_true(all(result$edges$out == "y"))
})

test_that("parse_formula handles formula with no lags", {
  m <- new_endogenmodel(y ~ a + b)
  class(m) <- c("deterministic", class(m))

  result <- parse_formula(m)
  expect_equal(result$outcome, "y")
  expect_true("a" %in% result$edges$`in`)
  expect_true("b" %in% result$edges$`in`)
  # No lag_ prefixes
  expect_false(any(grepl("^lag_", result$edges$`in`)))
})

test_that("parse_formula handles variable appearing both lagged and unlagged", {
  m <- new_endogenmodel(y ~ lag(x) + x)
  class(m) <- c("linear", class(m))

  result <- parse_formula(m)
  expect_equal(result$outcome, "y")
  # Must produce BOTH edges: lag_x -> y AND x -> y
  expect_true("lag_x" %in% result$edges$in.)
  expect_true("x" %in% result$edges$in.)
  expect_equal(nrow(result$edges), 2)
})

test_that("parse_formula handles three terms with mixed lag/unlagged", {
  m <- new_endogenmodel(y ~ z + lag(x) + x)
  class(m) <- c("linear", class(m))

  result <- parse_formula(m)
  in_vars <- result$edges$in.
  expect_true("z" %in% in_vars)
  expect_true("lag_x" %in% in_vars)
  expect_true("x" %in% in_vars)
  expect_equal(nrow(result$edges), 3)
})

test_that("parse_formula handles I() terms with multiple variables", {
  m <- new_endogenmodel(y ~ I(a * b) + lag(z))
  class(m) <- c("deterministic", class(m))

  result <- parse_formula(m)
  in_vars <- result$edges$in.
  expect_true("a" %in% in_vars)
  expect_true("b" %in% in_vars)
  expect_true("lag_z" %in% in_vars)
})

# ── parse_formula (independent models) ───────────────────────────────────

test_that("parse_formula handles independent (one-sided) model", {
  m <- new_endogenmodel(~ population)
  class(m) <- c("exogen", class(m))

  result <- parse_formula(m)
  expect_equal(result$outcome, "population")
  expect_equal(as.character(result$edges[, "in"]), "lag_population")
  expect_equal(as.character(result$edges[, "out"]), "population")
})

# ── update_dependency_graph ──────────────────────────────────────────────

test_that("update_dependency_graph adds vertices and edges", {
  g <- igraph::make_empty_graph(directed = TRUE)

  m <- new_endogenmodel(y ~ lag(x))
  class(m) <- c("linear", class(m))

  g <- update_dependency_graph(m, g)

  vnames <- igraph::V(g)$name
  expect_true("y" %in% vnames)
  expect_true("lag_x" %in% vnames)
  expect_true(igraph::are_adjacent(g, "lag_x", "y"))
})

test_that("update_dependency_graph does not duplicate existing vertices", {
  g <- igraph::make_empty_graph(directed = TRUE)

  m1 <- new_endogenmodel(y ~ lag(x))
  class(m1) <- c("linear", class(m1))
  g <- update_dependency_graph(m1, g)

  m2 <- new_endogenmodel(z ~ lag(x) + y)
  class(m2) <- c("linear", class(m2))
  g <- update_dependency_graph(m2, g)

  # lag_x and y already existed; should not be duplicated
  vnames <- igraph::V(g)$name
  expect_equal(sum(vnames == "lag_x"), 1)
  expect_equal(sum(vnames == "y"), 1)
  # New vertices added
  expect_true("z" %in% vnames)
})

# ── get_execution_order ──────────────────────────────────────────────────

test_that("get_execution_order returns correct topological order", {
  g <- igraph::make_empty_graph(directed = TRUE)

  # y depends on lag(x), z depends on y
  m1 <- new_endogenmodel(y ~ lag(x))
  class(m1) <- c("linear", class(m1))
  g <- update_dependency_graph(m1, g)

  m2 <- new_endogenmodel(z ~ y)
  class(m2) <- c("linear", class(m2))
  g <- update_dependency_graph(m2, g)

  order <- get_execution_order(g)

  # y must come before z in execution order
  expect_true(which(order == "y") < which(order == "z"))
})

test_that("get_execution_order places exogenous variables first", {
  g <- igraph::make_empty_graph(directed = TRUE)

  # exogenous model for population
  m1 <- new_endogenmodel(~ population)
  class(m1) <- c("exogen", class(m1))
  g <- update_dependency_graph(m1, g)

  # dependent model using population
  m2 <- new_endogenmodel(gdp ~ population)
  class(m2) <- c("linear", class(m2))
  g <- update_dependency_graph(m2, g)

  order <- get_execution_order(g)
  expect_true(which(order == "population") < which(order == "gdp"))
})

# ── .classify_term_vars (per-variable lag context) ───────────────────────

test_that(".classify_term_vars puts a bare variable in plain", {
  res <- .classify_term_vars(quote(x))
  expect_equal(res$plain, "x")
  expect_length(res$lagged, 0)
})

test_that(".classify_term_vars puts a lag()-wrapped variable in lagged", {
  res <- .classify_term_vars(quote(lag(x)))
  expect_equal(res$lagged, "x")
  expect_length(res$plain, 0)
})

test_that(".classify_term_vars classifies per-variable inside one I() term", {
  # A same-period driver sharing a term with a lagged self-reference must NOT
  # be swept into the lag bucket — this is the deterministic-model shape.
  res <- .classify_term_vars(quote(I(abs(lag(gdppc) * (1 + gdppc_grwt)))))
  expect_equal(res$lagged, "gdppc")
  expect_equal(res$plain, "gdppc_grwt")
})

test_that(".classify_term_vars reports a both-ways variable in both buckets", {
  res <- .classify_term_vars(quote(lag(x) - x))
  expect_equal(res$lagged, "x")
  expect_equal(res$plain, "x")
})

test_that(".classify_term_vars resolves namespaced lag and ignores inner fns", {
  res <- .classify_term_vars(quote(lag(zoo::rollsumr(w, 5))))
  expect_equal(res$lagged, "w")
  expect_length(res$plain, 0)
})

# ── .edges_from_formula (per-variable edges) ─────────────────────────────

test_that(".edges_from_formula keeps a same-period driver in a mixed I() term", {
  # Regression: I(abs(lag(gdppc)*(1+gdppc_grwt))) must yield gdppc_grwt -> gdppc
  # (same-period), NOT lag_gdppc_grwt -> gdppc. The latter is stripped at
  # execution-order time and silently drops the ordering constraint.
  e <- .edges_from_formula(gdppc ~ I(abs(lag(gdppc) * (1 + gdppc_grwt))), "gdppc")
  expect_true("gdppc_grwt" %in% e[["in."]])      # same-period
  expect_true("lag_gdppc" %in% e[["in."]])        # prior-period self-reference
  expect_false("lag_gdppc_grwt" %in% e[["in."]])  # the bug
  expect_equal(nrow(e), 2L)
})

test_that(".edges_from_formula emits both edges for a both-ways variable", {
  e <- .edges_from_formula(y ~ I(lag(x) - x), "y")
  expect_setequal(e[["in."]], c("lag_x", "x"))
  expect_true(all(e$out == "y"))
})

test_that(".edges_from_formula never emits a spurious empty lag_ vertex", {
  # paste0('lag_', character(0)) returns 'lag_' (zero-length recycled to ''),
  # which would inject a bogus vertex. The guard must prevent it.
  for (f in list(y ~ a + b, y ~ lag(x) + z, y ~ I(a * b) + lag(z))) {
    ins <- .edges_from_formula(f, "y")[["in."]]
    expect_false("lag_" %in% ins)
    expect_false("" %in% ins)
  }
})

test_that(".edges_from_formula returns a zero-row frame for an intercept-only formula", {
  e <- .edges_from_formula(~ 1, "y")
  expect_equal(nrow(e), 0L)
  expect_identical(names(e), c("in.", "out"))
})

# ── parse_formula heterolm variance robustness ───────────────────────────

test_that("parse_formula folds in heterolm variance-formula dependencies", {
  m <- new_endogenmodel(y ~ lag(x))
  m$variance_formula <- ~ lag(z) + w
  class(m) <- c("heterolm", class(m))

  in_vars <- parse_formula(m)$edges[["in."]]
  expect_true("lag_x" %in% in_vars)  # mean
  expect_true("lag_z" %in% in_vars)  # variance, lagged
  expect_true("w" %in% in_vars)      # variance, same-period
})

test_that("parse_formula tolerates an intercept-only heterolm variance", {
  m <- new_endogenmodel(y ~ lag(x))
  m$variance_formula <- ~ 1
  class(m) <- c("heterolm", class(m))

  expect_no_error(res <- parse_formula(m))
  expect_equal(res$outcome, "y")
  expect_true("lag_x" %in% res$edges[["in."]])
})

# ── get_execution_order: same-period ordering & cycle detection ───────────

test_that("get_execution_order orders a same-period producer before its consumer", {
  # deterministic gdppc reads same-period gdppc_grwt (produced by the linear
  # model); gdppc_grwt must be computed first.
  g <- igraph::make_empty_graph(directed = TRUE)

  det <- new_endogenmodel(gdppc ~ I(abs(lag(gdppc) * (1 + gdppc_grwt))))
  class(det) <- c("deterministic", class(det))
  g <- update_dependency_graph(det, g)

  lin <- new_endogenmodel(gdppc_grwt ~ lag(yjbest) + lag(log(gdppc)))
  class(lin) <- c("linear", class(lin))
  g <- update_dependency_graph(lin, g)

  order <- get_execution_order(g)
  expect_true(which(order == "gdppc_grwt") < which(order == "gdppc"))
})

test_that("get_execution_order errors on an unlagged self-reference", {
  g <- igraph::make_empty_graph(directed = TRUE)
  m <- new_endogenmodel(gdppc ~ I(gdppc * 1.02))
  class(m) <- c("deterministic", class(m))
  g <- update_dependency_graph(m, g)

  expect_error(get_execution_order(g), "same-period dependency cycle")
})

test_that("get_execution_order errors on a two-model same-period cycle", {
  g <- igraph::make_empty_graph(directed = TRUE)
  m1 <- new_endogenmodel(a ~ b); class(m1) <- c("linear", class(m1))
  m2 <- new_endogenmodel(b ~ a); class(m2) <- c("linear", class(m2))
  g <- update_dependency_graph(m1, g)
  g <- update_dependency_graph(m2, g)

  err <- expect_error(get_execution_order(g), "same-period dependency cycle")
  expect_match(conditionMessage(err), "a")
  expect_match(conditionMessage(err), "b")
})

test_that("get_execution_order allows a lagged self-reference (no cycle)", {
  g <- igraph::make_empty_graph(directed = TRUE)
  m <- new_endogenmodel(gdppc ~ I(lag(gdppc) * 1.02))
  class(m) <- c("deterministic", class(m))
  g <- update_dependency_graph(m, g)

  expect_no_error(order <- get_execution_order(g))
  expect_true("gdppc" %in% order)
})
# ── lag_ prefix anchoring (execution order) ───────────────────────────────

test_that("variables merely containing 'lag_' are not stripped as lag vertices", {
  # `flag_war` contains the substring "lag_": an unanchored match would strip
  # it from the graph and silently drop its model from the execution order.
  g <- igraph::make_empty_graph(directed = TRUE)
  m1 <- new_endogenmodel(flag_war ~ lag(x))
  class(m1) <- c("linear", class(m1))
  m2 <- new_endogenmodel(x ~ lag(flag_war))
  class(m2) <- c("linear", class(m2))
  g <- update_dependency_graph(m1, g)
  g <- update_dependency_graph(m2, g)

  ord <- get_execution_order(g)
  expect_true("flag_war" %in% ord)
  expect_true("x" %in% ord)
})

test_that("a system with a 'lag_'-substring outcome simulates without NA", {
  set.seed(42)
  dt <- data.table::data.table(
    u = rep(c("a", "b"), each = 20), t = rep(1:20, 2),
    flag_war = stats::rnorm(40), x = stats::rnorm(40)
  )
  models <- list(
    build_model("linear", formula = flag_war ~ lag(x)),
    build_model("linear", formula = x ~ lag(flag_war))
  )
  sys <- setup_system(models, dt, train_start = 1, test_start = 16,
                      horizon = 5, groupvar = "u", timevar = "t",
                      inner_sims = 2)
  expect_setequal(sys$execution_order, c("flag_war", "x"))
  fit <- fit_system(sys, nsim = 2)
  res <- simulate_system(fit)
  expect_false(anyNA(res$flag_war))
  expect_false(anyNA(res$x))
})

# ── shift() lag-context classification ────────────────────────────────────

test_that(".classify_term_vars treats backward shift() like lag()", {
  v <- .classify_term_vars(quote(shift(x)))
  expect_equal(v$lagged, "x")
  expect_equal(v$plain, character(0))

  v2 <- .classify_term_vars(quote(shift(x, 2L, type = "lag")))
  expect_equal(v2$lagged, "x")

  # A forward shift reads the future: stays same-period (cannot be sequenced).
  v3 <- .classify_term_vars(quote(shift(x, type = "lead")))
  expect_equal(v3$plain, "x")
  expect_equal(v3$lagged, character(0))
})

test_that("mutually shift()-lagged models get a valid execution order", {
  g <- igraph::make_empty_graph(directed = TRUE)
  m1 <- new_endogenmodel(y ~ shift(x))
  class(m1) <- c("linear", class(m1))
  m2 <- new_endogenmodel(x ~ shift(y))
  class(m2) <- c("linear", class(m2))
  g <- update_dependency_graph(m1, g)
  g <- update_dependency_graph(m2, g)

  expect_no_error(ord <- get_execution_order(g))
  expect_setequal(ord, c("y", "x"))
})
