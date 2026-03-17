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
