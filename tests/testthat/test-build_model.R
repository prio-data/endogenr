# Tests for build_model() and new_endogenmodel() --------------------------

test_that("new_endogenmodel creates a list with class endogenmodel", {
  m <- new_endogenmodel(y ~ x)
  expect_s3_class(m, "endogenmodel")
  expect_equal(m$formula, y ~ x)
})

test_that("build_model returns a function with the model type as class", {
  types <- c("deterministic", "linear", "exogen")
  formulas <- list(
    deterministic = gdp ~ I(abs(gdppc * population)),
    linear        = y ~ lag(x),
    exogen        = ~ population
  )

  for (type in types) {
    f <- build_model(type, formula = formulas[[type]])
    expect_true(is.function(f), info = paste("type:", type))
    expect_true(type %in% class(f), info = paste("type:", type))
  }
})

test_that("build_model errors on unknown model type", {
  expect_error(
    build_model("nonexistent", formula = y ~ x),
    "Unknown model type"
  )
})

test_that("build_model('linear') works with and without boot", {
  f1 <- build_model("linear", formula = y ~ lag(x), boot = "resid")
  expect_true(is.function(f1))
  expect_true("linear" %in% class(f1))

  f2 <- build_model("linear", formula = y ~ lag(x), boot = "wild")
  expect_true("linear" %in% class(f2))

  f3 <- build_model("linear", formula = y ~ lag(x))
  expect_true("linear" %in% class(f3))
})

test_that("build_model('exogen') has correct class", {
  f <- build_model("exogen", formula = ~ population)
  expect_true("exogen" %in% class(f))
})

test_that("build_model('univariate_fable') has correct class", {
  f <- build_model(
    "univariate_fable",
    formula = dem ~ error("A") + trend("N") + season("N"),
    method = "ets"
  )
  expect_true("univariate_fable" %in% class(f))
})

test_that("build_model('parametric_distribution') has correct class", {
  f <- build_model(
    "parametric_distribution",
    formula = ~ gdppc_grwt,
    distribution = "norm",
    start = list(mean = 0, sd = 1)
  )
  expect_true("parametric_distribution" %in% class(f))
})

test_that("build_model preserves original function class alongside type", {
  f <- build_model("deterministic", formula = y ~ I(x))
  expect_true("deterministic" %in% class(f))
  expect_true("purrr_function_partial" %in% class(f))
})
