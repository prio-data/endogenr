# Tests for build_model(), fit_model(), and new_endogenmodel() ---------------

test_that("new_endogenmodel creates a list with class endogenmodel", {
  m <- new_endogenmodel(y ~ x)
  expect_s3_class(m, "endogenmodel")
  expect_equal(m$formula, y ~ x)
})

test_that("build_model returns an endogenr_spec with the model type as class", {
  types <- c("deterministic", "linear", "exogen")
  formulas <- list(
    deterministic = gdp ~ I(abs(gdppc * population)),
    linear        = y ~ lag(x),
    exogen        = ~ population
  )

  for (type in types) {
    spec <- build_model(type, formula = formulas[[type]])
    expect_s3_class(spec, "endogenr_spec")
    expect_s3_class(spec, paste0(type, "_spec"))
    expect_equal(spec$type, type)
    expect_equal(spec$formula, formulas[[type]])
  }
})

test_that("build_model errors on unknown model type", {
  expect_error(
    build_model("nonexistent", formula = y ~ x),
    "Unknown model type"
  )
})

test_that("build_model('linear') captures boot in args", {
  s1 <- build_model("linear", formula = y ~ lag(x), boot = "resid")
  expect_equal(s1$args$boot, "resid")

  s2 <- build_model("linear", formula = y ~ lag(x), boot = "wild")
  expect_equal(s2$args$boot, "wild")

  s3 <- build_model("linear", formula = y ~ lag(x))
  expect_null(s3$args$boot)
})

test_that("build_model('exogen') has correct class", {
  spec <- build_model("exogen", formula = ~ population)
  expect_s3_class(spec, "exogen_spec")
})

test_that("build_model('univariate_fable') captures method in args", {
  spec <- build_model(
    "univariate_fable",
    formula = dem ~ error("A") + trend("N") + season("N"),
    method = "ets"
  )
  expect_s3_class(spec, "univariate_fable_spec")
  expect_equal(spec$args$method, "ets")
})

test_that("build_model('parametric_distribution') captures distribution args", {
  spec <- build_model(
    "parametric_distribution",
    formula = ~ gdppc_grwt,
    distribution = "norm",
    start = list(mean = 0, sd = 1)
  )
  expect_s3_class(spec, "parametric_distribution_spec")
  expect_equal(spec$args$distribution, "norm")
  expect_equal(spec$args$start, list(mean = 0, sd = 1))
})

test_that("build_model spec has type, formula, and args fields", {
  spec <- build_model("deterministic", formula = y ~ I(x))
  expect_true(all(c("type", "formula", "args") %in% names(spec)))
  expect_equal(spec$type, "deterministic")
})

test_that("build_model stores and validates bounds", {
  expect_equal(build_model("linear", y ~ lag(x), bounds = c(-1, 1))$bounds, c(-1, 1))
  expect_null(build_model("linear", y ~ lag(x))$bounds)
  expect_error(build_model("linear", y ~ lag(x), bounds = c(1, -1)), "lower <= upper")
  expect_error(build_model("linear", y ~ lag(x), bounds = 1), "length-2")
  # bounds is a named formal after ...; it must not leak into spec$args
  expect_null(build_model("linear", y ~ lag(x), bounds = c(0, 1))$args$bounds)
})

# ── fit_model tests ─────────────────────────────────────────────────────────

test_that("fit_model.linear_spec produces a fitted linear model", {
  dt <- data.table::as.data.table(expand.grid(gwcode = c(1, 2), year = 2000:2010))
  set.seed(1)
  dt$x <- rnorm(nrow(dt))
  dt$y <- 0.5 * dt$x + rnorm(nrow(dt), sd = 0.2)
  data.table::setkeyv(dt, c("gwcode", "year"))

  spec <- build_model("linear", formula = y ~ lag(x))
  ctx <- panel_context(unit = "gwcode", time = "year")

  fitted <- fit_model(spec, data = dt, ctx = ctx)
  expect_s3_class(fitted, "linear")
  expect_s3_class(fitted, "endogenmodel")
  expect_true(!is.null(fitted$fitted))
  expect_true(!is.null(fitted$coefs))
})

test_that("fit_model.deterministic_spec produces a fitted deterministic model", {
  spec <- build_model("deterministic", formula = y ~ I(x * 2))
  ctx <- panel_context(unit = "gwcode", time = "year")

  fitted <- fit_model(spec, ctx = ctx)
  expect_s3_class(fitted, "deterministic")
  expect_equal(fitted$outcome, "y")
})
