# Tests for R/validate.R ---------------------------------------------------

# ── Helper ──────────────────────────────────────────────────────────────────

make_panel <- function() {
  dt <- data.table::as.data.table(
    expand.grid(gwcode = c(1, 2), year = 2000:2010)
  )
  set.seed(1)
  dt$x <- rnorm(nrow(dt))
  dt$y <- 0.5 * dt$x + rnorm(nrow(dt), sd = 0.2)
  data.table::setkeyv(dt, c("gwcode", "year"))
  dt
}

make_ctx <- function() {
  panel_context(unit = "gwcode", time = "year")
}

# ── validate_panel ──────────────────────────────────────────────────────────

test_that("validate_panel passes on well-formed data", {
  dt <- make_panel()
  ctx <- make_ctx()
  expect_true(validate_panel(dt, ctx, test_start = 2008))
})

test_that("validate_panel errors on missing columns", {
  dt <- make_panel()
  ctx <- panel_context(unit = "country", time = "year")
  expect_error(validate_panel(dt, ctx, test_start = 2008), "Missing required columns.*country")
})

test_that("validate_panel errors on non-numeric time", {
  dt <- make_panel()
  dt$year <- as.character(dt$year)
  ctx <- make_ctx()
  expect_error(validate_panel(dt, ctx, test_start = 2008), "must be numeric")
})

test_that("validate_panel errors on non-integer time", {
  dt <- make_panel()
  dt$year <- dt$year + 0.5
  ctx <- make_ctx()
  expect_error(validate_panel(dt, ctx, test_start = 2008), "integer values")
})

test_that("validate_panel errors on gap in time series", {
  dt <- make_panel()
  # Remove year 2005 for unit 1

  dt <- dt[!(gwcode == 1 & year == 2005)]
  ctx <- make_ctx()
  expect_error(validate_panel(dt, ctx, test_start = 2008), "Non-contiguous.*gap after time step 2004")
})

test_that("validate_panel accepts late entry (unbalanced panel)", {
  dt <- make_panel()
  # Unit 2 enters in 2003
  dt <- dt[!(gwcode == 2 & year < 2003)]
  ctx <- make_ctx()
  expect_true(validate_panel(dt, ctx, test_start = 2008))
})

test_that("validate_panel accepts early exit (unbalanced panel)", {
  dt <- make_panel()
  # Unit 2 exits after 2005; unit 1 spans the full range
  dt <- dt[!(gwcode == 2 & year > 2005)]
  ctx <- make_ctx()
  expect_true(validate_panel(dt, ctx, test_start = 2008))
})

test_that("validate_panel errors when no unit covers the forecast origin", {
  dt <- make_panel()
  # Every unit ends before test_start - 1 = 2007
  dt <- dt[year <= 2005]
  ctx <- make_ctx()
  expect_error(validate_panel(dt, ctx, test_start = 2008),
               "No units have data at the forecast origin")
})

test_that("validate_panel errors on NA initial state", {
  dt <- make_panel()
  ctx <- make_ctx()
  # Set y to NA at the last pre-forecast time step for unit 1
  dt[gwcode == 1 & year == 2007, y := NA]
  expect_error(
    validate_panel(dt, ctx, test_start = 2008, model_outcomes = "y"),
    "Missing initial state.*2007.*outcome 'y'"
  )
})

test_that("validate_panel errors on unsorted data", {
  dt <- make_panel()
  # Reverse order
  dt <- dt[nrow(dt):1]
  ctx <- make_ctx()
  expect_error(validate_panel(dt, ctx, test_start = 2008), "must be sorted")
})

test_that("validate_panel skips initial state check when model_outcomes is NULL", {
  dt <- make_panel()
  ctx <- make_ctx()
  dt[gwcode == 1 & year == 2007, y := NA]
  # Should pass - no outcomes to check

  expect_true(validate_panel(dt, ctx, test_start = 2008))
})

# ── validate_system_closure ─────────────────────────────────────────────────

test_that("validate_system_closure passes for a closed system", {
  models <- list(
    build_model("exogen", formula = ~x),
    build_model("linear", formula = y ~ lag(x))
  )
  expect_true(validate_system_closure(models, c("gwcode", "year", "x", "y")))
})

test_that("validate_system_closure errors on missing RHS variable", {
  models <- list(
    build_model("linear", formula = y ~ lag(z))
  )
  expect_error(
    validate_system_closure(models, c("gwcode", "year", "x", "y")),
    "missing from the input data.*z"
  )
})

test_that("validate_system_closure errors on duplicate outcomes", {
  models <- list(
    build_model("linear", formula = y ~ lag(x)),
    build_model("linear", formula = y ~ lag(x))
  )
  expect_error(
    validate_system_closure(models, c("gwcode", "year", "x", "y")),
    "Duplicate model outcomes.*y"
  )
})

test_that("validate_system_closure accepts a predictor produced by another model", {
  models <- list(
    build_model("exogen", formula = ~x),
    build_model("linear", formula = y ~ lag(x)),
    build_model("linear", formula = z ~ lag(y))
  )
  # x via exogen, y via model 2; z's lag(y) and the lagged x are all produced.
  expect_true(validate_system_closure(models, c("gwcode", "year", "x", "y", "z")))
})

test_that("validate_system_closure errors when an outcome column is missing from data", {
  models <- list(
    build_model("deterministic", formula = gdp ~ I(abs(gdppc * population)))
  )
  expect_error(
    validate_system_closure(models, c("gwcode", "year", "gdppc", "population")),
    "missing from the input data.*gdp"
  )
})

test_that("validate_system_closure errors when an exogen variable is missing from data", {
  models <- list(
    build_model("exogen", formula = ~population)
  )
  expect_error(
    validate_system_closure(models, c("gwcode", "year")),
    "missing from the input data.*population"
  )
})

test_that("validate_system_closure handles independent models", {
  models <- list(
    build_model("exogen", formula = ~x),
    build_model("linear", formula = y ~ lag(x))
  )
  expect_true(validate_system_closure(models, c("gwcode", "year", "x", "y")))
})

test_that("validate_system_closure errors on an unmodeled same-period predictor", {
  models <- list(
    build_model("exogen", formula = ~x),
    build_model("linear", formula = y ~ region + lag(x))
  )
  expect_error(
    validate_system_closure(models, c("gwcode", "year", "x", "y", "region")),
    "not produced by any model: region"
  )
})

test_that("validate_system_closure errors on an unmodeled lagged predictor", {
  # lag(psecprop) where psecprop has no producing model: NA from forecast step 2.
  models <- list(
    build_model("linear", formula = y ~ lag(psecprop))
  )
  expect_error(
    validate_system_closure(models, c("gwcode", "year", "y", "psecprop")),
    "not produced by any model: psecprop"
  )
})

test_that("validate_system_closure flags factor()-wrapped predictors", {
  models <- list(
    build_model("exogen", formula = ~x),
    build_model("linear", formula = y ~ -1 + factor(region) + lag(x))
  )
  expect_error(
    validate_system_closure(models, c("gwcode", "year", "x", "y", "region")),
    "not produced by any model.*region"
  )
})

test_that("validate_system_closure accepts a predictor produced by exogen", {
  models <- list(
    build_model("exogen", formula = ~region),
    build_model("linear", formula = y ~ region)
  )
  expect_true(validate_system_closure(models, c("gwcode", "year", "y", "region")))
})

test_that("validate_system_closure accepts a same-period predictor produced by another model", {
  models <- list(
    build_model("linear", formula = z ~ lag(z)),
    build_model("linear", formula = y ~ z)
  )
  expect_true(validate_system_closure(models, c("gwcode", "year", "y", "z")))
})

test_that("validate_system_closure accepts a lagged self-reference", {
  models <- list(
    build_model("deterministic", formula = gdppc ~ I(abs(lag(gdppc) * 1.02)))
  )
  expect_true(validate_system_closure(models, c("gwcode", "year", "gdppc")))
})

test_that("validate_system_closure flags an unmodeled heterolm variance predictor", {
  models <- list(
    build_model("exogen", formula = ~x),
    build_model("heterolm", formula = y ~ lag(x), variance = ~ w)
  )
  expect_error(
    validate_system_closure(models, c("gwcode", "year", "x", "y", "w")),
    "not produced by any model: w"
  )
})


# ── Integration: setup_system with bad input ────────────────────────────────

test_that("setup_system catches non-contiguous time series", {
  dt <- make_panel()
  dt <- dt[!(gwcode == 1 & year == 2005)]
  models <- list(build_model("linear", formula = y ~ lag(x)))
  expect_error(
    setup_system(
      models = models, data = dt,
      train_start = 2000, test_start = 2008, horizon = 2,
      groupvar = "gwcode", timevar = "year", inner_sims = 1
    ),
    "Non-contiguous"
  )
})

test_that("setup_system catches unclosed system", {
  dt <- make_panel()
  models <- list(build_model("linear", formula = y ~ lag(missing_var)))
  expect_error(
    setup_system(
      models = models, data = dt,
      train_start = 2000, test_start = 2008, horizon = 2,
      groupvar = "gwcode", timevar = "year", inner_sims = 1
    ),
    "missing from the input data.*missing_var"
  )
})

test_that("setup_system errors on an unmodeled same-period predictor with guidance", {
  dt <- make_panel()
  dt$region <- ifelse(dt$gwcode == 1, "A", "B")
  models <- list(build_model("linear", formula = y ~ region + lag(x)))
  expect_error(
    setup_system(
      models = models, data = dt,
      train_start = 2000, test_start = 2008, horizon = 2,
      groupvar = "gwcode", timevar = "year", inner_sims = 1
    ),
    "not produced by any model: region"
  )
})
