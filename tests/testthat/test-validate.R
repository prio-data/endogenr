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

test_that("validate_panel errors on unbalanced panel", {
  dt <- make_panel()
  # Add extra year for unit 1
  extra <- data.table::data.table(gwcode = 1, year = 2011, x = 0, y = 0)
  dt <- rbind(dt, extra)
  data.table::setkeyv(dt, c("gwcode", "year"))
  ctx <- make_ctx()
  expect_error(validate_panel(dt, ctx, test_start = 2008), "Unbalanced panel")
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
    "not closed.*z"
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

test_that("validate_system_closure accepts RHS var provided by another model", {
  models <- list(
    build_model("linear", formula = y ~ lag(x)),
    build_model("linear", formula = z ~ lag(y))
  )
  # y is an outcome of model 1, so z ~ lag(y) should be fine even if y isn't in data_columns
  expect_true(validate_system_closure(models, c("gwcode", "year", "x")))
})

test_that("validate_system_closure handles independent models", {
  models <- list(
    build_model("exogen", formula = ~x),
    build_model("linear", formula = y ~ lag(x))
  )
  expect_true(validate_system_closure(models, c("gwcode", "year", "x")))
})

# ── Integration: setup_simulator with bad input ────────────────────────────

test_that("setup_simulator catches non-contiguous time series", {
  dt <- make_panel()
  dt <- dt[!(gwcode == 1 & year == 2005)]
  models <- list(build_model("linear", formula = y ~ lag(x)))
  expect_error(
    setup_simulator(
      models = models, data = dt,
      train_start = 2000, test_start = 2008, horizon = 2,
      groupvar = "gwcode", timevar = "year", inner_sims = 1
    ),
    "Non-contiguous"
  )
})

test_that("setup_simulator catches unclosed system", {
  dt <- make_panel()
  models <- list(build_model("linear", formula = y ~ lag(missing_var)))
  expect_error(
    setup_simulator(
      models = models, data = dt,
      train_start = 2000, test_start = 2008, horizon = 2,
      groupvar = "gwcode", timevar = "year", inner_sims = 1
    ),
    "not closed.*missing_var"
  )
})
