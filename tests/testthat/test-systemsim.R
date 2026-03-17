# Tests for R/systemsim.R --------------------------------------------------

# ── Helper: minimal tsibble ──────────────────────────────────────────────

make_tiny_tsibble <- function() {
  df <- expand.grid(gwcode = c(1, 2), year = 2000:2015)
  set.seed(1)
  df$x <- rnorm(nrow(df))
  df$y <- 0.5 * df$x + rnorm(nrow(df), sd = 0.2)
  tsibble::tsibble(df, key = "gwcode", index = "year")
}

# ── prepare_simulation_data ──────────────────────────────────────────────

test_that("prepare_simulation_data returns a tsibble with correct dimensions", {
  ts <- make_tiny_tsibble()

  sim_data <- prepare_simulation_data(
    data       = ts,
    groupvar   = "gwcode",
    timevar    = "year",
    train_start = 2000,
    test_start  = 2010,
    horizon     = 3,
    inner_sims  = 2
  )

  expect_s3_class(sim_data, "tbl_ts")

  # Should span 2000–2012 (train_start to test_start + horizon - 1)
  expect_equal(range(sim_data$year), c(2000, 2012))

  # 2 groups × 13 years × 2 inner_sims = 52 rows
  expect_equal(nrow(sim_data), 2 * 13 * 2)

  # Training rows should have non-NA y; forecast rows should be NA
  train_rows <- sim_data[sim_data$year < 2010, ]
  expect_false(all(is.na(train_rows$y)))

  forecast_rows <- sim_data[sim_data$year >= 2010, ]
  expect_true(all(is.na(forecast_rows$y)))
})

test_that("prepare_simulation_data includes sim column as key", {

  ts <- make_tiny_tsibble()
  sim_data <- prepare_simulation_data(
    data       = ts,
    groupvar   = "gwcode",
    timevar    = "year",
    train_start = 2000,
    test_start  = 2010,
    horizon     = 1,
    inner_sims  = 5
  )

  expect_true("sim" %in% names(sim_data))
  expect_equal(sort(unique(sim_data$sim)), 1:5)
  expect_true("sim" %in% tsibble::key_vars(sim_data))
})

# ── process_independent_models ───────────────────────────────────────────

test_that("process_independent_models validates test_start", {
  ts <- make_tiny_tsibble()
  sim_data <- prepare_simulation_data(
    data = ts, groupvar = "gwcode", timevar = "year",
    train_start = 2000, test_start = 2010, horizon = 1, inner_sims = 1
  )

  # Should error because test_start is character
  expect_error(
    process_independent_models(sim_data, list(NULL), "2010", 1, 1),
    "numeric"
  )
})

# ── setup_simulator (integration) ────────────────────────────────────────

test_that("setup_simulator returns expected structure with linear model", {
  ts <- make_tiny_tsibble()

  model_system <- list(
    build_model("linear", formula = y ~ lag(x))
  )

  result <- setup_simulator(
    models      = model_system,
    data        = ts,
    train_start = 2000,
    test_start  = 2010,
    horizon     = 3,
    groupvar    = "gwcode",
    timevar     = "year",
    inner_sims  = 2
  )

  expect_type(result, "list")
  expect_true(all(c(
    "simulation_data", "models", "fitted_models",
    "test_start", "horizon", "execution_order",
    "groupvar", "timevar", "inner_sims"
  ) %in% names(result)))

  expect_equal(result$test_start, 2010)
  expect_equal(result$horizon, 3)
  expect_equal(result$inner_sims, 2)
  expect_s3_class(result$simulation_data, "tbl_ts")
})

test_that("setup_simulator fitted_models have coefs and gof for linear", {
  ts <- make_tiny_tsibble()

  model_system <- list(
    build_model("linear", formula = y ~ lag(x))
  )

  result <- setup_simulator(
    models      = model_system,
    data        = ts,
    train_start = 2000,
    test_start  = 2010,
    horizon     = 1,
    groupvar    = "gwcode",
    timevar     = "year",
    inner_sims  = 1
  )

  fm <- result$fitted_models[[1]]
  expect_true(!is.null(fm$coefs))
  expect_true(!is.null(fm$gof))
  expect_true("estimate" %in% names(fm$coefs))
  expect_true("r.squared" %in% names(fm$gof))
})

test_that("setup_simulator builds correct execution_order", {
  ts <- make_tiny_tsibble()

  model_system <- list(
    build_model("linear", formula = y ~ lag(x))
  )

  result <- setup_simulator(
    models      = model_system,
    data        = ts,
    train_start = 2000,
    test_start  = 2010,
    horizon     = 1,
    groupvar    = "gwcode",
    timevar     = "year",
    inner_sims  = 1
  )

  expect_true("y" %in% result$execution_order)
})

test_that("setup_simulator works with deterministic model", {
  ts <- make_tiny_tsibble()

  model_system <- list(
    build_model("deterministic", formula = y ~ I(x * 2))
  )

  result <- setup_simulator(
    models      = model_system,
    data        = ts,
    train_start = 2000,
    test_start  = 2010,
    horizon     = 1,
    groupvar    = "gwcode",
    timevar     = "year",
    inner_sims  = 1
  )

  expect_type(result, "list")
  expect_true("y" %in% result$execution_order)
})

test_that("setup_simulator handles multiple models with dependency ordering", {
  df <- expand.grid(gwcode = c(1, 2), year = 2000:2015)
  set.seed(1)
  df$x <- rnorm(nrow(df))
  df$y <- 0.5 * df$x + rnorm(nrow(df), sd = 0.2)
  df$z <- 0.3 * df$y + rnorm(nrow(df), sd = 0.1)
  ts <- tsibble::tsibble(df, key = "gwcode", index = "year")

  model_system <- list(
    build_model("linear", formula = y ~ lag(x)),
    build_model("linear", formula = z ~ lag(y))
  )

  result <- setup_simulator(
    models      = model_system,
    data        = ts,
    train_start = 2000,
    test_start  = 2010,
    horizon     = 1,
    groupvar    = "gwcode",
    timevar     = "year",
    inner_sims  = 1
  )

  # y must come before z in execution order
  ord <- result$execution_order
  expect_true(which(ord == "y") < which(ord == "z"))
})

# ── simulate_endogenr (minimal integration) ──────────────────────────────

test_that("simulate_endogenr runs and returns expected output", {
  ts <- make_tiny_tsibble()

  model_system <- list(
    build_model("linear", formula = y ~ lag(x))
  )

  setup <- setup_simulator(
    models      = model_system,
    data        = ts,
    train_start = 2000,
    test_start  = 2013,
    horizon     = 2,
    groupvar    = "gwcode",
    timevar     = "year",
    inner_sims  = 2
  )

  set.seed(42)
  res <- simulate_endogenr(nsim = 1, simulator_setup = setup, parallel = FALSE)

  expect_s3_class(res, "tbl_df")
  expect_true(".sim" %in% names(res))
  expect_true("y" %in% names(res))
  expect_true("gwcode" %in% names(res))
  expect_true("year" %in% names(res))

  # Should have rows for both training and forecast periods
  expect_true(all(c(2013, 2014) %in% res$year))
})

# ── sim_to_dist ──────────────────────────────────────────────────────────

test_that("sim_to_dist nests simulations into distributional objects", {
  ts <- make_tiny_tsibble()

  setup <- setup_simulator(
    models      = list(build_model("linear", formula = y ~ lag(x))),
    data        = ts,
    train_start = 2000,
    test_start  = 2013,
    horizon     = 2,
    groupvar    = "gwcode",
    timevar     = "year",
    inner_sims  = 2
  )

  set.seed(42)
  res <- simulate_endogenr(nsim = 1, simulator_setup = setup, parallel = FALSE)
  res_ts <- tsibble::tsibble(res, key = c("gwcode", ".sim"), index = "year") |>
    dplyr::filter(year >= 2013)

  dist_result <- sim_to_dist(res_ts, "y")

  expect_s3_class(dist_result, "tbl_ts")
  # y column should now be distributional
  expect_true(inherits(dist_result$y[[1]], "distribution"))
})
