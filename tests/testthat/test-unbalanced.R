# Unbalanced panels: ragged entry/exit through setup, fit, and simulate ------
#
# Uses sim_panel_ragged() from helper-dgp.R. Units without a row at
# test_start - 1 are training-only; origin units must supply the history
# depth the model formulas require.

test_that("setup_system handles late entry and early exit", {
  dt <- sim_panel_ragged(units = 4, n_time = 40, seed = 1,
                         enter = c("3" = 15), exit = c("4" = 20))

  expect_message(
    setup <- setup_system(
      list(build_model("linear", formula = y ~ lag(y) + lag(x), boot = "resid"),
           build_model("exogen", formula = ~x)),
      dt, train_start = 1, test_start = 30, horizon = 5,
      groupvar = "unit", timevar = "time", inner_sims = 1),
    "training only")

  # The exited unit is not simulated but stays in the training data.
  expect_false(4 %in% unique(setup$simulation_data$unit))
  expect_true(4 %in% unique(setup$train_data$unit))

  # The late entrant is simulated; its grid rows span the full time range with
  # NA before entry (inert: only t >= test_start is written by the engine).
  sim3 <- setup$simulation_data[unit == 3]
  expect_equal(min(sim3$time), 1L)
  expect_true(sim3[time < 15, all(is.na(y))])
  expect_false(sim3[time >= 15 & time < 30, anyNA(y)])
})

test_that("setup_system errors when an origin unit enters too late for the lag depth", {
  dt <- sim_panel_ragged(units = 4, n_time = 40, seed = 2,
                         enter = c("4" = 29))   # enters at test_start - 1
  expect_error(
    setup_system(
      list(build_model("linear", formula = y ~ lag(lag(y)) + lag(x)),
           build_model("exogen", formula = ~x)),
      dt, train_start = 1, test_start = 30, horizon = 2,
      groupvar = "unit", timevar = "time", inner_sims = 1),
    "enter too late")
})

test_that("ragged panel end-to-end: forecasts cover origin units with no NAs", {
  future::plan(future::sequential)
  on.exit(future::plan(future::sequential), add = TRUE)
  dt <- sim_panel_ragged(units = 6, n_time = 40, seed = 3,
                         enter = c("5" = 15), exit = c("6" = 20))

  suppressMessages(
    setup <- setup_system(
      list(build_model("linear", formula = y ~ lag(y) + lag(x), boot = "resid"),
           build_model("exogen", formula = ~x)),
      dt, train_start = 1, test_start = 30, horizon = 5,
      groupvar = "unit", timevar = "time", inner_sims = 2))
  set.seed(13)
  fit <- fit_system(setup, nsim = 2)
  res <- simulate_system(fit)

  fc <- res[time >= 30]
  # Exactly (origin units) x horizon x (nsim * inner_sims) draws, no NA.
  expect_equal(nrow(fc), 5L * 5L * 2L * 2L)
  expect_setequal(unique(fc$unit), 1:5)
  expect_false(anyNA(fc$y))

  acc <- get_accuracy(res, "y", dt, test_start = 30)
  expect_true(all(is.finite(acc$crps)))
})

test_that("random-window fits succeed on ragged data", {
  future::plan(future::sequential)
  on.exit(future::plan(future::sequential), add = TRUE)
  dt <- sim_panel_ragged(units = 6, n_time = 40, seed = 4,
                         enter = c("5" = 20), exit = c("6" = 25))

  suppressMessages(
    setup <- setup_system(
      list(build_model("linear", formula = y ~ lag(y) + lag(x), boot = "resid"),
           build_model("exogen", formula = ~x)),
      dt, train_start = 1, test_start = 30, horizon = 3,
      groupvar = "unit", timevar = "time", inner_sims = 1, min_window = 8))
  set.seed(21)
  fit <- fit_system(setup, nsim = 6)

  # Windows predating unit 5's entry pool whatever rows exist; every draw
  # still produces finite estimates.
  co <- get_coefficients(fit)
  expect_true(all(is.finite(co$estimate)))
  expect_equal(data.table::uniqueN(co$.draw), 6L)
})
