library(testthat)
library(data.table)

test_that("center() stores the training mean and reapplies it at predict", {
  set.seed(1)
  dt <- data.table::CJ(unit = 1:6, time = 1:20)
  dt[, x := stats::rnorm(.N, 10, 1)]
  dt[, z := stats::rnorm(.N)]
  dt[, y := 0.5 * x + z + stats::rnorm(.N, 0, 0.3)]
  ctx <- panel_context(unit = "unit", time = "time")

  train <- dt[time <= 15L]
  m <- linearmodel(y ~ center(lag(x)) * z, data = train, ctx = ctx)

  m_train <- mean(train[order(unit, time), data.table::shift(x), by = unit]$V1,
                  na.rm = TRUE)
  centers <- Filter(function(e) identical(endogenr:::.pt_call_name(e), "center"),
                    m$ts_map)
  expect_length(centers, 1L)
  expect_equal(centers[[1L]]$center, m_train)
  expect_setequal(m$coefs$term,
                  c("(Intercept)", "center_lag_x", "z", "center_lag_x:z"))

  # Shift x in the conditioning period far from the training mean: if the mean
  # were recomputed from prediction-window data the result would differ.
  sim <- data.table::copy(dt)
  sim[, sim := 1L]
  sim[time == 15L, x := x + 100]
  pred <- predict(m, data = sim, t = 16L,
                  ctx = panel_context(unit = "unit", time = "time", sim = "sim"),
                  what = "expectation")
  pred <- pred[order(unit)]

  b  <- stats::coef(m$fitted)
  xc <- sim[time == 15L][order(unit), x] - m_train
  zz <- sim[time == 16L][order(unit), z]
  expected <- b[["(Intercept)"]] + b[["center_lag_x"]] * xc + b[["z"]] * zz +
    b[["center_lag_x:z"]] * xc * zz
  expect_equal(pred$y, unname(expected))
})

test_that("center() misuse errors clearly", {
  set.seed(2)
  dt <- data.table::CJ(unit = 1:3, time = 1:10)
  dt[, x := stats::rnorm(.N)]
  dt[, y := stats::rnorm(.N)]
  ctx <- panel_context(unit = "unit", time = "time")

  expect_error(linearmodel(y ~ lag(center(x)), data = dt, ctx = ctx),
               regexp = "inside the time-series function lag")
  expect_error(linearmodel(center(y) ~ x, data = dt, ctx = ctx),
               regexp = "right-hand side")
  expect_error(build_model("deterministic", formula = w ~ I(center(x))),
               regexp = "only supported in estimated model formulas")
  expect_error(center(1:3), regexp = "resolved by endogenr")
})
