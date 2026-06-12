# Fit parity with a plain pooled regression ---------------------------------
#
# Pins the contract that endogenr's estimation sample and coefficients match a
# plain lm() on the same training window:
#   1. boot = NULL, full window  -> coefficients equal lm() at 1e-8, even with
#      NAs in unrelated panel columns (model-column pruning).
#   2. boot = "resid"            -> every draw estimates on lm()'s sample
#      (regression for the whole-frame na.omit defect) and centers on lm().
#   3. call LHS (log(y2))        -> bootstrap draws actually vary
#      (regression for the .boot_y response assignment).
#   4. get_train_window()        -> end always <= test_start - 1; min_window
#      equal to the full range returns the full window.

make_junk_panel <- function() {
  dt <- sim_panel_ar1(units = 30, n_time = 40, seed = 1)
  set.seed(99)
  dt[, junk := stats::rnorm(.N)]
  dt[sample(.N, 300), junk := NA]
  dt
}

ref_lm <- function(dt, test_start) {
  ref_dt <- data.table::copy(dt)
  ref_dt[, `:=`(lag_y = data.table::shift(y),
                lag_x = data.table::shift(x)), by = unit]
  stats::lm(y ~ lag_y + lag_x, ref_dt[time < test_start])
}

test_that("full-window fit equals plain lm() despite NAs in unrelated columns", {
  future::plan(future::sequential)
  on.exit(future::plan(future::sequential), add = TRUE)
  dt <- make_junk_panel()

  setup <- setup_system(
    list(build_model("linear", formula = y ~ lag(y) + lag(x)),
         build_model("exogen", formula = ~x)),
    dt, train_start = 1, test_start = 35, horizon = 3,
    groupvar = "unit", timevar = "time", inner_sims = 1)   # min_window = NULL
  fit <- fit_system(setup, nsim = 1)

  co  <- get_coefficients(fit)
  est <- stats::setNames(co$estimate, co$term)
  refc <- stats::coef(ref_lm(dt, 35))

  expect_setequal(names(est), names(refc))
  nm <- sort(names(refc))
  expect_equal(est[nm], refc[nm], tolerance = 1e-8)
})

test_that("bootstrap draws estimate on lm()'s sample and center on its fit", {
  future::plan(future::sequential)
  on.exit(future::plan(future::sequential), add = TRUE)
  dt <- make_junk_panel()
  ref <- ref_lm(dt, 35)

  setup <- setup_system(
    list(build_model("linear", formula = y ~ lag(y) + lag(x), boot = "resid"),
         build_model("exogen", formula = ~x)),
    dt, train_start = 1, test_start = 35, horizon = 3,
    groupvar = "unit", timevar = "time", inner_sims = 1)
  set.seed(7)
  fit <- fit_system(setup, nsim = 50)

  # Every draw refits on exactly lm()'s estimation sample: complete cases on
  # the model's own columns, untouched by NAs in the unrelated junk column.
  nobs_draws <- vapply(fit$fitted_draws,
                       function(d) stats::nobs(d[[1]]$fitted), numeric(1))
  expect_true(all(nobs_draws == stats::nobs(ref)))

  lagx <- get_coefficients(fit)[term == "lag_x", estimate]
  expect_lt(abs(mean(lagx) - unname(stats::coef(ref)["lag_x"])), 0.02)
})

test_that("call-LHS bootstrap draws vary across refits", {
  future::plan(future::sequential)
  on.exit(future::plan(future::sequential), add = TRUE)
  dt <- sim_panel_ar1(units = 30, n_time = 40, seed = 2)
  dt[, y2 := exp(as.numeric(scale(y)))]

  setup <- setup_system(
    list(build_model("linear", formula = log(y2) ~ lag(x), boot = "resid"),
         build_model("exogen", formula = ~x)),
    dt, train_start = 1, test_start = 35, horizon = 2,
    groupvar = "unit", timevar = "time", inner_sims = 1)
  set.seed(11)
  fit <- fit_system(setup, nsim = 8)

  lagx <- get_coefficients(fit)[term == "lag_x", estimate]
  expect_gt(stats::sd(lagx), 0)
})

test_that("get_train_window caps end at the forecast origin", {
  set.seed(42)
  for (i in 1:200) {
    w <- get_train_window(1965, 2015, min_window = 10)
    expect_lte(w$end, 2014)                  # always <= test_start - 1
    expect_gte(w$end - w$start + 1, 10)      # at least min_window long
  }
  # min_window equal to the full range returns the full window (no sampling).
  w <- get_train_window(1965, 2015, min_window = 50)
  expect_equal(w$start, 1965)
  expect_equal(w$end, 2014)
})
