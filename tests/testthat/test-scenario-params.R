# Tests for TWFE scenario parameters (Steps 1-6) ----------------------------
#
# Covers:
#   1. hist_mean() semantics and .required_history registration
#   2. mundlak_means() spec builder
#   3. mundlak_means() updates in-sim (expanding and windowed)
#   4. Time-FE blocker fix (factor(time) in linear model)
#   5. Policy changes ensemble variance
#   6. setup_param() overview and print
#   7. Coefficient overrides (scalar, vector, function, guarded)
#   8. No regression for non-scenario models

future::plan(future::sequential)

# ── hist_mean: unit tests ────────────────────────────────────────────────────

test_that("hist_mean: expanding mean matches cumulative mean", {
  expect_equal(hist_mean(c(1, 2, 3, 4)), c(1, 1.5, 2, 2.5))
})

test_that("hist_mean: windowed mean (window = 2)", {
  expect_equal(hist_mean(c(1, 2, 3, 4), window = 2), c(1, 1.5, 2.5, 3.5))
})

test_that("hist_mean: NA propagation", {
  expect_equal(hist_mean(c(1, NA, 3)), c(1, 1, 2))
})

test_that("hist_mean: all NA with min_obs = 1 returns all NA", {
  out <- hist_mean(c(NA_real_, NA_real_), min_obs = 1L)
  expect_true(all(is.na(out)))
})

test_that("hist_mean: min_obs defers output until enough obs", {
  out <- hist_mean(c(1, 2, 3), min_obs = 2L)
  expect_true(is.na(out[1]))
  expect_equal(out[2], 1.5)
  expect_equal(out[3], 2)
})

test_that("hist_mean: validates window argument", {
  expect_error(hist_mean(1:3, window = 0),  "`window`")
  expect_error(hist_mean(1:3, window = -1), "`window`")
  expect_error(hist_mean(1:3, window = "a"), "`window`")
})

test_that("hist_mean: validates min_obs argument", {
  expect_error(hist_mean(1:3, min_obs = 0),  "`min_obs`")
  expect_error(hist_mean(1:3, min_obs = -1), "`min_obs`")
})

test_that(".required_history treats hist_mean as Inf", {
  expect_equal(.required_history(y ~ hist_mean(x)), Inf)
  expect_equal(.required_history(y ~ hist_mean(x, window = 3)), Inf)
})

# ── mundlak_means: build and validate ────────────────────────────────────────

test_that("mundlak_means returns list of deterministic specs", {
  specs <- mundlak_means("x")
  expect_true(is.list(specs))
  expect_equal(length(specs), 1L)
  expect_true(inherits(specs[[1L]], "deterministic_spec"))
  expect_match(deparse(specs[[1L]]$formula), "m_x")
})

test_that("mundlak_means: named vars produce custom output names", {
  specs <- mundlak_means(c(my_mean = "x"))
  f     <- specs[[1L]]$formula
  expect_equal(as.character(f[[2L]]), "my_mean")
})

test_that("mundlak_means: multiple vars produce multiple specs", {
  specs <- mundlak_means(c("x", "z"))
  expect_equal(length(specs), 2L)
  expect_equal(as.character(specs[[1L]]$formula[[2L]]), "m_x")
  expect_equal(as.character(specs[[2L]]$formula[[2L]]), "m_z")
})

test_that("mundlak_means: per-variable windows accepted", {
  specs <- mundlak_means(c("x", "z"), window = c(m_x = 3, m_z = Inf))
  expect_match(deparse(specs[[1L]]$formula), "window = 3")
  expect_match(deparse(specs[[2L]]$formula), "window = Inf")
})

test_that("mundlak_means: validates window entries", {
  expect_error(mundlak_means("x", window = 0),   "`window`")
  expect_error(mundlak_means("x", window = "a"), "`window`")
})

# ── Mundlak means update in-sim ─────────────────────────────────────────────

test_that("mundlak_means expanding mean updates during simulation", {
  skip_on_cran()
  dt <- sim_panel_ar1(units = 6L, n_time = 20L, seed = 42)
  # Pre-populate m_x with expanding mean so validate_panel sees a non-NA
  # initial state at test_start - 1 = 15; the deterministic model overwrites
  # it each forecast step.
  dt[, m_x := hist_mean(x), by = unit]

  sys <- setup_system(
    c(
      list(build_model("exogen", formula = ~x)),
      mundlak_means("x"),
      list(build_model("linear", formula = y ~ lag(y) + m_x))
    ),
    data        = dt,
    train_start = 1, test_start = 16, horizon = 4,
    groupvar    = "unit", timevar   = "time",
    inner_sims  = 2L
  )
  fit <- fit_system(sys, nsim = 2L)
  res <- simulate_system(fit)

  # Forecast rows should be finite
  fcast <- res[res$time >= 16]
  expect_false(any(is.na(fcast$y)))

  # At time 16 the deterministic model has consumed x[1..16] (x at t=16 is
  # available from the exogen model), so m_x[16] = mean(x[1..16]).
  t16 <- unique(res[res$time == 16, .(unit, m_x)])
  for (u in unique(dt$unit)) {
    all_x <- dt[dt$unit == u, x]          # x[1..20]
    expected_m_x <- mean(all_x[1:16])     # hist_mean at position 16
    expect_true(all(abs(t16[t16$unit == u, ]$m_x - expected_m_x) < 1e-10))
  }
})

test_that("mundlak_means windowed (window = 3) matches trailing 3-obs mean", {
  skip_on_cran()
  dt <- sim_panel_ar1(units = 6L, n_time = 20L, seed = 7)
  dt[, m_x := hist_mean(x, window = 3), by = unit]

  sys <- setup_system(
    c(
      list(build_model("exogen", formula = ~x)),
      mundlak_means("x", window = 3),
      list(build_model("linear", formula = y ~ lag(y) + m_x))
    ),
    data        = dt,
    train_start = 1, test_start = 16, horizon = 2,
    groupvar    = "unit", timevar   = "time",
    inner_sims  = 2L
  )
  fit <- fit_system(sys, nsim = 2L)
  res <- simulate_system(fit)

  # At time 16 the trailing 3-obs mean of x[1..16] = mean(x[14:16])
  t16 <- unique(res[res$time == 16, .(unit, m_x)])
  for (u in unique(dt$unit)) {
    all_x <- dt[dt$unit == u, x]
    last3 <- all_x[14:16]
    expect_true(all(abs(t16[t16$unit == u, ]$m_x - mean(last3)) < 1e-10))
  }
})

# ── Time-FE blocker fix ──────────────────────────────────────────────────────

test_that("linear model with factor(time) simulates without error", {
  skip_on_cran()
  dt  <- sim_panel_common_shock(units = 12L, n_time = 30L, seed = 1L)
  sys <- setup_system(
    list(build_model("linear", formula = y ~ lag(y) + factor(time))),
    data        = dt,
    train_start = 1, test_start = 25, horizon = 5,
    groupvar    = "unit", timevar   = "time",
    inner_sims  = 3L
  )
  fit <- fit_system(sys, nsim = 2L)

  expect_no_error(res <- simulate_system(fit))
  fcast <- res[res$time >= 25]
  expect_false(any(is.na(fcast$y)))
  expect_true(all(is.finite(fcast$y)))
})

# ── Policy changes ensemble variance ────────────────────────────────────────

test_that("resample policy adds between-trajectory variance vs fixed(0)", {
  skip_on_cran()
  dt  <- sim_panel_common_shock(units = 10L, n_time = 30L, seed = 2L)
  sys <- setup_system(
    list(build_model("linear", formula = y ~ lag(y) + factor(time))),
    data        = dt,
    train_start = 1, test_start = 25, horizon = 3,
    groupvar    = "unit", timevar   = "time",
    inner_sims  = 3L
  )
  fit <- fit_system(sys, nsim = 30L)

  # Resample run (default)
  set.seed(10)
  res_rs  <- simulate_system(fit)

  # Fixed-zero run
  dp <- setup_param(fit)
  dp$y$policy <- list(type = "fixed", value = 0)
  set.seed(10)
  res_fix <- simulate_system(fit, scenario_params = dp)

  # For each forecast year, compute across-.sim variance of the unit-mean
  var_rs  <- vapply(25:27, function(yr) {
    vals <- tapply(res_rs[res_rs$time == yr, y],
                   res_rs[res_rs$time == yr, .sim], mean)
    var(vals)
  }, numeric(1))

  var_fix <- vapply(25:27, function(yr) {
    vals <- tapply(res_fix[res_fix$time == yr, y],
                   res_fix[res_fix$time == yr, .sim], mean)
    var(vals)
  }, numeric(1))

  expect_true(any(var_rs > var_fix))
})

# ── setup_param overview ─────────────────────────────────────────────────────

test_that("setup_param returns expected structure for a TWFE model", {
  skip_on_cran()
  dt  <- sim_panel_common_shock(units = 10L, n_time = 20L, seed = 3L)
  sys <- setup_system(
    list(build_model("linear", formula = y ~ lag(y) + factor(time))),
    data        = dt,
    train_start = 1, test_start = 16, horizon = 4,
    groupvar    = "unit", timevar   = "time",
    inner_sims  = 2L
  )
  fit <- fit_system(sys, nsim = 2L)
  dp  <- setup_param(fit)

  expect_s3_class(dp, "endogenr_scenario_params")
  expect_true("y" %in% names(dp))
  expect_equal(dp$y$type, "linear")

  # n_levels equals the unique factor(time) levels in the fitted lm.
  # lag(y) drops time=1 from the fitting data, so levels = 2..(test_start-1).
  n_fit_levels <- length(unique(dt[dt$time >= 2 & dt$time < 16, time]))
  expect_equal(dp$y$time_fe$n_levels, n_fit_levels)

  expect_false(is.null(dp$y$time_fe))
  expect_true(all(c("term", "timevar", "n_levels", "effects", "summary") %in%
                    names(dp$y$time_fe)))
})

test_that("print.endogenr_scenario_params runs without error", {
  skip_on_cran()
  dt  <- sim_panel_common_shock(units = 6L, n_time = 20L, seed = 4L)
  sys <- setup_system(
    list(build_model("linear", formula = y ~ lag(y) + factor(time))),
    data        = dt,
    train_start = 1, test_start = 16, horizon = 3,
    groupvar    = "unit", timevar   = "time",
    inner_sims  = 2L
  )
  fit <- fit_system(sys, nsim = 2L)
  dp  <- setup_param(fit)

  expect_no_error(capture.output(print(dp)))
})

test_that("distribution policy simulates without error", {
  skip_on_cran()
  dt  <- sim_panel_common_shock(units = 6L, n_time = 20L, seed = 5L)
  sys <- setup_system(
    list(build_model("linear", formula = y ~ lag(y) + factor(time))),
    data        = dt,
    train_start = 1, test_start = 16, horizon = 3,
    groupvar    = "unit", timevar   = "time",
    inner_sims  = 2L
  )
  fit <- fit_system(sys, nsim = 4L)
  dp  <- setup_param(fit)
  dp$y$policy$type <- "distribution"

  expect_no_error(simulate_system(fit, scenario_params = dp))
})

# ── Coefficient overrides ────────────────────────────────────────────────────

test_that("scalar coefficient override: predict.linear matches hand-built check", {
  skip_on_cran()
  dt  <- sim_panel_ar1(units = 8L, n_time = 20L, seed = 6L)
  dt[, m_x := hist_mean(x), by = unit]

  sys <- setup_system(
    c(
      list(build_model("exogen", formula = ~x)),
      mundlak_means("x"),
      list(build_model("linear", formula = y ~ lag(y) + m_x))
    ),
    data        = dt,
    train_start = 1, test_start = 16, horizon = 4,
    groupvar    = "unit", timevar   = "time",
    inner_sims  = 2L
  )
  fit <- fit_system(sys, nsim = 2L)

  dp <- setup_param(fit)
  v  <- 999
  dp$y$coefficients$m_x <- v

  model <- Filter(function(m) !is.null(m$outcome) && m$outcome == "y",
                  fit$fitted_draws[[1L]])[[1L]]
  ctx   <- fit$ctx
  sd    <- fit$simulation_data

  # Use t=14 (in-sample; all predictors populated): test_start=14 → h=1.
  # This isolates the override mechanism from missing-data in forecast rows.
  out_ov <- predict(model, data = sd, t = 14L, ctx = ctx,
                    what = "expectation", scenario = dp, test_start = 14L)

  # Hand-built reference: manually swap coefficient and predict
  fit_hand              <- model$fitted
  fit_hand$coefficients[["m_x"]] <- v
  model_hand            <- model
  model_hand$fitted     <- fit_hand
  out_hand <- predict(model_hand, data = sd, t = 14L, ctx = ctx,
                      what = "expectation", scenario = NULL, test_start = NULL)

  expect_true(isTRUE(all.equal(out_ov$y, out_hand$y, tolerance = 1e-10)))
  # Sanity: override must differ from unmodified prediction
  out_base <- predict(model, data = sd, t = 14L, ctx = ctx,
                      what = "expectation", scenario = NULL, test_start = NULL)
  expect_false(isTRUE(all.equal(out_ov$y, out_base$y, tolerance = 1e-6)))
})

test_that("length-horizon vector override yields per-step beta*", {
  skip_on_cran()
  dt  <- sim_panel_ar1(units = 6L, n_time = 20L, seed = 7L)
  dt[, m_x := hist_mean(x), by = unit]

  sys <- setup_system(
    c(
      list(build_model("exogen", formula = ~x)),
      mundlak_means("x"),
      list(build_model("linear", formula = y ~ lag(y) + m_x))
    ),
    data        = dt,
    train_start = 1, test_start = 16, horizon = 4,
    groupvar    = "unit", timevar   = "time",
    inner_sims  = 2L
  )
  fit   <- fit_system(sys, nsim = 2L)
  dp    <- setup_param(fit)
  betas <- c(10, 20, 30, 40)
  dp$y$coefficients$m_x <- betas

  model <- Filter(function(m) !is.null(m$outcome) && m$outcome == "y",
                  fit$fitted_draws[[1L]])[[1L]]
  ctx   <- fit$ctx
  sd    <- fit$simulation_data

  # t=14 → h=1 (beta=10), t=15 → h=2 (beta=20); both in-sample, all populated.
  out1 <- predict(model, data = sd, t = 14L, ctx = ctx,
                  what = "expectation", scenario = dp, test_start = 14L)
  out2 <- predict(model, data = sd, t = 15L, ctx = ctx,
                  what = "expectation", scenario = dp, test_start = 14L)

  mk_ref <- function(b, t_val) {
    fc <- model$fitted
    fc$coefficients[["m_x"]] <- b
    mm <- model; mm$fitted <- fc
    predict(mm, data = sd, t = t_val, ctx = ctx,
            what = "expectation", scenario = NULL, test_start = NULL)
  }
  ref1 <- mk_ref(10, 14L)
  ref2 <- mk_ref(20, 15L)

  expect_true(isTRUE(all.equal(out1$y, ref1$y, tolerance = 1e-10)))
  expect_true(isTRUE(all.equal(out2$y, ref2$y, tolerance = 1e-10)))
  # Different betas → different predictions
  expect_false(isTRUE(all.equal(ref1$y, ref2$y, tolerance = 1e-6)))
})

test_that("function override is honoured at each step", {
  skip_on_cran()
  dt  <- sim_panel_ar1(units = 6L, n_time = 20L, seed = 8L)
  dt[, m_x := hist_mean(x), by = unit]

  sys <- setup_system(
    c(
      list(build_model("exogen", formula = ~x)),
      mundlak_means("x"),
      list(build_model("linear", formula = y ~ lag(y) + m_x))
    ),
    data        = dt,
    train_start = 1, test_start = 16, horizon = 3,
    groupvar    = "unit", timevar   = "time",
    inner_sims  = 2L
  )
  fit   <- fit_system(sys, nsim = 2L)
  dp    <- setup_param(fit)
  dp$y$coefficients$m_x <- function(h, beta_hat) beta_hat * (1 + h)

  model <- Filter(function(m) !is.null(m$outcome) && m$outcome == "y",
                  fit$fitted_draws[[1L]])[[1L]]
  ctx   <- fit$ctx
  sd    <- fit$simulation_data
  b0    <- stats::coef(model$fitted)[["m_x"]]

  # h=1 at t=14 → beta* = b0 * 2; h=2 at t=15 → beta* = b0 * 3
  out_h1 <- predict(model, data = sd, t = 14L, ctx = ctx,
                    what = "expectation", scenario = dp, test_start = 14L)
  out_h2 <- predict(model, data = sd, t = 15L, ctx = ctx,
                    what = "expectation", scenario = dp, test_start = 14L)

  ref1 <- { fc <- model$fitted; fc$coefficients[["m_x"]] <- b0 * 2
             mm <- model; mm$fitted <- fc
             predict(mm, data = sd, t = 14L, ctx = ctx,
                     what = "expectation", scenario = NULL, test_start = NULL) }
  ref2 <- { fc <- model$fitted; fc$coefficients[["m_x"]] <- b0 * 3
             mm <- model; mm$fitted <- fc
             predict(mm, data = sd, t = 15L, ctx = ctx,
                     what = "expectation", scenario = NULL, test_start = NULL) }

  expect_true(isTRUE(all.equal(out_h1$y, ref1$y, tolerance = 1e-10)))
  expect_true(isTRUE(all.equal(out_h2$y, ref2$y, tolerance = 1e-10)))
})

test_that("coefficient override on non-linear model errors in simulate_system", {
  skip_on_cran()
  dt  <- sim_panel_ar1(units = 6L, n_time = 20L, seed = 9L)
  sys <- setup_system(
    list(
      build_model("exogen",  formula = ~x),
      build_model("glm",     formula = y ~ lag(y) + x,
                  family = stats::gaussian())
    ),
    data        = dt,
    train_start = 1, test_start = 16, horizon = 3,
    groupvar    = "unit", timevar   = "time",
    inner_sims  = 2L
  )
  fit <- fit_system(sys, nsim = 2L)

  # Manually craft a scenario_params-like object with coef override on glm
  dp <- structure(
    list(y = list(
      type         = "glm_endogenr",
      time_fe      = NULL,
      policy       = list(type = "resample"),
      coef_table   = NULL,
      estimates    = c(`(Intercept)` = 0),
      coefficients = list(`(Intercept)` = 0.5)
    )),
    class = "endogenr_scenario_params"
  )

  expect_error(simulate_system(fit, scenario_params = dp),
               "only supported for `linear` models")
})

# ── No regression for non-scenario linear models ─────────────────────────────

test_that("linear model without factor(time) and no overrides is unchanged by scenario=NULL", {
  skip_on_cran()
  dt  <- sim_panel_ar1(units = 8L, n_time = 20L, seed = 11L)
  sys <- setup_system(
    list(
      build_model("exogen", formula = ~x),
      build_model("linear", formula = y ~ lag(y) + x)
    ),
    data        = dt,
    train_start = 1, test_start = 16, horizon = 3,
    groupvar    = "unit", timevar   = "time",
    inner_sims  = 2L
  )
  fit <- fit_system(sys, nsim = 4L)

  set.seed(42)
  res_base <- simulate_system(fit)
  set.seed(42)
  res_sc   <- simulate_system(fit, scenario_params = NULL)

  expect_equal(res_base$y, res_sc$y)
})
