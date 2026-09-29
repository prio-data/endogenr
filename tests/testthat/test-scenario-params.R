# Tests for the baked scenario-parameter framework ----------------------------
#
# Covers:
#   1. hist_mean() semantics and .required_history registration
#   2. mundlak_means() spec builder
#   3. mundlak_means() updates in-sim (expanding and windowed)
#   4. Linear model with factor(time) simulates without error
#   5. fe_resample (default) adds between-trajectory variance vs fe_fixed(0)
#   6. setup_param() structure — new baked contract
#   7. print.endogenr_scenario_params runs without error
#   8. fe_distribution simulates without error
#   9. Coefficient overrides (scalar, vector, function, guarded)
#  10. Non-linear override errors in coef_override()
#  11. No regression: scenario=NULL equals no scenario arg
#  12. Determinism: same seed -> identical baked values
#  13. Convergence: fe_converge shrinks values toward target
#  14. Unit-FE persist identity: default == explicit fe_persist(active=FALSE)
#  15. fe_converge on unit FE changes simulation results

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

  fcast <- res[res$time >= 16]
  expect_false(any(is.na(fcast$y)))

  t16 <- unique(res[res$time == 16, .(unit, m_x)])
  for (u in unique(dt$unit)) {
    all_x    <- dt[dt$unit == u, x]
    expected <- mean(all_x[1:16])
    expect_true(all(abs(t16[t16$unit == u, ]$m_x - expected) < 1e-10))
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

# ── default time-FE sampler vs fe_fixed(0) ───────────────────────────────────

test_that("default time-FE sampler adds between-trajectory variance relative to fe_fixed(0)", {
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

  # Default run: internal setup_param bakes fe_ar under caller's seed.
  set.seed(10)
  res_rs <- simulate_system(fit)

  # Fixed-zero run: same simulation randomness, no time-FE variance.
  dp <- setup_param(fit)
  dp$y$time_fe <- fe_fixed(dp$y$time_fe, value = 0)
  set.seed(10)
  res_fix <- simulate_system(fit, scenario_params = dp)

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

# ── setup_param structure — baked contract ───────────────────────────────────

test_that("setup_param returns baked structure for a TWFE model", {
  skip_on_cran()
  nsim <- 3L; inner_sims <- 2L; horizon <- 4L
  dt   <- sim_panel_common_shock(units = 10L, n_time = 20L, seed = 3L)
  sys  <- setup_system(
    list(build_model("linear", formula = y ~ lag(y) + factor(time))),
    data        = dt,
    train_start = 1, test_start = 16, horizon = horizon,
    groupvar    = "unit", timevar   = "time",
    inner_sims  = inner_sims
  )
  fit <- fit_system(sys, nsim = nsim)
  dp  <- setup_param(fit)

  expect_s3_class(dp, "endogenr_scenario_params")
  expect_true("y" %in% names(dp))
  expect_equal(dp$y$type, "linear")

  # time_fe block
  expect_equal(dp$y$time_fe$kind, "time")
  expect_false(is.null(dp$y$time_fe$term))
  n_fit_levels <- length(unique(dt[dt$time >= 2 & dt$time < 16, time]))
  expect_equal(length(dp$y$time_fe$levels), n_fit_levels)

  # Baked values shape: c(nsim, inner_sims, horizon)
  expect_equal(dim(dp$y$time_fe$values), c(nsim, inner_sims, horizon))

  # dims attribute
  d <- attr(dp, "dims")
  expect_equal(d$nsim, nsim)
  expect_equal(d$inner_sims, inner_sims)
  expect_equal(d$horizon, horizon)
})

# ── print runs without error ─────────────────────────────────────────────────

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

# ── fe_distribution ──────────────────────────────────────────────────────────

test_that("fe_distribution simulates without error", {
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
  dp$y$time_fe <- fe_distribution(dp$y$time_fe, sd = 0.01)

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
  dp$y$coef <- coef_override(dp$y$coef, "m_x", v)

  model <- Filter(function(m) !is.null(m$outcome) && m$outcome == "y",
                  fit$fitted_draws[[1L]])[[1L]]
  ctx   <- fit$ctx
  sd    <- fit$simulation_data

  # Use the slice for draw 1 as the scenario argument to predict.linear
  scen_slice <- endogenr:::.slice_scenario(dp, 1L)

  out_ov <- predict(model, data = sd, t = 14L, ctx = ctx,
                    what = "expectation", scenario = scen_slice, test_start = 14L)

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
  dp$y$coef <- coef_override(dp$y$coef, "m_x", betas)

  model <- Filter(function(m) !is.null(m$outcome) && m$outcome == "y",
                  fit$fitted_draws[[1L]])[[1L]]
  ctx   <- fit$ctx
  sd    <- fit$simulation_data

  scen_slice <- endogenr:::.slice_scenario(dp, 1L)

  # t=14 → h=1 (beta=10), t=15 → h=2 (beta=20); both in-sample, all populated.
  out1 <- predict(model, data = sd, t = 14L, ctx = ctx,
                  what = "expectation", scenario = scen_slice, test_start = 14L)
  out2 <- predict(model, data = sd, t = 15L, ctx = ctx,
                  what = "expectation", scenario = scen_slice, test_start = 14L)

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
  dp$y$coef <- coef_override(dp$y$coef, "m_x", function(h, beta_hat) beta_hat * (1 + h))

  model <- Filter(function(m) !is.null(m$outcome) && m$outcome == "y",
                  fit$fitted_draws[[1L]])[[1L]]
  ctx   <- fit$ctx
  sd    <- fit$simulation_data
  b0    <- stats::coef(model$fitted)[["m_x"]]

  scen_slice <- endogenr:::.slice_scenario(dp, 1L)

  # h=1 at t=14 → beta* = b0 * 2; h=2 at t=15 → beta* = b0 * 3
  out_h1 <- predict(model, data = sd, t = 14L, ctx = ctx,
                    what = "expectation", scenario = scen_slice, test_start = 14L)
  out_h2 <- predict(model, data = sd, t = 15L, ctx = ctx,
                    what = "expectation", scenario = scen_slice, test_start = 14L)

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

# ── Non-linear override errors in coef_override() ────────────────────────────

test_that("coef_override on non-linear model's coef block errors", {
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
  dp  <- setup_param(fit)

  # dp$y$coef has linear = FALSE (it's a glm)
  expect_error(coef_override(dp$y$coef, "(Intercept)", 0.5), "linear")
})

# ── No regression for non-scenario linear models ─────────────────────────────

test_that("linear model without FE and no overrides is unchanged by scenario=NULL", {
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

# ── Determinism: same seed → identical baked values ──────────────────────────

test_that("setup_param is deterministic under the same seed", {
  skip_on_cran()
  dt  <- sim_panel_common_shock(units = 8L, n_time = 20L, seed = 12L)
  sys <- setup_system(
    list(build_model("linear", formula = y ~ lag(y) + factor(time))),
    data        = dt,
    train_start = 1, test_start = 16, horizon = 3,
    groupvar    = "unit", timevar   = "time",
    inner_sims  = 2L
  )
  fit <- fit_system(sys, nsim = 3L)

  set.seed(1)
  a <- setup_param(fit)
  set.seed(1)
  b <- setup_param(fit)

  expect_identical(a$y$time_fe$values, b$y$time_fe$values)
})

# ── Convergence: fe_converge shrinks toward target ───────────────────────────

test_that("fe_converge(to=0, path='linear') yields values shrinking toward 0", {
  skip_on_cran()
  dt  <- sim_panel_common_shock(units = 8L, n_time = 20L, seed = 13L)
  sys <- setup_system(
    list(build_model("linear", formula = y ~ lag(y) + factor(time))),
    data        = dt,
    train_start = 1, test_start = 16, horizon = 4,
    groupvar    = "unit", timevar   = "time",
    inner_sims  = 2L
  )
  fit <- fit_system(sys, nsim = 5L)
  set.seed(1)
  dp  <- setup_param(fit)
  blk <- dp$y$time_fe

  # Apply convergence toward 0
  blk_conv <- fe_converge(blk, to = 0, path = "linear")

  # Mean absolute value should decrease monotonically toward 0
  h1_abs <- mean(abs(blk_conv$values[,, 1L]))
  hH_abs <- mean(abs(blk_conv$values[,, 4L]))
  expect_true(hH_abs < h1_abs)
})

# ── Unit-FE persist identity ─────────────────────────────────────────────────

test_that("unit-FE model: default == explicit fe_persist (active=FALSE)", {
  skip_on_cran()
  dt  <- sim_panel_ar1(units = 8L, n_time = 20L, seed = 14L)
  sys <- setup_system(
    list(
      build_model("exogen", formula = ~x),
      build_model("linear", formula = y ~ lag(y) + x + factor(unit))
    ),
    data        = dt,
    train_start = 1, test_start = 16, horizon = 3,
    groupvar    = "unit", timevar   = "time",
    inner_sims  = 2L
  )
  fit <- fit_system(sys, nsim = 4L)

  # Explicit params with default fe_persist (active = FALSE)
  sp <- setup_param(fit)
  expect_false(isTRUE(sp$y$unit_fe$active))

  set.seed(7)
  r0 <- simulate_system(fit)
  set.seed(7)
  r1 <- simulate_system(fit, scenario_params = setup_param(fit))

  expect_equal(r0$y, r1$y)
})

# ── fe_converge on unit FE changes results ───────────────────────────────────

test_that("fe_converge on unit FE changes simulation results vs default", {
  skip_on_cran()
  dt  <- sim_panel_ar1(units = 8L, n_time = 20L, seed = 15L)
  sys <- setup_system(
    list(
      build_model("exogen", formula = ~x),
      build_model("linear", formula = y ~ lag(y) + x + factor(unit))
    ),
    data        = dt,
    train_start = 1, test_start = 16, horizon = 3,
    groupvar    = "unit", timevar   = "time",
    inner_sims  = 2L
  )
  fit <- fit_system(sys, nsim = 4L)

  dp <- setup_param(fit)
  dp$y$unit_fe <- fe_converge(dp$y$unit_fe, to = 0, path = "linear")
  expect_true(isTRUE(dp$y$unit_fe$active))

  set.seed(7)
  r_default  <- simulate_system(fit)
  set.seed(7)
  r_converge <- simulate_system(fit, scenario_params = dp)

  # Convergence toward zero changes unit-FE contribution → results must differ
  expect_false(isTRUE(all.equal(r_default$y, r_converge$y)))
})

# ── Time-FE parameter uncertainty + AR(1) sampler ────────────────────────────

.twfe_fit <- function(inner_sims = 2L, horizon = 4L, nsim = 3L, seed = 21L) {
  dt  <- sim_panel_common_shock(units = 10L, n_time = 20L, seed = seed)
  sys <- setup_system(
    list(build_model("linear", formula = y ~ lag(y) + factor(time))),
    data        = dt,
    train_start = 1, test_start = 16, horizon = horizon,
    groupvar    = "unit", timevar   = "time",
    inner_sims  = inner_sims
  )
  fit_system(sys, nsim = nsim)
}

test_that("inner_sims = 1 or horizon = 1 does not drop FE slice dimensions", {
  skip_on_cran()
  for (args in list(list(inner_sims = 1L, horizon = 4L),
                    list(inner_sims = 2L, horizon = 1L))) {
    fit <- do.call(.twfe_fit, args)
    set.seed(1)
    expect_no_error(res <- simulate_system(fit))
    fc <- res[res$time >= 16]
    expect_true(nrow(fc) > 0L)
    expect_true(all(is.finite(fc$y)))
  }

  dt  <- sim_panel_ar1(units = 8L, n_time = 20L, seed = 22L)
  sys <- setup_system(
    list(
      build_model("exogen", formula = ~x),
      build_model("linear", formula = y ~ lag(y) + x + factor(unit))
    ),
    data        = dt,
    train_start = 1, test_start = 16, horizon = 3,
    groupvar    = "unit", timevar   = "time",
    inner_sims  = 1L
  )
  fit <- fit_system(sys, nsim = 2L)
  set.seed(2)
  sp <- setup_param(fit)
  sp$y$unit_fe <- fe_resample(sp$y$unit_fe)
  expect_no_error(res <- simulate_system(fit, scenario_params = sp))
  expect_true(all(is.finite(res[res$time >= 16]$y)))

  expect_error(fe_ar(sp$y$unit_fe), "supports time fixed effects only")
})

test_that("time-FE vcov matches lm standard errors with a zero reference row", {
  skip_on_cran()
  fit <- .twfe_fit(nsim = 1L)
  m   <- fit$fitted_models[[1L]]
  V   <- m$time_fe$vcov
  lv  <- names(m$time_fe$effects)
  expect_equal(dimnames(V), list(lv, lv))

  se  <- m$coefs$std.error[match(paste0("factor(time)", lv[-1L]), m$coefs$term)]
  expect_equal(unname(sqrt(diag(V))[-1L]), se, tolerance = 1e-8)
  expect_true(all(V[1L, ] == 0))
  expect_true(all(V[, 1L] == 0))
})

test_that("fe_ar continues the AR(1) recursion from the origin, stepping gaps", {
  skip_on_cran()
  fit <- .twfe_fit()
  set.seed(3)
  blk <- setup_param(fit)$y$time_fe
  expect_identical(blk$heuristic$strategy, "ar1")
  H   <- blk$dims$horizon
  L   <- length(blk$effects_by_draw[[1L]])
  tau <- numeric(L + H + 2L)
  for (t in 2:length(tau)) tau[t] <- 1 + 0.5 * tau[t - 1L]

  with_levels <- function(blk, lv) {
    for (i in seq_along(blk$effects_by_draw)) {
      blk$effects_by_draw[[i]] <- stats::setNames(tau[seq_len(L)], lv)
      blk$vcov_by_draw[[i]]    <- matrix(0, L, L, dimnames = list(lv, lv))
    }
    blk
  }

  # Levels 2..15 end at the origin (test_start - 1 = 15): no gap.
  b0  <- with_levels(blk, names(blk$effects_by_draw[[1L]]))
  out <- fe_ar(b0)
  expected <- tau[L + seq_len(H)]
  for (i in seq_len(blk$dims$nsim)) for (s in seq_len(blk$dims$inner_sims)) {
    expect_equal(out$values[i, s, ], expected, tolerance = 1e-8)
  }
  expect_equal(as.vector(out$heuristic$rho), rep(0.5, length(out$heuristic$rho)),
               tolerance = 1e-8)

  # Last level 13: two gap years (14, 15) are stepped before the forecast.
  b2  <- with_levels(blk, as.character(seq(to = 13L, length.out = L)))
  out <- fe_ar(b2)
  expected <- tau[L + 2L + seq_len(H)]
  for (i in seq_len(blk$dims$nsim)) for (s in seq_len(blk$dims$inner_sims)) {
    expect_equal(out$values[i, s, ], expected, tolerance = 1e-8)
  }
})

test_that("time-FE samplers draw parameter uncertainty from the effect vcov", {
  skip_on_cran()
  fit <- .twfe_fit()
  set.seed(4)
  blk <- setup_param(fit)$y$time_fe
  lv  <- names(blk$effects_by_draw[[1L]])
  L   <- length(lv)
  blk$dims$inner_sims <- 400L

  zero_eff <- stats::setNames(numeric(L), lv)
  V <- diag(c(0, rep(0.25, L - 1L)))
  dimnames(V) <- list(lv, lv)
  for (i in seq_along(blk$effects_by_draw)) {
    blk$effects_by_draw[[i]] <- zero_eff
    blk$vcov_by_draw[[i]]    <- V
  }
  set.seed(1)
  # Mixture: reference level (0) w.p. 1/L, else N(0, 0.25).
  expect_equal(stats::sd(fe_resample(blk)$values), 0.5 * sqrt((L - 1) / L),
               tolerance = 0.05)

  for (i in seq_along(blk$vcov_by_draw)) blk$vcov_by_draw[[i]][] <- 0
  expect_true(all(fe_ar(blk)$values == 0))
})
