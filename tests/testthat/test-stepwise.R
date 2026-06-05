# Pillar 4 — step-wise dynamic multi-horizon forecasting --------------------
#
# The core of the engine: an exact deterministic recurrence, error propagation
# that grows with the horizon, a true model out-scoring a misspecified one, the
# spatial-lag transform, and reproducibility.

# ── exact deterministic recurrence over the full horizon ───────────────────

test_that("a deterministic system reproduces its closed-form recurrence", {
  future::plan(future::sequential)
  on.exit(future::plan(future::sequential), add = TRUE)

  set.seed(1)
  df <- data.table::as.data.table(expand.grid(unit = 1:3, time = 1:20))
  df[, x := stats::rnorm(.N)]                 # exogenous, present at all times
  df[, y := stats::rnorm(.N)]                 # training values

  # y_t = 0.5 * y_{t-1} + 0.3 * x_t  (x same-period, from exogen)
  models <- list(
    build_model("deterministic", formula = y ~ I(0.5 * lag(y) + 0.3 * x)),
    build_model("exogen", formula = ~x)
  )
  test_start <- 15L
  horizon <- 4L
  setup <- setup_system(models, df, train_start = 1, test_start = test_start,
                        horizon = horizon, groupvar = "unit", timevar = "time",
                        inner_sims = 2)
  set.seed(2)
  res <- simulate_system(fit_system(setup, nsim = 1))

  # Closed-form expectation per unit, iterating the recurrence on the known x.
  for (u in 1:3) {
    xu <- df[unit == u][order(time)]
    y_prev <- xu[time == test_start - 1L, y]
    expected <- numeric(horizon)
    for (h in seq_len(horizon)) {
      tt <- test_start + h - 1L
      y_now <- 0.5 * y_prev + 0.3 * xu[time == tt, x]
      expected[h] <- y_now
      y_prev <- y_now
    }
    got <- unique(res[unit == u & time >= test_start, .(time, y)])
    data.table::setorder(got, time)
    expect_equal(got$y, expected, tolerance = 1e-9,
                 info = paste("unit", u))
  }
})

# ── spatial lag: W . y matches a hand-computed neighbour average ───────────

test_that("predict.spatial_lag computes the weighted neighbour average", {
  skip_if_not_installed("sfdep")
  nb  <- list(2L, c(1L, 3L), 2L)              # chain 1-2-3
  wt  <- list(1, c(0.5, 0.5), 1)              # row-standardised
  ids <- c(10, 20, 30)
  mod <- spatial_lag_model(sl ~ lag(y), nb = nb, wt = wt, unit_ids = ids)
  expect_true(mod$use_lag)

  sim_ctx <- panel_context(unit = "unit", time = "time", sim = "sim")
  grid <- data.table::CJ(unit = ids, time = 1:3, sim = 1:2)
  grid[, y := unit + time + 100 * (sim - 1)]  # sim 2 offset by 100

  out <- predict(mod, t = 2L, data = grid, ctx = sim_ctx)
  data.table::setkey(out, sim, unit)

  # use_lag => reads t - 1 = 1; result stamped at t = 2.
  expect_true(all(out$time == 2L))
  # sim 1 at t=1: y = (11, 21, 31) => sl10 = 21, sl20 = .5*11+.5*31 = 21, sl30 = 21
  expect_equal(out[sim == 1 & unit == 10, sl], 21)
  expect_equal(out[sim == 1 & unit == 20, sl], 0.5 * 11 + 0.5 * 31)
  expect_equal(out[sim == 1 & unit == 30, sl], 21)
  # per-sim independence: sim 2 is the sim-1 answer shifted by the +100 offset.
  expect_equal(out[sim == 2 & unit == 20, sl], 0.5 * 111 + 0.5 * 131)
})

test_that("spatial_lag without lag() reads the current period", {
  skip_if_not_installed("sfdep")
  nb  <- list(2L, 1L)
  wt  <- list(1, 1)
  ids <- c(1, 2)
  mod <- spatial_lag_model(sl ~ y, nb = nb, wt = wt, unit_ids = ids)
  expect_false(mod$use_lag)

  sim_ctx <- panel_context(unit = "unit", time = "time", sim = "sim")
  grid <- data.table::CJ(unit = ids, time = 1:3, sim = 1L)
  grid[, y := unit * 10 + time]
  out <- predict(mod, t = 2L, data = grid, ctx = sim_ctx)
  # reads t = 2: y1 = 12, y2 = 22 => sl1 = 22, sl2 = 12
  expect_equal(out[unit == 1, sl], 22)
  expect_equal(out[unit == 2, sl], 12)
})

# ── spatial lag: island handling ───────────────────────────────────────────

test_that("integer(0) island no longer errors and uses island_default", {
  skip_if_not_installed("sfdep")
  nb  <- list(2L, c(1L, 3L), 2L, integer(0))    # unit 4 = island (manual encoding)
  wt  <- list(1, c(0.5, 0.5), 1, numeric(0))
  ids <- c(10, 20, 30, 40)
  mod <- spatial_lag_model(sl ~ lag(y), nb = nb, wt = wt,
                           unit_ids = ids, island_default = -1)

  sim_ctx <- panel_context(unit = "unit", time = "time", sim = "sim")
  grid <- data.table::CJ(unit = ids, time = 1:3, sim = 1L)
  # t=1: y[10]=20, y[20]=30, y[30]=40, y[40]=50  (9 + unit + time)
  grid[, y := 9 + unit + time]

  # Must not error
  out <- predict(mod, t = 2L, data = grid, ctx = sim_ctx)
  data.table::setkey(out, unit)

  # Connected units: reads t-1=1; y[10]=20, y[20]=30, y[30]=40
  # sl[10] = y[20] = 30; sl[20] = 0.5*y[10]+0.5*y[30] = 30; sl[30] = y[20] = 30
  expect_equal(out[unit == 10, sl], 30)
  expect_equal(out[unit == 20, sl], 0.5 * 20 + 0.5 * 40)
  expect_equal(out[unit == 30, sl], 30)
  # Island unit 40 must get island_default, not 0 or an error
  expect_equal(out[unit == 40, sl], -1)
})

test_that("sfdep 0L island uses island_default (not silent 0)", {
  skip_if_not_installed("sfdep")
  nb  <- list(2L, c(1L, 3L), 2L, 0L)            # unit 4 = island (sfdep/spdep encoding)
  wt  <- list(1, c(0.5, 0.5), 1, 0)
  ids <- c(10, 20, 30, 40)
  mod <- spatial_lag_model(sl ~ lag(y), nb = nb, wt = wt,
                           unit_ids = ids)        # island_default = NA (default)

  sim_ctx <- panel_context(unit = "unit", time = "time", sim = "sim")
  grid <- data.table::CJ(unit = ids, time = 1:3, sim = 1L)
  grid[, y := 9 + unit + time]

  out <- predict(mod, t = 2L, data = grid, ctx = sim_ctx)
  data.table::setkey(out, unit)

  # Connected units: y[10]=20, y[20]=30, y[30]=40 at t-1=1
  expect_equal(out[unit == 10, sl], 30)
  expect_equal(out[unit == 20, sl], 0.5 * 20 + 0.5 * 40)
  expect_equal(out[unit == 30, sl], 30)
  # Island: NA_real_, not 0
  expect_true(is.na(out[unit == 40, sl]))
})

test_that("is_island mask is correct at construction for both encodings", {
  skip_if_not_installed("sfdep")
  expected_mask <- c(FALSE, FALSE, FALSE, TRUE)

  # integer(0) encoding
  mod_empty <- spatial_lag_model(
    sl ~ lag(y),
    nb = list(2L, c(1L, 3L), 2L, integer(0)),
    wt = list(1, c(0.5, 0.5), 1, numeric(0)),
    unit_ids = c(10, 20, 30, 40)
  )
  expect_equal(mod_empty$is_island, expected_mask)

  # 0L encoding
  mod_zero <- spatial_lag_model(
    sl ~ lag(y),
    nb = list(2L, c(1L, 3L), 2L, 0L),
    wt = list(1, c(0.5, 0.5), 1, 0),
    unit_ids = c(10, 20, 30, 40)
  )
  expect_equal(mod_zero$is_island, expected_mask)

  # All-connected: normalisation must be a no-op
  mod_connected <- spatial_lag_model(
    sl ~ lag(y),
    nb = list(2L, c(1L, 3L), 2L),
    wt = list(1, c(0.5, 0.5), 1),
    unit_ids = c(10, 20, 30)
  )
  expect_equal(mod_connected$is_island, c(FALSE, FALSE, FALSE))
})

# ── reproducibility & shape ────────────────────────────────────────────────

test_that("same seed + sequential plan gives identical simulations", {
  future::plan(future::sequential)
  on.exit(future::plan(future::sequential), add = TRUE)

  dt <- sim_panel_ar1(units = 4, n_time = 18, seed = 9,
                      groupvar = "unit", timevar = "time")
  models <- list(
    build_model("linear", formula = y ~ lag(y) + lag(x), boot = "resid"),
    build_model("exogen", formula = ~x)
  )
  setup <- setup_system(models, dt, train_start = 1, test_start = 14,
                        horizon = 4, groupvar = "unit", timevar = "time",
                        inner_sims = 3, min_window = 6)

  run <- function() {
    set.seed(123)
    simulate_system(fit_system(setup, nsim = 2))
  }
  a <- run()
  b <- run()
  expect_equal(a$y, b$y)
  expect_equal(a$.sim, b$.sim)
  # .sim count = nsim * inner_sims; panel attributes stamped.
  expect_equal(length(unique(a$.sim)), 2L * 3L)
  expect_identical(attr(a, "panel_unit"), "unit")
  expect_identical(attr(a, "panel_time"), "time")
})

# ── error propagation grows with horizon (gated) ───────────────────────────

test_that("AR(1) forecast variance grows with horizon and with persistence", {
  skip_on_cran()
  skip_if_not_slow()
  future::plan(future::sequential)
  on.exit(future::plan(future::sequential), add = TRUE)

  horizon_var <- function(rho) {
    dt <- sim_panel_ar1(units = 12, n_time = 60, rho = rho, beta = 0,
                        sd = 1, seed = 100, groupvar = "unit", timevar = "time")
    setup <- setup_system(
      list(build_model("linear", formula = y ~ lag(y))),
      dt, train_start = 1, test_start = 50, horizon = 8,
      groupvar = "unit", timevar = "time", inner_sims = 400)
    set.seed(3)
    res <- simulate_system(fit_system(setup, nsim = 1))
    # Predictive spread = variance across sims WITHIN a unit, averaged over
    # units (cross-unit marginal variance would swamp the horizon signal).
    per_unit <- res[time >= 50, .(v = stats::var(y)), by = c("unit", "time")]
    per_unit[, .(v = mean(v)), by = "time"][order(time), v]
  }

  v_stat <- horizon_var(0.3)
  v_pers <- horizon_var(0.9)

  # A near-unit-root system fans out strongly with the horizon; a mean-reverting
  # one stays bounded (variance plateaus quickly).
  expect_gt(v_pers[length(v_pers)], 2 * v_pers[1])        # persistent: fans out
  expect_lt(v_stat[length(v_stat)], 1.6 * v_stat[1])      # stationary: bounded
  expect_gt(v_pers[length(v_pers)], 1.5 * v_stat[length(v_stat)])
})

# ── the true model out-scores a misspecified one (gated) ───────────────────

test_that("the correctly-specified system beats an omitted-covariate one", {
  skip_on_cran()
  skip_if_not_slow()
  future::plan(future::sequential)
  on.exit(future::plan(future::sequential), add = TRUE)

  test_start <- 30L
  score_gap <- replicate_seeds(10, function(s) {
    dt <- sim_panel_ar1(units = 20, n_time = 35, rho = 0.6, beta = 0.4,
                        sd = 1, seed = s, groupvar = "unit", timevar = "time")

    fit_score <- function(models) {
      setup <- setup_system(models, dt, train_start = 1, test_start = test_start,
                            horizon = 5, groupvar = "unit", timevar = "time",
                            inner_sims = 60)
      res <- simulate_system(fit_system(setup, nsim = 1))
      acc <- get_accuracy(res[time >= test_start], "y", dt,
                          ctx = panel_context("unit", "time"))
      mean(acc$crps, na.rm = TRUE)
    }

    true_crps <- fit_score(list(
      build_model("linear", formula = y ~ lag(y) + lag(x)),
      build_model("exogen", formula = ~x)))
    wrong_crps <- fit_score(list(
      build_model("linear", formula = y ~ lag(y))))      # omits lag(x)

    true_crps - wrong_crps          # negative => true model is better
  })

  # Aggregate ordering: the true model has lower CRPS on average, and wins in a
  # clear majority of replications (indicative, not per-cell).
  expect_lt(mean(score_gap), 0)
  expect_gt(mean(score_gap < 0), 0.7)
})

# ── long-horizon vs dynamic simulator parity (gated) ───────────────────────

test_that("long-horizon and dynamic accuracy stack and are comparable", {
  skip_on_cran()
  skip_if_not_slow()
  future::plan(future::sequential)
  on.exit(future::plan(future::sequential), add = TRUE)

  dt <- sim_panel_ar1(units = 20, n_time = 40, rho = 0.6, beta = 0.4,
                      sd = 1, seed = 7, groupvar = "unit", timevar = "time")
  test_start <- 33L
  horizons <- 1:4

  # Dynamic simulator, scored per (unit, horizon).
  setup <- setup_system(
    list(build_model("linear", formula = y ~ lag(y) + lag(x)),
         build_model("exogen", formula = ~x)),
    dt, train_start = 1, test_start = test_start, horizon = max(horizons),
    groupvar = "unit", timevar = "time", inner_sims = 80)
  set.seed(4)
  sim_res <- simulate_system(fit_system(setup, nsim = 1))
  sim_acc <- get_accuracy(sim_res[time >= test_start], "y", dt,
                          ctx = panel_context("unit", "time"),
                          test_start = test_start, by = c("unit", "horizon"))

  # Direct long-horizon reduced form, conditioned on the same origin info.
  lh_setup <- setup_long_horizon(
    dt, list(lh = lead_horizon(y) ~ lag(y) + lag(x)),
    horizons, "unit", "time", test_start = test_start)
  lh_fc <- forecast_long_horizon(lh_setup, dt, nsim = 1, inner_sims = 80)
  lh_acc <- get_lh_accuracy(lh_fc, dt, lh_setup)

  stacked <- compare_approaches(lh_acc, sim_acc)
  expect_true("approach" %in% names(stacked))
  expect_setequal(unique(stacked$approach), c("long_horizon", "simulation"))
  expect_true(all(is.finite(stacked$crps)))

  # Comparable accuracy on a correctly-specified DGP (broad band, indicative).
  ratio <- mean(sim_acc$crps, na.rm = TRUE) / mean(lh_acc$crps, na.rm = TRUE)
  expect_gt(ratio, 0.4)
  expect_lt(ratio, 2.5)
})
