# Pillar 1 — core data flow & scheduling ------------------------------------
#
# The load-bearing invariants of the simulation grid and the dynamic loop:
# grid shape, the (unit, sim, time) sort key the positional lag depends on,
# independent-then-dependent population, and per-(unit, sim) trajectory
# isolation (no cross-unit / cross-sim leakage).

# ── grid integrity ─────────────────────────────────────────────────────────

test_that("setup grid is units x times x inner_sims with NA-only forecast cells", {
  dt <- sim_panel_ar1(units = 4, n_time = 20, seed = 1,
                      groupvar = "unit", timevar = "time")
  models <- list(
    build_model("linear", formula = y ~ lag(x)),
    build_model("exogen", formula = ~x)
  )
  setup <- setup_system(models, dt, train_start = 1, test_start = 15,
                        horizon = 3, groupvar = "unit", timevar = "time",
                        inner_sims = 5)
  sd <- setup$simulation_data

  # 4 units x (1..17) times x 5 sims
  expect_equal(nrow(sd), 4L * 17L * 5L)
  expect_setequal(unique(sd$sim), 1:5)
  expect_equal(range(sd$time), c(1L, 17L))

  # training cells carry y; forecast cells (time >= test_start) are NA — and
  # this holds for every column (the grid is filled from training only; the
  # exogenous x forecast cells are populated later, during simulation).
  expect_false(anyNA(sd[time < 15, y]))
  expect_true(all(is.na(sd[time >= 15, y])))
  expect_false(anyNA(sd[time < 15, x]))
  expect_true(all(is.na(sd[time >= 15, x])))
})

test_that("simulation grid is sorted by (unit, sim, time)", {
  dt <- sim_panel_ar1(units = 3, n_time = 12, seed = 2,
                      groupvar = "unit", timevar = "time")
  models <- list(
    build_model("linear", formula = y ~ lag(x)),
    build_model("exogen", formula = ~x)
  )
  setup <- setup_system(models, dt, train_start = 1, test_start = 9,
                        horizon = 2, groupvar = "unit", timevar = "time",
                        inner_sims = 3)
  sd <- setup$simulation_data
  # The key the positional lag relies on.
  expect_equal(data.table::key(sd), c("unit", "sim", "time"))
  # Within every (unit, sim) block, time is strictly increasing.
  chk <- sd[, .(ok = !is.unsorted(time, strictly = TRUE)), by = c("unit", "sim")]
  expect_true(all(chk$ok))
})

test_that("execution order respects a same-period x -> y -> z chain", {
  # SAME-PERIOD dependencies force an ordering; purely-lagged chains do not
  # (every term reads an already-computed prior period, so order is free).
  dt <- sim_panel_ar1(units = 3, n_time = 16, seed = 3,
                      groupvar = "unit", timevar = "time")
  dt[, z := 0.3 * y + stats::rnorm(.N, sd = 0.1)]
  models <- list(
    build_model("linear", formula = z ~ y),       # same-period y
    build_model("linear", formula = y ~ x),        # same-period x
    build_model("exogen", formula = ~x)
  )
  setup <- setup_system(models, dt, train_start = 1, test_start = 12,
                        horizon = 2, groupvar = "unit", timevar = "time",
                        inner_sims = 1)
  ord <- setup$execution_order
  expect_true(which(ord == "x") < which(ord == "y"))
  expect_true(which(ord == "y") < which(ord == "z"))
})

# ── independent-then-dependent population ──────────────────────────────────

test_that("an independent output is populated before a dependent model reads it", {
  future::plan(future::sequential)
  on.exit(future::plan(future::sequential), add = TRUE)

  set.seed(10)
  df <- data.table::as.data.table(expand.grid(unit = 1:3, time = 1:20))
  df[, v := stats::rnorm(.N, mean = 5)]
  df[, w := v]                      # derived output; needs an initial state

  models <- list(
    build_model("parametric_distribution", formula = ~v, distribution = "norm"),
    build_model("deterministic", formula = w ~ I(lag(v)))
  )
  setup <- setup_system(models, df, train_start = 1, test_start = 15,
                        horizon = 4, groupvar = "unit", timevar = "time",
                        inner_sims = 4)
  set.seed(11)
  res <- simulate_system(fit_system(setup, nsim = 2))
  fc <- res[time >= 15]

  # v (independent) is fully populated across the whole forecast horizon, and
  # w (dependent) = lag(v) is therefore never NA.
  expect_false(anyNA(fc$v))
  expect_false(anyNA(fc$w))

  # w at t equals v at t-1 within the same (unit, sim) trajectory.
  data.table::setkey(res, unit, .sim, time)
  res[, v_lag := data.table::shift(v, 1L), by = c("unit", ".sim")]
  chk <- res[time >= 15]
  expect_equal(chk$w, chk$v_lag)
})

# ── trajectory isolation: no cross-unit / cross-sim leakage ────────────────

test_that("a deterministic recurrence stays within its own (unit, sim) series", {
  future::plan(future::sequential)
  on.exit(future::plan(future::sequential), add = TRUE)

  # Distinct per-unit levels; a pure deterministic growth recurrence.
  df <- data.table::as.data.table(expand.grid(unit = 1:4, time = 1:20))
  df[, y := unit * 100 + time]      # unit-specific, so any leak is detectable
  models <- list(build_model("deterministic", formula = y ~ I(lag(y) * 1.1)))

  setup <- setup_system(models, df, train_start = 1, test_start = 15,
                        horizon = 4, groupvar = "unit", timevar = "time",
                        inner_sims = 3)
  set.seed(1)
  res <- simulate_system(fit_system(setup, nsim = 1))

  # Deterministic => every sim of a unit is identical (no cross-sim noise).
  spread <- res[time >= 15, .(rng = diff(range(y))), by = c("unit", "time")]
  expect_true(all(spread$rng < 1e-9))

  # Each unit follows ITS OWN recurrence y_t = y_origin * 1.1^(t - origin).
  origin <- 14L
  for (u in 1:4) {
    y0 <- u * 100 + origin
    got <- unique(res[unit == u & time >= 15, .(time, y)])
    data.table::setorder(got, time)
    expect_equal(got$y, y0 * 1.1^(got$time - origin), tolerance = 1e-8)
  }
})

test_that("stochastic trajectories diverge across sims but share the grid keys", {
  future::plan(future::sequential)
  on.exit(future::plan(future::sequential), add = TRUE)

  dt <- sim_panel_ar1(units = 3, n_time = 20, rho = 0.7, seed = 5,
                      groupvar = "unit", timevar = "time")
  models <- list(
    build_model("linear", formula = y ~ lag(y)),
    build_model("exogen", formula = ~x)
  )
  setup <- setup_system(models, dt, train_start = 1, test_start = 15,
                        horizon = 4, groupvar = "unit", timevar = "time",
                        inner_sims = 6)
  set.seed(2)
  res <- simulate_system(fit_system(setup, nsim = 2))

  # .sim count = nsim * inner_sims
  expect_equal(length(unique(res$.sim)), 2L * 6L)

  # At a forecast cell, different sims give different draws (per-sim innovation).
  cell <- res[unit == 1 & time == 16]
  expect_gt(stats::sd(cell$y), 0)
})
