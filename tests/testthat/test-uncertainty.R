# Pillar 3 — uncertainty estimation -----------------------------------------
#
# Pins the predictive-draw shapes and the scoring kernel (deterministic), and
# MEASURES the calibration / independence concerns flagged in the roadmap
# (gated, indicative under broken assumptions).

# ── .score_draws kernel ────────────────────────────────────────────────────

test_that(".score_draws returns NA for NA truth or all-NA draws", {
  expect_true(is.na(.score_draws(list(stats::rnorm(50)), NA_real_)$crps))
  expect_true(is.na(.score_draws(list(rep(NA_real_, 50)), 1)$crps))
})

test_that(".score_draws CRPS is non-negative and equals MAE for a point mass", {
  set.seed(1)
  s <- .score_draws(list(stats::rnorm(500), stats::rnorm(500) + 2),
                    c(0.2, 1.5))
  expect_true(all(s$crps >= 0))

  # A degenerate point-mass forecast: CRPS == MAE == |point - truth|.
  pm <- .score_draws(list(rep(3, 100)), 5)
  expect_equal(pm$crps, pm$mae)
  expect_equal(pm$crps, 2)
})

test_that(".score_draws Winkler penalty triggers only outside the interval", {
  set.seed(1)
  d <- stats::rnorm(5000)
  width <- unname(diff(stats::quantile(d, c(0.05, 0.95))))
  inside  <- .score_draws(list(d), 0, level = 90)$winkler
  outside <- .score_draws(list(d), 10, level = 90)$winkler
  expect_equal(inside, width, tolerance = 1e-8)   # inside => width, no penalty
  expect_gt(outside, width)                        # outside => width + penalty
})

# ── getpi adds residual scale (PI strictly wider than the CI) ──────────────

test_that("getpi widens the parameter-only interval by the residual scale", {
  set.seed(1)
  dt <- data.frame(x = stats::rnorm(300))
  dt$y <- 1 + 2 * dt$x + stats::rnorm(300, sd = 2)
  fit <- stats::lm(y ~ x, dt)
  p <- stats::predict(fit, dt[1, , drop = FALSE], se.fit = TRUE)

  set.seed(2)
  pi <- getpi(p, nsamples = 5000)                          # parameter + residual
  set.seed(2)
  ci <- p$fit + p$se.fit * stats::rt(5000, p$df)           # parameter only
  pi_w <- diff(stats::quantile(pi, c(0.05, 0.95)))
  ci_w <- diff(stats::quantile(ci, c(0.05, 0.95)))
  expect_gt(pi_w, 3 * ci_w)              # residual scale dominates here
})

# ── coverage / calibration (gated; indicative) ─────────────────────────────

test_that("a calibrated linear baseline has ~nominal one-step coverage", {
  skip_on_cran()
  skip_if_not_slow()
  future::plan(future::sequential)
  on.exit(future::plan(future::sequential), add = TRUE)

  # DGP that SATISFIES the per-equation assumptions: independent rows, correct
  # functional form, homoscedastic. One-step forecast => coverage ~ nominal.
  dt <- sim_panel_pooled(units = 150, n_time = 15, b0 = 1, b1 = 2, b2 = -1,
                         sd = 1, seed = 1, groupvar = "unit", timevar = "time")
  models <- list(
    build_model("linear", formula = y ~ x1 + x2),
    build_model("exogen", formula = ~x1),
    build_model("exogen", formula = ~x2)
  )
  setup <- setup_system(models, dt, train_start = 1, test_start = 15,
                        horizon = 1, groupvar = "unit", timevar = "time",
                        inner_sims = 600)
  set.seed(2)
  res <- simulate_system(fit_system(setup, nsim = 1))

  fc <- res[time == 15]
  draws <- split(fc$y, fc$unit)
  truth_dt <- dt[time == 15][order(unit)]
  truth <- truth_dt$y
  draws <- draws[as.character(truth_dt$unit)]

  cov90 <- empirical_coverage(draws, truth, level = 90)
  expect_gt(cov90, 0.83)
  expect_lt(cov90, 0.97)

  # PIT ~ uniform: mean near 0.5, mass not piled in the tails.
  pit <- pit_values(draws, truth)
  expect_equal(mean(pit), 0.5, tolerance = 0.08)
})

test_that("boot + random window widens the PI vs the plain baseline", {
  skip_on_cran()
  skip_if_not_slow()
  future::plan(future::sequential)
  on.exit(future::plan(future::sequential), add = TRUE)

  dt <- sim_panel_pooled(units = 80, n_time = 16, b0 = 1, b1 = 2, b2 = -1,
                         sd = 1, seed = 3, groupvar = "unit", timevar = "time")
  mk <- function(boot, min_window, nsim) {
    models <- list(
      build_model("linear", formula = y ~ x1 + x2, boot = boot),
      build_model("exogen", formula = ~x1),
      build_model("exogen", formula = ~x2)
    )
    setup <- setup_system(models, dt, train_start = 1, test_start = 16,
                          horizon = 1, groupvar = "unit", timevar = "time",
                          inner_sims = 300, min_window = min_window)
    set.seed(4)
    res <- simulate_system(fit_system(setup, nsim = nsim))
    fc <- res[time == 16]
    mean_interval_width(split(fc$y, fc$unit), level = 90)
  }

  base_w <- mk(boot = NULL, min_window = NULL, nsim = 1)
  over_w <- mk(boot = "resid", min_window = 8, nsim = 30)
  # The random window + residual bootstrap stack extra parameter variance on
  # top of the residual scale -> measurably wider (documented over-dispersion).
  expect_gt(over_w, base_w)
})

# ── parametric_distribution draws match the fitted family ──────────────────

test_that("parametric_distribution draws match the fitted normal moments", {
  skip_on_cran()
  skip_if_not_slow()
  skip_if_not_installed("fitdistrplus")

  set.seed(1)
  dt <- data.table::as.data.table(expand.grid(unit = 1:5, time = 1:40))
  dt[, v := stats::rnorm(.N, mean = 3, sd = 2)]
  ctx <- panel_context(unit = "unit", time = "time", sim = "sim")
  m <- parametric_distribution_model(~v, distribution = "norm", data = dt,
                                     ctx = ctx)

  grid <- data.table::CJ(unit = 1:5, time = 41:44, sim = 1:400)
  grid[, v := NA_real_]
  set.seed(2)
  out <- predict(m, data = grid, ctx = ctx, test_start = 41, horizon = 4,
                 inner_sims = 400)
  fitted_mu <- unname(m$fitted$estimate["mean"])
  fitted_sd <- unname(m$fitted$estimate["sd"])
  expect_equal(mean(out$v), fitted_mu, tolerance = 0.1)
  expect_equal(stats::sd(out$v), fitted_sd, tolerance = 0.1)
})

# ── innovation independence: the documented SUR / common-shock gap ─────────

test_that("simulated innovations are cross-unit independent (SUR gap, documented)", {
  skip_on_cran()
  skip_if_not_slow()
  future::plan(future::sequential)
  on.exit(future::plan(future::sequential), add = TRUE)

  # DGP with a shock SHARED across units each period -> the truth has strong
  # cross-unit correlation. The simulator draws innovations independently per
  # (unit, sim) row, so it CANNOT reproduce that correlation. This test
  # documents the gap (it is an expected shortfall, not an engine bug).
  dt <- sim_panel_common_shock(units = 8, n_time = 40, rho = 0.4,
                               common_sd = 1.5, idio_sd = 0.3, seed = 1,
                               groupvar = "unit", timevar = "time")

  # Truth: detrended y is strongly cross-unit correlated (the common shock).
  wide <- data.table::dcast(dt, time ~ unit, value.var = "y")
  ymat <- as.matrix(wide[, -1])
  truth_cor <- mean(stats::cor(diff(ymat))[upper.tri(diag(ncol(ymat)))])
  expect_gt(truth_cor, 0.3)

  # Simulated ensemble: per-(unit) deviation from the ensemble mean at a fixed
  # forecast time, correlated ACROSS units over sims.
  setup <- setup_system(
    list(build_model("linear", formula = y ~ lag(y))),
    dt, train_start = 1, test_start = 35, horizon = 3,
    groupvar = "unit", timevar = "time", inner_sims = 400)
  set.seed(2)
  res <- simulate_system(fit_system(setup, nsim = 1))

  cell <- res[time == 36]
  cell[, dev := y - mean(y), by = "unit"]
  dev_wide <- data.table::dcast(cell, .sim ~ unit, value.var = "dev")
  dmat <- as.matrix(dev_wide[, -1])
  sim_cor <- mean(stats::cor(dmat)[upper.tri(diag(ncol(dmat)))])

  # Near zero: the engine's innovations are independent by construction.
  expect_lt(abs(sim_cor), 0.15)
})
