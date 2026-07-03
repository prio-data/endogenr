# Tests for the unified predictive-draw interface ----------------------------
#
# draw_predictive(model, newdata, n_param, n_innov) returns a numeric matrix,
# nrow(newdata) x (max(n_param, 1) * max(n_innov, 1)); column (j - 1) * K + k
# holds parameter draw j, innovation draw k. These tests pin the shape
# contract, the exact expectation equivalences (n_param = 0, n_innov = 0), the
# parameter/innovation variance separation, the warn-once behaviour of
# innovation-only families, and the fable path-coherence contract.

# Pooled cross-sectional fixture with unit/time/sim keys, a gaussian outcome
# (y, residual sd = 2) and a poisson count outcome (ycount).
draw_data <- function(seed = 1, units = 25L, n_time = 24L) {
  dt <- sim_panel_pooled(units = units, n_time = n_time, b0 = 1, b1 = 2,
                         b2 = -1, sd = 2, seed = seed,
                         groupvar = "unit", timevar = "time")
  dt[, sim := 1L]
  dt[, ycount := stats::rpois(.N, lambda = exp(0.4 + 0.3 * x1))]
  dt[]
}

draw_ctx <- function() panel_context(unit = "unit", time = "time", sim = "sim")

# ── shape contract ──────────────────────────────────────────────────────────

test_that("draw_predictive returns the nrow(newdata) x (P * K) matrix contract", {
  dt  <- draw_data()
  ctx <- draw_ctx()
  nd  <- dt[time == 1L]
  lin   <- linearmodel(y ~ x1, data = dt, ctx = ctx)
  glm_m <- glmmodel(ycount ~ x1, family = stats::poisson(), data = dt, ctx = ctx)

  set.seed(1)
  for (m in list(lin, glm_m)) {
    d32 <- draw_predictive(m, nd, n_param = 3, n_innov = 2)
    expect_true(is.matrix(d32) && is.numeric(d32))
    expect_equal(dim(d32), c(nrow(nd), 6L))

    d11 <- draw_predictive(m, nd)                    # defaults: one 1 x 1 block
    expect_true(is.matrix(d11))
    expect_equal(dim(d11), c(nrow(nd), 1L))

    d1x1 <- draw_predictive(m, nd[1L])               # a matrix even at 1 x 1
    expect_true(is.matrix(d1x1))
    expect_equal(dim(d1x1), c(1L, 1L))

    d05 <- draw_predictive(m, nd, n_param = 0, n_innov = 5)
    expect_equal(dim(d05), c(nrow(nd), 5L))

    expect_error(draw_predictive(m, nd, n_param = -1),
                 "non-negative finite number")
    expect_error(draw_predictive(m, nd, n_innov = NA_real_),
                 "non-negative finite number")
    expect_error(draw_predictive(m, nd, n_param = c(1, 2)),
                 "non-negative finite number")
  }
})

test_that(".lm_predictive_draws lays out columns as (param block j) x (innovation k)", {
  # A degenerate residual scale isolates the parameter component: within a
  # parameter block all K columns collapse onto that block's drawn mu, so the
  # (j - 1) * K + k layout is directly visible.
  p <- list(fit = 0, se.fit = 1, residual.scale = 1e-9, df = 1e6)
  set.seed(11)
  d <- .lm_predictive_draws(p, n_param = 2, n_innov = 3, param_scope = "draw")
  expect_equal(dim(d), c(1L, 6L))
  expect_lt(max(abs(d[1, 1:3] - d[1, 1])), 1e-6)   # block 1 shares one mu draw
  expect_lt(max(abs(d[1, 4:6] - d[1, 4])), 1e-6)   # block 2 likewise
  expect_gt(abs(d[1, 4] - d[1, 1]), 1e-3)          # blocks are distinct draws
})

# ── expectation equivalences (exact) ────────────────────────────────────────

test_that("(n_param = 0, n_innov = 0) is exactly the deterministic conditional mean", {
  dt  <- draw_data()
  ctx <- draw_ctx()
  nd  <- dt[time <= 2L]
  lin   <- linearmodel(y ~ x1, data = dt, ctx = ctx)
  glm_m <- glmmodel(ycount ~ x1, family = stats::poisson(), data = dt, ctx = ctx)

  lin_mean <- draw_predictive(lin, nd, n_param = 0, n_innov = 0)
  expect_equal(dim(lin_mean), c(nrow(nd), 1L))
  expect_equal(as.vector(lin_mean),
               unname(stats::predict(lin$fitted, newdata = nd, se.fit = TRUE)$fit))

  glm_mean <- draw_predictive(glm_m, nd, n_param = 0, n_innov = 0)
  link <- stats::predict(glm_m$fitted, newdata = nd, type = "link", se.fit = TRUE)
  expect_equal(as.vector(glm_mean),
               unname(glm_m$family$linkinv(link$fit)))
})

# ── parameter vs innovation separation, linear (gated) ──────────────────────

test_that("linear draws separate parameter and innovation uncertainty", {
  skip_on_cran()
  skip_if_not_slow()

  dt  <- draw_data()
  ctx <- draw_ctx()
  lin <- linearmodel(y ~ x1, data = dt, ctx = ctx)
  nd1 <- dt[1L]
  p   <- stats::predict(lin$fitted, newdata = nd1, se.fit = TRUE)
  t_sd <- sqrt(p$df / (p$df - 2))                  # sd of a standard t_df

  # Repeated engine-scope draws through the public generic: the per-row
  # marginal is fit + sqrt(se.fit^2 + scale^2) * t_df.
  set.seed(31)
  full <- replicate(3000, as.vector(
    draw_predictive(lin, nd1, n_param = 1, n_innov = 1, param_scope = "row")))
  expect_equal(stats::sd(full),
               sqrt(p$se.fit^2 + p$residual.scale^2) * t_sd,
               tolerance = 0.1)

  # Parameter-only draws (n_innov = 0) spread on the se.fit scale ...
  set.seed(32)
  par_only <- as.vector(draw_predictive(lin, nd1, n_param = 400, n_innov = 0,
                                        param_scope = "draw"))
  expect_equal(stats::sd(par_only), p$se.fit * t_sd, tolerance = 0.15)

  # ... innovation-only draws (n_param = 0) on the residual scale, far wider.
  set.seed(33)
  innov_only <- as.vector(draw_predictive(lin, nd1, n_param = 0, n_innov = 400))
  expect_equal(stats::sd(innov_only), p$residual.scale, tolerance = 0.15)
  expect_gt(stats::sd(innov_only), 3 * stats::sd(par_only))
})

# ── innovation-only families warn once that n_param is ignored ──────────────

test_that("heterolm draws warn once that n_param is ignored, then stay silent", {
  skip_if_not_installed("heterolm")
  dt  <- draw_data()
  ctx <- draw_ctx()
  m   <- heterolmmodel(y ~ x1, variance = ~1, data = dt, ctx = ctx)
  nd  <- dt[time == 1L]

  options(.endogenr_no_param_draw_heterolm = NULL)
  set.seed(41)
  expect_warning(
    d1 <- draw_predictive(m, nd, n_param = 1, n_innov = 1),
    "carry no parameter-uncertainty draw")
  # one-time-per-session: a second call is silent
  expect_silent(d2 <- draw_predictive(m, nd, n_param = 1, n_innov = 1))
  expect_equal(dim(d1), c(nrow(nd), 1L))
  expect_true(all(is.finite(d1)) && all(is.finite(d2)))
})

test_that("gamlss draws warn once that n_param is ignored, then stay silent", {
  skip_if_no_gamlss()
  dt  <- draw_data()
  ctx <- draw_ctx()
  m   <- gamlssmodel(y ~ x1, data = dt, ctx = ctx)
  nd  <- dt[time == 1L]

  options(.endogenr_no_param_draw_gamlss = NULL)
  set.seed(42)
  expect_warning(
    d1 <- draw_predictive(m, nd, n_param = 1, n_innov = 1),
    "carry no parameter-uncertainty draw")
  expect_silent(d2 <- draw_predictive(m, nd, n_param = 1, n_innov = 1))
  expect_equal(dim(d1), c(nrow(nd), 1L))
  expect_true(all(is.finite(d1)) && all(is.finite(d2)))
})

# ── default method refuses non-stochastic families ──────────────────────────

test_that("draw_predictive.default rejects models without a predictive distribution", {
  ctx <- draw_ctx()
  m <- deterministicmodel(w ~ I(lag(v)), ctx = ctx)
  expect_error(draw_predictive(m, data.frame(v = 1:3)),
               "no stochastic predictive distribution")
})

# ── parametric_distribution draws ───────────────────────────────────────────

test_that("parametric_distribution rejects n_innov = 0 and draws from the fit", {
  skip_if_not_installed("fitdistrplus")
  set.seed(61)
  dt <- data.table::data.table(unit = 1L, time = 1:500,
                               v = stats::rnorm(500, mean = 1, sd = 2))
  m <- parametric_distribution_model(~v, distribution = "norm", data = dt,
                                     ctx = draw_ctx())

  nd <- data.table::data.table(row = 1:50)
  expect_error(draw_predictive(m, nd, n_param = 0, n_innov = 0),
               "not available")
  set.seed(62)
  d <- draw_predictive(m, nd, n_param = 0, n_innov = 3)
  expect_equal(dim(d), c(50L, 3L))
  expect_true(all(is.finite(d)))
})

test_that("parametric_distribution parameter draws dominate the column-mean spread", {
  skip_on_cran()
  skip_if_not_slow()
  skip_if_not_installed("fitdistrplus")

  set.seed(63)
  dt <- data.table::data.table(unit = 1L, time = 1:500,
                               v = stats::rnorm(500, mean = 1, sd = 2))
  m <- parametric_distribution_model(~v, distribution = "norm", data = dt,
                                     ctx = draw_ctx())
  nd <- data.table::data.table(row = seq_len(20000))

  # Each n_param column shares ONE MVN(estimate, vcov) parameter draw, so its
  # 20000-row mean moves with the drawn mu (sd ~ 2 / sqrt(500) ~ 0.089). At
  # the point estimates a column mean carries only innovation noise
  # (~ 2 / sqrt(20000) ~ 0.014); n_innov = 60 in one call is exactly 60
  # independent point-estimate innovation columns.
  set.seed(64)
  d_param <- draw_predictive(m, nd, n_param = 60, n_innov = 1)
  set.seed(65)
  d_point <- draw_predictive(m, nd, n_param = 0, n_innov = 60)
  expect_equal(dim(d_param), c(20000L, 60L))
  expect_equal(dim(d_point), c(20000L, 60L))
  expect_gt(stats::sd(colMeans(d_param)), 3 * stats::sd(colMeans(d_point)))
})

# ── univariate_fable: coherent generate() paths mapped onto newdata ─────────

# Strongly-trending 2-unit panel; a fitted ETS(A,A,N) yields smooth coherent
# paths (same construction as test-univariate-fable.R).
fable_draw_panel <- function(units = 2L, t_n = 40L, seed = 71) {
  set.seed(seed)
  rows <- lapply(seq_len(units), function(u) {
    data.table::data.table(unit = u, time = seq_len(t_n),
                           y = cumsum(0.8 + stats::rnorm(t_n, sd = 0.4)))
  })
  data.table::rbindlist(rows)
}

test_that("univariate_fable draws are generate() paths joined onto newdata", {
  skip_if_no_fable()
  dt  <- fable_draw_panel()
  ctx <- draw_ctx()
  m <- univariate_fable_model(y ~ error("A") + trend("A") + season("N"),
                              data = dt, method = "ets", ctx = ctx)
  nd <- data.table::CJ(unit = 1:2, time = 41:50)

  set.seed(72)
  d <- draw_predictive(m, nd, n_param = 0, n_innov = 4, ctx = ctx, horizon = 10)
  expect_true(is.matrix(d))
  expect_equal(dim(d), c(nrow(nd), 4L))
  expect_true(all(is.finite(d)))

  # The columns ARE fabletools::generate() paths under the same seed, joined
  # onto the newdata rows by (unit, time).
  set.seed(72)
  manual <- data.table::as.data.table(dplyr::as_tibble(
    fabletools::generate(m$fitted, h = 10, times = 4)))
  for (k in 1:4) {
    exp_k <- manual[.rep == as.character(k)][nd, on = c("unit", "time")]$.sim
    expect_equal(d[, k], exp_k)
  }

  # Rows outside the generated h-step window come back NA.
  nd_out <- rbind(nd, data.table::data.table(unit = 1L, time = 60L))
  set.seed(73)
  d_out <- draw_predictive(m, nd_out, n_param = 0, n_innov = 2, ctx = ctx,
                           horizon = 10)
  expect_true(all(is.na(d_out[nrow(nd_out), ])))
  expect_true(all(is.finite(d_out[-nrow(nd_out), ])))

  # Refusals: no conditional mean; ctx is mandatory.
  expect_error(draw_predictive(m, nd, n_param = 0, n_innov = 0, ctx = ctx),
               "not available for univariate_fable")
  expect_error(draw_predictive(m, nd, n_param = 0, n_innov = 1),
               "ctx is required")
})

test_that("univariate_fable draw columns are temporally coherent paths", {
  skip_if_no_fable()
  dt  <- fable_draw_panel()
  ctx <- draw_ctx()
  m <- univariate_fable_model(y ~ error("A") + trend("A") + season("N"),
                              data = dt, method = "ets", ctx = ctx)
  nd <- data.table::CJ(unit = 1:2, time = 41:50)

  set.seed(74)
  d <- draw_predictive(m, nd, n_param = 0, n_innov = 40, ctx = ctx, horizon = 10)

  # Deviations from the per-(unit, time) ensemble mean remove the shared
  # deterministic trend; a coherent ETS path integrates its innovations, so
  # within-column deviations stay strongly lag-1 autocorrelated (independent
  # per-horizon stitching would give ~0 — see test-univariate-fable.R).
  dev <- d - rowMeans(d)
  acf1 <- function(v) stats::cor(v[-length(v)], v[-1])
  per_path <- sapply(1:2, function(u) apply(dev[nd$unit == u, ], 2, acf1))
  expect_gt(mean(per_path), 0.4)
})
