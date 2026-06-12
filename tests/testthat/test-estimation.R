# Pillar 2 — model estimation -----------------------------------------------
#
# Window selection, bootstrap mechanics, per-unit materialisation, the
# get_coefficients contract, and (gated) parameter recovery on known DGPs.

# ── get_train_window ───────────────────────────────────────────────────────

test_that("get_train_window returns the full range when min_window is NULL", {
  w <- get_train_window(1990, 2010, min_window = NULL)
  expect_equal(w$start, 1990)
  expect_equal(w$end, 2009)          # test_start - 1
  expect_null(w$window)
})

test_that("get_train_window honours min_window bounds within the training span", {
  set.seed(1)
  for (i in 1:300) {
    w <- get_train_window(1990, 2010, min_window = 8)
    expect_gte(w$start, 1990)
    expect_lte(w$end, 2009)                  # always <= test_start - 1
    expect_gte(w$end - w$start + 1, 8)       # at least min_window long
  }
})

test_that("get_train_window returns the full window when min_window equals the range", {
  w <- get_train_window(1990, 2010, min_window = 20)
  expect_equal(w$start, 1990)
  expect_equal(w$end, 2009)
})

test_that("get_train_window rejects an impossible min_window", {
  expect_error(get_train_window(1990, 2010, min_window = 50), "smaller or equal")
  expect_error(get_train_window(1990, 2010, min_window = 0), "1 or larger")
})

# ── bootstrap mechanics ────────────────────────────────────────────────────

test_that("bootstraplm recovers the OLS estimate on the response scale", {
  set.seed(1)
  dt <- data.table::data.table(x = stats::rnorm(400))
  dt[, y := 1 + 2 * x + stats::rnorm(.N, sd = 1.5)]
  ols <- stats::lm(y ~ x, dt)

  set.seed(2)
  bx_resid <- replicate(300, stats::coef(bootstraplm(y ~ x, dt, "resid"))["x"])
  bx_wild  <- replicate(300, stats::coef(bootstraplm(y ~ x, dt, "wild"))["x"])

  expect_equal(mean(bx_resid), unname(stats::coef(ols)["x"]), tolerance = 0.05)
  expect_equal(mean(bx_wild),  unname(stats::coef(ols)["x"]), tolerance = 0.05)
  expect_true(all(is.finite(bx_resid)) && all(is.finite(bx_wild)))
})

test_that("bootstrapglm reconstructs y on the link scale and recovers coefs", {
  set.seed(1)
  dt <- data.table::data.table(x = stats::rnorm(400))
  dt[, y := 0.5 + 1.2 * x + stats::rnorm(.N, sd = 1)]
  g <- stats::glm(y ~ x, dt, family = stats::gaussian())

  set.seed(2)
  bx <- replicate(300, stats::coef(bootstrapglm(y ~ x, dt, stats::gaussian(), "resid"))["x"])
  expect_equal(mean(bx), unname(stats::coef(g)["x"]), tolerance = 0.05)

  # log-link family: reconstruction is on the link scale (recovers the slope).
  dt[, yc := stats::rpois(.N, exp(0.3 + 0.4 * x))]
  gp <- stats::glm(yc ~ x, dt, family = stats::poisson())
  set.seed(3)
  bxp <- suppressWarnings(
    replicate(200, stats::coef(bootstrapglm(yc ~ x, dt, stats::poisson(), "resid"))["x"]))
  expect_equal(mean(bxp), unname(stats::coef(gp)["x"]), tolerance = 0.05)
})

# ── per-unit, time-ordered materialisation (predict path) ──────────────────

test_that("materialize_formula lags within (unit, sim), with no cross-group bleed", {
  ctx <- panel_context(unit = "unit", time = "time", sim = "sim")
  grid <- data.table::CJ(unit = c(1, 2), sim = c(1, 2), time = 1:3)
  # distinct value per (unit, sim, time) so any bleed is detectable
  grid[, x := unit * 1000 + sim * 100 + time]
  grid[, y := 0]

  mat <- materialize_formula(y ~ lag(x), grid, ctx)
  lagcol <- setdiff(names(mat), c("unit", "sim", "time", "y"))
  data.table::setkey(mat, unit, sim, time)

  # first time step of every (unit, sim) is NA; later steps read t-1 of the
  # SAME (unit, sim) — never another group's value.
  expect_true(all(is.na(mat[time == 1, get(lagcol)])))
  expect_equal(mat[time == 2, get(lagcol)], grid[time == 1, x])
  expect_equal(mat[time == 3, get(lagcol)], grid[time == 2, x])
})

# ── get_coefficients: shared vs refit ──────────────────────────────────────

test_that("get_coefficients: a shared fit is identical across draws (full window)", {
  future::plan(future::sequential)
  on.exit(future::plan(future::sequential), add = TRUE)
  dt <- sim_panel_ar1(units = 6, n_time = 30, seed = 1,
                      groupvar = "unit", timevar = "time")
  setup <- setup_system(
    list(build_model("linear", formula = y ~ lag(x)),
         build_model("exogen", formula = ~x)),
    dt, train_start = 1, test_start = 25, horizon = 3,
    groupvar = "unit", timevar = "time", inner_sims = 1)   # no min_window

  set.seed(2)
  co <- get_coefficients(fit_system(setup, nsim = 3))
  xrows <- co[term == "lag_x"]
  # one shared fit => identical estimate across the 3 draws
  expect_equal(length(unique(xrows$estimate)), 1L)
  # window is the full training range
  expect_true(all(xrows$.window_start == 1L & xrows$.window_end == 24L))
})

test_that("get_coefficients: random-window refit varies and records its window", {
  future::plan(future::sequential)
  on.exit(future::plan(future::sequential), add = TRUE)
  dt <- sim_panel_ar1(units = 6, n_time = 30, seed = 1,
                      groupvar = "unit", timevar = "time")
  setup <- setup_system(
    list(build_model("linear", formula = y ~ lag(x), boot = "resid"),
         build_model("exogen", formula = ~x)),
    dt, train_start = 1, test_start = 25, horizon = 3,
    groupvar = "unit", timevar = "time", inner_sims = 1, min_window = 8)

  set.seed(2)
  co <- get_coefficients(fit_system(setup, nsim = 5))
  xrows <- co[term == "lag_x"]
  # refit per draw => estimates vary
  expect_gt(length(unique(xrows$estimate)), 1L)
  # recorded windows are within [train_start, test_start - 1] and >= min_window
  expect_true(all(xrows$.window_start >= 1L))
  expect_true(all(xrows$.window_end <= 24L))
  expect_true(all(xrows$.window_end - xrows$.window_start + 1L >= 8L))
})

# ── parameter recovery (gated) ─────────────────────────────────────────────

test_that("linear recovers pooled regression coefficients", {
  skip_on_cran()
  skip_if_not_slow()
  dt <- sim_panel_pooled(units = 40, n_time = 30, b0 = 1, b1 = 2, b2 = -1,
                         sd = 1, seed = 1, groupvar = "unit", timevar = "time")
  ctx <- panel_context(unit = "unit", time = "time")
  m <- linearmodel(y ~ x1 + x2, data = dt, ctx = ctx)
  est <- m$coefs
  expect_equal(est$estimate[est$term == "(Intercept)"], 1, tolerance = 0.1)
  expect_equal(est$estimate[est$term == "x1"], 2, tolerance = 0.1)
  expect_equal(est$estimate[est$term == "x2"], -1, tolerance = 0.1)
})

test_that("glm (poisson) recovers the log-mean coefficients", {
  skip_on_cran()
  skip_if_not_slow()
  set.seed(1)
  dt <- data.table::as.data.table(expand.grid(unit = 1:30, time = 1:30))
  dt[, x := stats::rnorm(.N)]
  dt[, y := stats::rpois(.N, exp(0.5 + 0.8 * x))]
  ctx <- panel_context(unit = "unit", time = "time")
  m <- glmmodel(y ~ x, family = stats::poisson(), data = dt, ctx = ctx)
  est <- m$coefs
  expect_equal(est$estimate[est$term == "(Intercept)"], 0.5, tolerance = 0.1)
  expect_equal(est$estimate[est$term == "x"], 0.8, tolerance = 0.1)
})

test_that("heterolm recovers the mean slope and the variance-model sign", {
  skip_on_cran()
  skip_if_not_slow()
  skip_if_not_installed("heterolm")
  set.seed(1)
  dt <- data.table::as.data.table(expand.grid(unit = 1:20, time = 1:30))
  dt[, x := stats::rnorm(.N)]
  dt[, z := stats::rnorm(.N)]
  dt[, y := 1 + 2 * x + stats::rnorm(.N, sd = sqrt(exp(-1 + 1.2 * z)))]
  ctx <- panel_context(unit = "unit", time = "time")
  m <- heterolmmodel(y ~ x, variance = ~ z, data = dt, ctx = ctx)

  co <- data.table::as.data.table(m$fitted$coefficients)
  mean_x <- co[equation == "mean" & term == "x", estimate]
  var_z  <- co[equation == "variance" & term == "z", estimate]
  expect_equal(mean_x, 2, tolerance = 0.15)
  expect_gt(var_z, 0)                       # variance rises with z (correct sign)
  expect_equal(var_z, 1.2, tolerance = 0.4)

  # predicted sigma is increasing in z
  pr <- predict(m$fitted, newdata = as.data.frame(dt), type = "response")
  expect_gt(stats::cor(pr$sigma, dt$z), 0.5)
})
