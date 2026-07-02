# Tests for the gamlss model type
#
# Tier 1 (always-on): formula decomposition, spec construction, graph/closure
#   integration — no gamlss package required.
# Tier 2 (skip_if_no_gamlss): constructor, predict, end-to-end simulate.
# Tier 3 (slow): sigma.formula calibration, non-Gaussian family draws.

library(testthat)

# ── Helper: small balanced panel ─────────────────────────────────────────────

.make_panel_gamlss <- function(units = 6L, n_time = 25L, seed = 1L) {
  set.seed(seed)
  dt <- data.table::CJ(unit = seq_len(units), year = seq_len(n_time))
  dt[, region := factor(ifelse(unit %% 2 == 0, "A", "B"))]
  dt[, x := stats::rnorm(.N)]
  dt[, z := abs(stats::rnorm(.N))]  # positive covariate for sigma
  dt[, y := 0.5 * stats::rnorm(.N)]
  data.table::setkey(dt, unit, year)
  dt[]
}

# ============================================================================
# TIER 1 — No gamlss required
# ============================================================================

# ── .gamlss_decompose ────────────────────────────────────────────────────────

test_that(".gamlss_decompose: grouping var from random() is a graph predictor", {
  f    <- y ~ pb(lag(x)) + random(region)
  sig  <- ~ cs(lag(z))
  dec  <- endogenr:::.gamlss_decompose(f, sig)

  expect_equal(dec$outcome, "y")
  grhs_vars <- all.vars(rlang::f_rhs(dec$graph_formula))
  # pb/cs stripped to underlying vars; grouping factor is an ordinary predictor
  expect_setequal(grhs_vars, c("x", "z", "region"))
})

test_that(".gamlss_decompose: ra() grouping var is a graph predictor", {
  dec  <- endogenr:::.gamlss_decompose(y ~ x + ra(g))
  grhs <- all.vars(rlang::f_rhs(dec$graph_formula))
  expect_true("x" %in% grhs)
  expect_true("g" %in% grhs)
})

test_that(".gamlss_decompose: re(fixed=~x1, random=~1|g) puts x1 and g in graph", {
  dec  <- endogenr:::.gamlss_decompose(y ~ re(fixed = ~x1, random = ~1 | g))
  grhs <- all.vars(rlang::f_rhs(dec$graph_formula))
  expect_true("x1" %in% grhs)
  expect_true("g"  %in% grhs)
})

test_that(".gamlss_decompose: y ~ random(g) yields graph RHS 'g' (grouping is predictor)", {
  dec <- endogenr:::.gamlss_decompose(y ~ random(g))
  expect_equal(dec$outcome, "y")
  expect_equal(all.vars(rlang::f_rhs(dec$graph_formula)), "g")
})

test_that(".gamlss_decompose accumulates predictors across all four formulas", {
  f    <- y ~ pb(lag(x1))
  sig  <- ~ lag(x2)
  nu   <- ~ x3
  tau  <- ~ x4
  dec  <- endogenr:::.gamlss_decompose(f, sig, nu, tau)
  grhs <- all.vars(rlang::f_rhs(dec$graph_formula))
  expect_setequal(grhs, c("x1", "x2", "x3", "x4"))
})

test_that(".gamlss_decompose deduplicates pred_terms across formulas", {
  f    <- y ~ lag(x)
  sig  <- ~ lag(x)    # same predictor in sigma
  dec  <- endogenr:::.gamlss_decompose(f, sig)
  grhs <- all.vars(rlang::f_rhs(dec$graph_formula))
  expect_equal(sum(grhs == "x"), 1L)
})

test_that(".gamlss_decompose: intercept-only mu formula yields a variable-free graph RHS", {
  dec <- endogenr:::.gamlss_decompose(gdppc_grwt ~ 1)
  expect_equal(dec$outcome, "gdppc_grwt")
  expect_equal(all.vars(rlang::f_rhs(dec$graph_formula)), character(0))
})

test_that("validate_system_closure accepts an intercept-only gamlss model", {
  # Regression: graph RHS for `y ~ 1` must not introduce a phantom "1" predictor.
  spec <- list(type = "gamlss", formula = y ~ 1, args = list())
  expect_no_error(validate_system_closure(list(spec), data_columns = "y"))
})

# ── build_model ───────────────────────────────────────────────────────────────

test_that("build_model('gamlss', ...) returns correctly classed spec", {
  f    <- y ~ pb(lag(x)) + random(region)
  spec <- build_model("gamlss",
                      formula       = f,
                      sigma.formula = ~ lag(z),
                      family        = gamlss.dist::NO())

  expect_s3_class(spec, "gamlss_spec")
  expect_s3_class(spec, "endogenr_spec")
  expect_equal(spec$type, "gamlss")
  expect_identical(spec$formula, f)
  expect_true(inherits(spec$args[["sigma.formula"]], "formula"))
})

test_that("build_model('gamlss', ...) stores all four parameter formulas", {
  spec <- build_model("gamlss",
                      formula       = y ~ x,
                      sigma.formula = ~ z1,
                      nu.formula    = ~ z2,
                      tau.formula   = ~ z3)
  expect_true(inherits(spec$args[["sigma.formula"]], "formula"))
  expect_true(inherits(spec$args[["nu.formula"]], "formula"))
  expect_true(inherits(spec$args[["tau.formula"]], "formula"))
})

# ── setup_system / graph-level ────────────────────────────────────────────────

test_that("setup_system succeeds with gamlss spec when grouping var has a producer", {
  dt <- .make_panel_gamlss()

  system <- list(
    build_model("gamlss",
                formula = y ~ pb(lag(x)) + random(region),
                family  = gamlss.dist::NO()),
    build_model("exogen", formula = ~x),
    build_model("exogen", formula = ~region)   # carries the grouping factor forward
  )

  sys <- setup_system(
    models      = system,
    data        = dt,
    train_start = 1L, test_start = 22L, horizon = 3L,
    groupvar    = "unit", timevar = "year", inner_sims = 2L
  )

  expect_s3_class(sys, "endogenr_system_setup")
  expect_true("y" %in% sys$execution_order)
  # region is an ordinary predictor produced by the exogen, so it IS a vertex
  # and must be scheduled before the gamlss outcome that reads it.
  expect_true("region" %in% sys$execution_order)
  expect_lt(match("region", sys$execution_order), match("y", sys$execution_order))
})

test_that("setup_system errors when a gamlss grouping var has no producer", {
  dt <- .make_panel_gamlss()

  system <- list(
    build_model("gamlss",
                formula = y ~ lag(x) + random(region),
                family  = gamlss.dist::NO()),
    build_model("exogen", formula = ~x)
    # no model produces `region` -> must error at setup
  )

  expect_error(
    setup_system(system, dt, train_start = 1L, test_start = 22L, horizon = 3L,
                 groupvar = "unit", timevar = "year", inner_sims = 2L),
    regexp = "region"
  )
})

test_that("setup_system: grouping by the panel key needs no producer", {
  dt <- .make_panel_gamlss()

  system <- list(
    build_model("gamlss",
                formula = y ~ lag(x) + random(unit),  # group by the unit key
                family  = gamlss.dist::NO()),
    build_model("exogen", formula = ~x)
  )

  # `unit` is a panel key — always present in every grid row, so no exogen needed
  expect_no_error(
    setup_system(system, dt, train_start = 1L, test_start = 22L, horizon = 3L,
                 groupvar = "unit", timevar = "year", inner_sims = 2L)
  )
})

test_that("setup_system keeps sigma.formula predictor in keep_cols", {
  dt <- .make_panel_gamlss()

  system <- list(
    build_model("gamlss",
                formula       = y ~ lag(x),
                sigma.formula = ~ lag(z)),
    build_model("exogen", formula = ~x),
    build_model("exogen", formula = ~z)
  )

  # Should not error — z is needed and exogen model supplies it
  expect_no_error(
    setup_system(system, dt, train_start = 1L, test_start = 22L,
                 horizon = 3L, groupvar = "unit", timevar = "year",
                 inner_sims = 2L)
  )
})

test_that("validate_system_closure rejects unmodeled fixed predictor in gamlss", {
  dt <- .make_panel_gamlss()
  dt[, foo := stats::rnorm(.N)]

  system <- list(
    build_model("gamlss",
                formula = y ~ lag(foo) + random(region),
                family  = gamlss.dist::NO()),
    build_model("exogen", formula = ~region)
    # no model produces foo
  )

  expect_error(
    setup_system(system, dt, train_start = 1L, test_start = 22L, horizon = 3L,
                 groupvar = "unit", timevar = "year", inner_sims = 2L),
    regexp = "foo"
  )
})

test_that("validate_system_closure rejects missing sigma.formula predictor", {
  dt <- .make_panel_gamlss()

  system <- list(
    build_model("gamlss",
                formula       = y ~ lag(x),
                sigma.formula = ~ lag(missing_var))
    # no model produces missing_var AND it's not in dt
  )

  expect_error(
    setup_system(system, dt, train_start = 1L, test_start = 22L, horizon = 3L,
                 groupvar = "unit", timevar = "year", inner_sims = 2L),
    regexp = "missing_var"
  )
})

# ============================================================================
# TIER 2 — Requires gamlss
# ============================================================================

test_that("gamlssmodel() constructor works for Normal with pb smoother", {
  skip_if_no_gamlss()
  dt  <- .make_panel_gamlss(units = 8L, n_time = 20L)
  ctx <- panel_context(unit = "unit", time = "year")

  m <- gamlssmodel(
    formula = y ~ pb(lag(x)),
    data    = dt,
    ctx     = ctx
  )

  expect_s3_class(m, "endogenr_gamlss")
  expect_equal(m$outcome, "y")
  expect_true(is.numeric(m$required_history))
  expect_gte(m$required_history, 1L)
  expect_false(is.null(m$fitted))
  expect_false(is.null(m$gamlss_data))
  expect_true(is.data.frame(m$gamlss_data))
  expect_false(is.null(m$graph_formula))
})

test_that("gamlssmodel() with sigma.formula stores correct history depth", {
  skip_if_no_gamlss()
  dt  <- .make_panel_gamlss(units = 6L, n_time = 20L)
  ctx <- panel_context(unit = "unit", time = "year")

  m <- gamlssmodel(
    formula       = y ~ lag(x),
    sigma.formula = ~ lag(z),
    data          = dt,
    ctx           = ctx
  )

  expect_gte(m$required_history, 1L)
  # coefs should be present (mu coefficients)
  expect_false(is.null(m$coefs))
  expect_true(is.data.frame(m$coefs))
})

test_that("predict.gamlss_endogenr returns correct structure: pi mode", {
  skip_if_no_gamlss()

  dt  <- .make_panel_gamlss(units = 6L, n_time = 20L)
  ctx_fit <- panel_context(unit = "unit", time = "year")

  m <- gamlssmodel(
    formula = y ~ pb(lag(x)),
    data    = dt,
    ctx     = ctx_fit
  )

  # Build a sim grid (add sim column)
  sim_ctx <- panel_context(unit = "unit", time = "year", sim = "sim")
  dt2 <- data.table::copy(dt)
  dt2[, sim := 1L]

  pred <- predict(m, data = dt2, t = 18L, ctx = sim_ctx, what = "pi")

  expected_cols <- c("unit", "sim", "year", "y")
  expect_true(all(expected_cols %in% names(pred)))
  expect_true(all(pred$year == 18L))
  expect_equal(nrow(pred), length(unique(dt$unit)))
  expect_false(anyNA(pred$y))
})

test_that("predict.gamlss_endogenr what='expectation' returns mu (no NA)", {
  skip_if_no_gamlss()

  dt  <- .make_panel_gamlss(units = 6L, n_time = 20L)
  ctx_fit <- panel_context(unit = "unit", time = "year")

  m <- gamlssmodel(
    formula = y ~ pb(lag(x)),
    data    = dt,
    ctx     = ctx_fit
  )

  sim_ctx <- panel_context(unit = "unit", time = "year", sim = "sim")
  dt2 <- data.table::copy(dt)
  dt2[, sim := 1L]

  pred_e <- predict(m, data = dt2, t = 18L, ctx = sim_ctx, what = "expectation")
  expect_false(anyNA(pred_e$y))
  expect_equal(nrow(pred_e), length(unique(dt$unit)))
})

test_that("end-to-end: fit_system + simulate_system with gamlss spec", {
  skip_if_no_gamlss()

  dt <- .make_panel_gamlss(units = 6L, n_time = 25L)

  system <- list(
    build_model("gamlss",
                formula = y ~ pb(lag(x)),
                family  = gamlss.dist::NO()),
    build_model("exogen", formula = ~x)
  )

  sys <- setup_system(
    models      = system,
    data        = dt,
    train_start = 1L, test_start = 22L, horizon = 3L,
    groupvar    = "unit", timevar = "year", inner_sims = 4L
  )
  fit <- fit_system(sys, nsim = 4L)
  sim <- simulate_system(fit)

  forecast_rows <- sim[sim$year >= 22L]
  expect_gt(nrow(forecast_rows), 0L)
  expect_false(anyNA(forecast_rows$y))
  expect_true(all(c("unit", "year", ".sim", "y") %in% names(sim)))
})

test_that("gamlssmodel with random() term runs end-to-end", {
  skip_if_no_gamlss()

  dt <- .make_panel_gamlss(units = 6L, n_time = 25L)

  system <- list(
    build_model("gamlss",
                formula = y ~ lag(x) + random(region),
                family  = gamlss.dist::NO()),
    build_model("exogen", formula = ~x),
    build_model("exogen", formula = ~region)   # carry the grouping factor forward
  )

  sys <- setup_system(
    models      = system,
    data        = dt,
    train_start = 1L, test_start = 22L, horizon = 3L,
    groupvar    = "unit", timevar = "year", inner_sims = 2L
  )
  fit <- fit_system(sys, nsim = 2L)
  sim <- simulate_system(fit)

  forecast_rows <- sim[sim$year >= 22L]
  expect_false(anyNA(forecast_rows$y))
})

# ============================================================================
# TIER 3 — Calibration (slow)
# ============================================================================

test_that("sigma.formula heteroscedasticity: spread larger for high-z units (slow)", {
  skip_if_no_gamlss()
  skip_if_not_slow()

  # DGP: variance of y increases with z
  set.seed(42)
  n_u <- 12L; n_t <- 40L
  dt <- data.table::CJ(unit = seq_len(n_u), year = seq_len(n_t))
  dt[, x := stats::rnorm(.N)]
  dt[, z := c(rep(0.2, (.N %/% 2)), rep(2.0, .N - (.N %/% 2)))]
  dt[, y := stats::rnorm(.N, mean = 0.5 * x, sd = exp(0.2 + 0.8 * z))]
  data.table::setkey(dt, unit, year)

  system <- list(
    build_model("gamlss",
                formula       = y ~ lag(x),
                sigma.formula = ~ lag(z),
                family        = gamlss.dist::NO()),
    build_model("exogen", formula = ~x),
    build_model("exogen", formula = ~z)
  )

  sys <- setup_system(system, dt, train_start = 1L, test_start = 36L,
                      horizon = 4L, groupvar = "unit", timevar = "year",
                      inner_sims = 100L)
  fit <- fit_system(sys, nsim = 50L)
  sim <- simulate_system(fit)

  fc <- sim[sim$year >= 36L]
  # Low-z units (z=0.2) vs high-z units (z=2.0): SD should be substantially larger
  low_z_units  <- dt[dt$z < 1, unique(unit)]
  high_z_units <- dt[dt$z > 1, unique(unit)]
  sd_low  <- stats::sd(fc[fc$unit %in% low_z_units,  "y"][[1L]])
  sd_high <- stats::sd(fc[fc$unit %in% high_z_units, "y"][[1L]])
  expect_gt(sd_high, sd_low * 1.5)
})

test_that("BCT (skew) family: fit_system/simulate_system runs and draws are finite (slow)", {
  skip_if_no_gamlss()
  skip_if_not_slow()

  set.seed(123)
  n_u <- 8L; n_t <- 35L
  dt <- data.table::CJ(unit = seq_len(n_u), year = seq_len(n_t))
  dt[, x := abs(stats::rnorm(.N)) + 1]
  # Positive response — BCT needs y > 0
  dt[, y := exp(0.3 * x + stats::rnorm(.N, sd = 0.5))]
  data.table::setkey(dt, unit, year)

  system <- list(
    build_model("gamlss",
                formula = y ~ lag(x),
                family  = gamlss.dist::BCCG()),
    build_model("exogen", formula = ~x)
  )

  sys <- setup_system(system, dt, train_start = 1L, test_start = 32L,
                      horizon = 3L, groupvar = "unit", timevar = "year",
                      inner_sims = 20L)
  fit <- fit_system(sys, nsim = 10L)
  sim <- simulate_system(fit)

  fc <- sim[sim$year >= 32L]
  expect_false(anyNA(fc$y))
  expect_true(all(is.finite(fc$y)))
  # BCCG is positive-support: draws should be > 0
  expect_true(all(fc$y > 0))
})
