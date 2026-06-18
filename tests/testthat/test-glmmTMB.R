# Tests for the glmmTMB model type
#
# Tier 1 (always-on): formula decomposition, spec construction, graph/closure
#   integration — no glmmTMB package required.
# Tier 2 (skip_if_no_glmmTMB): constructor, predict, end-to-end simulate.
# Tier 3 (slow): dispformula calibration, count-family draws.

library(testthat)

# ── Helper: small balanced panel with a grouping column ──────────────────────

.make_panel_glmmtmb <- function(units = 6L, n_time = 25L, seed = 1L) {
  set.seed(seed)
  dt <- data.table::CJ(unit = seq_len(units), year = seq_len(n_time))
  dt[, region := factor(ifelse(unit %% 2 == 0, "A", "B"))]
  dt[, x := stats::rnorm(.N)]
  dt[, y := stats::rnorm(.N)]
  data.table::setkey(dt, unit, year)
  dt[]
}

# ============================================================================
# TIER 1 — No glmmTMB required
# ============================================================================

test_that(".glmmTMB_decompose: grouping factor is a graph predictor", {
  f <- gdppc_grwt ~ lag(dem) + lag(log(gdppc)) + lag(yj01best) +
         (1 + lag(dem) | region)
  disp <- ~ lag(dem) + lag(yj01best)
  dec <- endogenr:::.glmmTMB_decompose(f, disp)

  expect_equal(dec$outcome, "gdppc_grwt")

  # Effect vars AND the grouping factor are all ordinary predictors (edges)
  grhs_vars <- all.vars(rlang::f_rhs(dec$graph_formula))
  expect_setequal(grhs_vars, c("dem", "gdppc", "yj01best", "region"))
})

test_that(".glmmTMB_decompose: cov-struct coords and grouping are graph predictors", {
  f <- y ~ x + ar1(times + 0 | g)
  dec <- endogenr:::.glmmTMB_decompose(f)

  # Coordinates and grouping are ordinary predictors (must be produced/carried)
  grhs_vars <- all.vars(rlang::f_rhs(dec$graph_formula))
  expect_setequal(grhs_vars, c("x", "times", "g"))
})

test_that(".glmmTMB_decompose: intercept-only RE (1 | region) -> region predictor", {
  f <- y ~ x + (1 | region)
  dec <- endogenr:::.glmmTMB_decompose(f)

  grhs <- all.vars(rlang::f_rhs(dec$graph_formula))
  expect_true("region" %in% grhs)
  expect_true("x" %in% grhs)
})

test_that(".glmmTMB_decompose: y ~ (1 | g) yields graph RHS 'g' (grouping is predictor)", {
  f <- y ~ (1 | g)
  dec <- endogenr:::.glmmTMB_decompose(f)

  expect_equal(dec$outcome, "y")
  expect_equal(all.vars(rlang::f_rhs(dec$graph_formula)), "g")
})

test_that("build_model('glmmTMB', ...) returns correctly classed spec", {
  f <- y ~ lag(x) + (1 | region)
  spec <- build_model("glmmTMB", formula = f,
                      dispformula = ~ lag(x),
                      family = stats::gaussian())

  expect_s3_class(spec, "glmmTMB_spec")
  expect_s3_class(spec, "endogenr_spec")
  expect_equal(spec$type, "glmmTMB")
  expect_identical(spec$formula, f)
  expect_true(inherits(spec$args$dispformula, "formula"))
  expect_true(is.function(spec$args$family) || inherits(spec$args$family, "family"))
})

test_that("setup_system succeeds with glmmTMB spec when grouping var has a producer", {
  dt <- .make_panel_glmmtmb()

  system <- list(
    build_model("glmmTMB",
                formula  = y ~ lag(x) + (1 | region),
                family   = stats::gaussian()),
    build_model("exogen", formula = ~x),
    build_model("exogen", formula = ~region)   # carries the grouping factor forward
  )

  sys <- setup_system(
    models     = system,
    data       = dt,
    train_start = 1L, test_start = 22L, horizon = 3L,
    groupvar   = "unit", timevar = "year", inner_sims = 2L
  )

  expect_s3_class(sys, "endogenr_system_setup")
  expect_true("y" %in% sys$execution_order)
  # region is an ordinary predictor produced by the exogen, so it IS a vertex
  # and is scheduled before the glmmTMB outcome that reads it.
  expect_true("region" %in% sys$execution_order)
  expect_lt(match("region", sys$execution_order), match("y", sys$execution_order))
})

test_that("setup_system errors when a glmmTMB grouping var has no producer", {
  dt <- .make_panel_glmmtmb()

  system <- list(
    build_model("glmmTMB",
                formula = y ~ lag(x) + (1 | region),
                family  = stats::gaussian()),
    build_model("exogen", formula = ~x)
    # no model produces `region` -> must error at setup
  )

  expect_error(
    setup_system(system, dt, train_start = 1L, test_start = 22L, horizon = 3L,
                 groupvar = "unit", timevar = "year", inner_sims = 2L),
    regexp = "region"
  )
})

test_that("validate_system_closure rejects truly unmodeled fixed predictor in glmmTMB", {
  dt <- .make_panel_glmmtmb()
  dt[, foo := stats::rnorm(.N)]

  system <- list(
    build_model("glmmTMB",
                formula = y ~ lag(foo) + (1 | region),
                family  = stats::gaussian()),
    build_model("exogen", formula = ~region)
    # no model produces foo
  )

  expect_error(
    setup_system(system, dt, train_start = 1L, test_start = 22L, horizon = 3L,
                 groupvar = "unit", timevar = "year", inner_sims = 2L),
    regexp = "foo"
  )
})

test_that(".is_trivial_rhs identifies ~1 and ~0 formulas", {
  expect_true(endogenr:::.is_trivial_rhs(~1))
  expect_true(endogenr:::.is_trivial_rhs(~0))
  expect_true(endogenr:::.is_trivial_rhs(NULL))
  expect_false(endogenr:::.is_trivial_rhs(~ x))
  expect_false(endogenr:::.is_trivial_rhs(~ lag(x)))
})

test_that(".glmmTMB_response_draw returns correct length and sign", {
  mu <- c(1, 2, 3)
  disp <- c(0.5, 0.5, 0.5)
  set.seed(1)

  d_gauss <- endogenr:::.glmmTMB_response_draw(mu, "gaussian", disp)
  expect_length(d_gauss, 3L)

  d_pois <- endogenr:::.glmmTMB_response_draw(mu, "poisson", disp)
  expect_true(all(d_pois >= 0))
  expect_true(all(d_pois == as.integer(d_pois)))  # non-negative integers

  d_nb2 <- endogenr:::.glmmTMB_response_draw(mu, "nbinom2", disp)
  expect_true(all(d_nb2 >= 0))

  d_nb1 <- endogenr:::.glmmTMB_response_draw(mu, "nbinom1", disp)
  expect_true(all(d_nb1 >= 0))

  d_gamma <- endogenr:::.glmmTMB_response_draw(mu, "Gamma", disp)
  expect_true(all(d_gamma > 0))

  d_beta <- endogenr:::.glmmTMB_response_draw(c(0.3, 0.5, 0.7), "beta", disp)
  expect_true(all(d_beta >= 0 & d_beta <= 1))

  # Unsupported family warns and returns mu
  expect_warning(
    d_unk <- endogenr:::.glmmTMB_response_draw(mu, "genpois", disp),
    regexp = "no response-scale predictive draw"
  )
  expect_equal(d_unk, mu)
})

test_that(".glmmTMB_response_draw: new continuous families match target moments", {
  set.seed(1); n <- 1e5
  dt_ <- endogenr:::.glmmTMB_response_draw(rep(5, n), "t", rep(2, n),
                                           c("Student-t df" = 10))
  expect_equal(mean(dt_), 5, tolerance = 0.05)
  expect_equal(stats::sd(dt_), 2 * sqrt(10 / 8), tolerance = 0.1)   # disp*sqrt(df/(df-2))

  dl <- endogenr:::.glmmTMB_response_draw(rep(10, n), "lognormal", rep(4, n))
  expect_true(all(dl > 0))
  expect_equal(mean(dl), 10, tolerance = 0.2)
  expect_equal(stats::sd(dl), 4, tolerance = 0.3)

  ds <- endogenr:::.glmmTMB_response_draw(rep(0, n), "skewnormal", rep(1, n),
                                          c("Skewnormal shape" = 6))
  expect_equal(mean(ds), 0, tolerance = 0.05)
  expect_gt(mean(((ds - mean(ds)) / stats::sd(ds))^3), 0.2)         # right-skewed
})

test_that(".glmmTMB_response_draw: truncated counts are positive integers", {
  set.seed(2); n <- 1e5
  for (fam_disp in list(c("truncated_poisson", NA), c("truncated_nbinom2", 2),
                        c("truncated_nbinom1", 0.5))) {
    fam  <- fam_disp[[1]]
    disp <- if (is.na(fam_disp[[2]])) rep(1, n) else rep(as.numeric(fam_disp[[2]]), n)
    d <- endogenr:::.glmmTMB_response_draw(rep(0.7, n), fam, disp)
    expect_true(all(d >= 1), info = fam)
    expect_true(all(d == as.integer(d)), info = fam)
  }
})

test_that(".glmmTMB_response_draw: tweedie draw (requires tweedie pkg)", {
  skip_if_not_installed("tweedie")
  set.seed(3); n <- 1e5
  out <- endogenr:::.glmmTMB_response_draw(rep(5, n), "tweedie", rep(2, n),
                                           c("Tweedie power" = 1.5))
  expect_true(all(out >= 0))
  expect_equal(mean(out), 5, tolerance = 0.2)
})

test_that(".glmmTMB_linkinv returns a function", {
  skip_if_no_glmmTMB()
  # Build a tiny model to test linkinv extraction
  dt <- .make_panel_glmmtmb(units = 4L, n_time = 10L)
  fit <- glmmTMB::glmmTMB(y ~ x, data = as.data.frame(dt), family = stats::gaussian())
  linkinv <- endogenr:::.glmmTMB_linkinv(fit)
  expect_true(is.function(linkinv))
  # Identity link for gaussian: linkinv(1.5) == 1.5
  expect_equal(linkinv(1.5), 1.5)
})

# ============================================================================
# TIER 2 — Requires glmmTMB
# ============================================================================

test_that("glmmTMBmodel() constructor works for gaussian with random intercept", {
  skip_if_no_glmmTMB()
  dt <- .make_panel_glmmtmb(units = 8L, n_time = 20L)
  ctx <- panel_context(unit = "unit", time = "year")

  m <- glmmTMBmodel(
    formula     = y ~ lag(x) + (1 | region),
    dispformula = ~1,
    ziformula   = ~0,
    family      = stats::gaussian(),
    data        = dt,
    ctx         = ctx
  )

  expect_s3_class(m, "glmmTMB_endogenr")
  expect_equal(m$outcome, "y")
  expect_true(is.numeric(m$required_history))
  expect_gte(m$required_history, 1L)
  expect_true(!is.null(m$fitted))
  expect_true(!is.null(m$coefs))
  expect_true(is.data.frame(m$coefs))
  expect_true(nrow(m$coefs) >= 1L)
  expect_true(all(c("term", "estimate", "std.error") %in% names(m$coefs)))
  expect_true(!is.null(m$graph_formula))
})

test_that("predict.glmmTMB_endogenr returns correct structure at time t", {
  skip_if_no_glmmTMB()

  dt <- .make_panel_glmmtmb(units = 6L, n_time = 20L)
  ctx_fit  <- panel_context(unit = "unit", time = "year")

  m <- glmmTMBmodel(
    formula = y ~ lag(x) + (1 | region),
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
  expect_true(!anyNA(pred$y))

  # what = "expectation"
  pred_e <- predict(m, data = dt2, t = 18L, ctx = sim_ctx, what = "expectation")
  expect_true(!anyNA(pred_e$y))
})

test_that("end-to-end: fit_system + simulate_system with glmmTMB spec", {
  skip_if_no_glmmTMB()

  dt <- .make_panel_glmmtmb(units = 6L, n_time = 25L)

  system <- list(
    build_model("glmmTMB",
                formula = y ~ lag(x) + (1 | region),
                family  = stats::gaussian()),
    build_model("exogen", formula = ~x),
    build_model("exogen", formula = ~region)   # carry the grouping factor forward
  )

  sys <- setup_system(
    models      = system,
    data        = dt,
    train_start = 1L, test_start = 22L, horizon = 3L,
    groupvar    = "unit", timevar = "year", inner_sims = 3L
  )
  fit <- fit_system(sys, nsim = 4L)
  sim <- simulate_system(fit)

  # Forecast window rows should be fully populated
  forecast_rows <- sim[sim$year >= 22L]
  expect_gt(nrow(forecast_rows), 0L)
  expect_false(anyNA(forecast_rows$y))

  # Correct column set
  expect_true(all(c("unit", "year", ".sim", "y") %in% names(sim)))
})

test_that("trivial dispformula caches constant dispersion (disp_const)", {
  skip_if_no_glmmTMB()

  dt  <- .make_panel_glmmtmb(units = 6L, n_time = 20L)
  ctx <- panel_context(unit = "unit", time = "year")

  # Default dispformula (~1): disp_const is the scalar sigma(fitted). (The per-row
  # equivalence disp == sigma(fitted) is pinned by the ar1 disp test below, which
  # materialises the predictor columns the fit expects.)
  m <- glmmTMBmodel(y ~ lag(x) + (1 | region), data = dt, ctx = ctx)
  expect_true(is.numeric(m$disp_const) && length(m$disp_const) == 1L)
  expect_equal(m$disp_const, as.numeric(stats::sigma(m$fitted)))

  # Cached family + link-inverse match the live ones.
  expect_equal(m$fam_name, family(m$fitted)$family)

  # Non-trivial dispformula: no constant cache -> predict method takes per-row path.
  m_h <- glmmTMBmodel(y ~ lag(x), dispformula = ~ x, data = dt, ctx = ctx)
  expect_null(m_h$disp_const)
})

test_that("dispformula: per-row dispersion reaches predict draw (slow)", {
  skip_if_no_glmmTMB()
  skip_if_not_slow()

  # DGP: variance increases in x. x is a UNIT-LEVEL covariate (constant within
  # unit) so the heteroscedastic signal persists across the whole forecast
  # horizon — at every step lag(x) stays high for high-x units. (An IID-per-row
  # x would make lag(x) random noise after the first forecast step, diluting the
  # signal the test is meant to detect.)
  set.seed(7)
  n_u <- 10L; n_t <- 30L
  dt <- data.table::CJ(unit = seq_len(n_u), year = seq_len(n_t))
  unit_x <- seq(-1.5, 1.5, length.out = n_u)
  dt[, x := unit_x[unit]]
  dt[, y := stats::rnorm(.N, mean = 0.5 * x, sd = exp(0.5 + 0.8 * x))]
  data.table::setkey(dt, unit, year)

  system <- list(
    build_model("glmmTMB",
                formula     = y ~ lag(x),
                dispformula = ~ lag(x),
                family      = stats::gaussian()),
    build_model("exogen", formula = ~x)
  )

  sys <- setup_system(system, dt, train_start = 1L, test_start = 27L,
                      horizon = 3L, groupvar = "unit", timevar = "year",
                      inner_sims = 50L)
  fit <- fit_system(sys, nsim = 20L)
  sim <- simulate_system(fit)

  # Units with high x should have wider intervals than units with low x
  # Compare mean interval width for top vs bottom x units
  origin <- dt[year == 26L, .(unit, x)]
  high_x <- origin[x > stats::quantile(x, 0.75)]$unit
  low_x  <- origin[x < stats::quantile(x, 0.25)]$unit

  high_sim <- sim[unit %in% high_x & year >= 27L]
  low_sim  <- sim[unit %in% low_x  & year >= 27L]

  sd_high <- stats::sd(high_sim$y)
  sd_low  <- stats::sd(low_sim$y)

  # High-x draws should be more spread
  expect_gt(sd_high, sd_low * 1.2)
})

test_that("count family: nbinom2 draws are non-negative integers (slow)", {
  skip_if_no_glmmTMB()
  skip_if_not_slow()

  set.seed(11)
  n_u <- 8L; n_t <- 35L
  dt <- data.table::CJ(unit = seq_len(n_u), year = seq_len(n_t))
  dt[, x := stats::rnorm(.N)]
  # NB2 DGP: mean = exp(0.5 + 0.7 * x), size = 3
  dt[, y := stats::rnbinom(.N, mu = exp(0.5 + 0.7 * dt$x), size = 3)]
  data.table::setkey(dt, unit, year)

  system <- list(
    build_model("glmmTMB",
                formula = y ~ lag(x),
                family  = glmmTMB::nbinom2()),
    build_model("exogen", formula = ~x)
  )

  sys <- setup_system(system, dt, train_start = 1L, test_start = 30L,
                      horizon = 5L, groupvar = "unit", timevar = "year",
                      inner_sims = 100L)
  fit <- fit_system(sys, nsim = 30L)
  sim <- simulate_system(fit)

  fc <- sim[year >= 30L]
  expect_true(!anyNA(fc$y))
  expect_true(all(fc$y >= 0))
  expect_true(all(fc$y == as.integer(fc$y)))

  # Marginal mean recovery: should be within 50% of truth (loose; small sample)
  grand_mu <- mean(fc$y)
  true_mu  <- mean(exp(0.5 + 0.7 * dt[year < 30L]$x))
  expect_lt(abs(grand_mu / true_mu - 1), 0.5)
})

# ============================================================================
# TIER 2/3 — Temporal covariance structures (ar1/ou/…) forecasting
# ============================================================================
#
# A temporal cov-struct RE (e.g. ar1(years + 0 | gwcode)) must be forecast with
# the WHOLE forecast-so-far block in one prediction so glmmTMB applies the
# correct phi^k correlation decay. predict.glmmTMB_endogenr does this for any
# model whose $has_covstruct flag is TRUE; these tests pin that behaviour and
# the related dispersion-prediction (interval-scale) fix.

# Balanced panel with a genuine per-group AR(1) latent over the time coordinate,
# so the estimated correlation phi and residual sigma are solidly non-degenerate
# and the phi^k identity is discriminating.
.make_panel_ar1 <- function(n_gw = 8L, years = 2000:2018,
                            phi_true = 0.7, seed = 99L) {
  set.seed(seed)
  dt <- data.table::CJ(gwcode = seq_len(n_gw), year = years)
  re_path <- function(n) {
    e <- stats::rnorm(n, sd = 0.5)
    r <- numeric(n)
    r[1L] <- e[1L]
    for (i in seq.int(2L, n)) r[i] <- phi_true * r[i - 1L] + e[i]
    r
  }
  dt[, lat := re_path(.N), by = gwcode]
  dt[, y := 1 + lat + stats::rnorm(.N, sd = 0.6)]
  dt[, lat := NULL]
  dt[, years := factor(year)]
  dt[, gwcode := factor(gwcode)]
  data.table::setkey(dt, gwcode, year)
  dt[]
}

test_that("glmmTMBmodel: has_covstruct + last_train_time set at fit", {
  skip_if_no_glmmTMB()

  dt  <- .make_panel_ar1()
  ctx <- panel_context(unit = "gwcode", time = "year")
  m_ar <- glmmTMBmodel(y ~ ar1(years + 0 | gwcode),
                       data = dt[year <= 2015L], ctx = ctx)
  expect_true(isTRUE(m_ar$has_covstruct))
  expect_equal(m_ar$last_train_time, 2015)

  # A plain random intercept is NOT a covariance structure.
  dt2  <- .make_panel_glmmtmb(units = 6L, n_time = 20L)
  ctx2 <- panel_context(unit = "unit", time = "year")
  m_re <- glmmTMBmodel(y ~ lag(x) + (1 | region), data = dt2, ctx = ctx2)
  expect_false(isTRUE(m_re$has_covstruct))
  expect_equal(m_re$last_train_time, max(dt2$year))
})

test_that("glmmTMB ar1: step-by-step predict yields phi^k decay (block path)", {
  skip_if_no_glmmTMB()

  dt  <- .make_panel_ar1()
  ctx <- panel_context(unit = "gwcode", time = "year")
  m   <- glmmTMBmodel(y ~ ar1(years + 0 | gwcode),
                      data = dt[year <= 2015L], ctx = ctx)

  phi <- attr(glmmTMB::VarCorr(m$fitted)$cond$gwcode, "correlation")[1, 2]
  re  <- glmmTMB::ranef(m$fitted)$cond$gwcode
  b0  <- glmmTMB::fixef(m$fitted)$cond[[1]]
  expect_match(colnames(re)[ncol(re)], "2015")   # last col = last training year
  re_last <- as.numeric(re[, ncol(re)])
  expect_gt(abs(phi), 0.1)                        # AR1 signal present

  sim_ctx <- panel_context(unit = "gwcode", time = "year", sim = "sim")
  g <- data.table::copy(dt)
  g[, sim := 1L]

  # Conditional mean at forecast step k must be intercept + phi^k * RE(last train).
  for (k in 1:3) {
    p <- predict(m, data = g, t = 2015L + k, ctx = sim_ctx, what = "expectation")
    data.table::setorder(p, gwcode)
    expect_equal(p$y - b0, phi^k * re_last, tolerance = 1e-4)
  }

  # Discrimination vs the old single-step bug (which returned phi^1 every step):
  # the step-2 random effect must differ from the step-1 random effect.
  p1 <- predict(m, data = g, t = 2016L, ctx = sim_ctx, what = "expectation")
  p2 <- predict(m, data = g, t = 2017L, ctx = sim_ctx, what = "expectation")
  data.table::setorder(p1, gwcode)
  data.table::setorder(p2, gwcode)
  expect_false(isTRUE(all.equal(p2$y, p1$y)))
})

test_that("glmmTMB ar1: disp predict needs allow.new.levels (interval-scale fix)", {
  skip_if_no_glmmTMB()

  dt  <- .make_panel_ar1()
  ctx <- panel_context(unit = "gwcode", time = "year")
  m   <- glmmTMBmodel(y ~ ar1(years + 0 | gwcode),
                      data = dt[year <= 2015L], ctx = ctx)
  blk <- as.data.frame(dt[year >= 2016L & year <= 2018L])

  # Without the flag glmmTMB cannot place new cov-struct levels -> it errors,
  # which the predict method's tryCatch swallowed, falling back to dispersion = 1
  # and inflating gaussian predictive intervals. The fix passes the flag so the
  # true model sigma reaches the response draw.
  expect_error(stats::predict(m$fitted, newdata = blk, type = "disp"))
  d <- stats::predict(m$fitted, newdata = blk, type = "disp",
                      allow.new.levels = TRUE)
  expect_equal(as.numeric(d), rep(sigma(m$fitted), nrow(blk)), tolerance = 1e-6)
})

test_that("end-to-end: ar1 cov-struct simulate runs, no NA, no new-level warnings", {
  skip_if_no_glmmTMB()

  dt <- .make_panel_ar1()
  dt[, x := stats::rnorm(.N)]   # real predictor -> lag(x) materialised in the block

  system <- list(
    build_model("glmmTMB", formula = y ~ lag(x) + ar1(years + 0 | gwcode),
                family = stats::gaussian()),
    build_model("exogen", formula = ~x),
    build_model("exogen", formula = ~years)
  )

  warns <- character(0)
  sim <- withCallingHandlers(
    {
      sys <- setup_system(system, dt, train_start = 2000L, test_start = 2016L,
                          horizon = 3L, groupvar = "gwcode", timevar = "year",
                          inner_sims = 3L)
      f <- fit_system(sys, nsim = 4L)
      simulate_system(f)
    },
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  fc <- sim[year >= 2016L]
  expect_gt(nrow(fc), 0L)
  expect_false(anyNA(fc$y))
  expect_true(all(c("gwcode", "year", ".sim", "y") %in% names(sim)))
  expect_false(any(grepl("new random effect levels", warns)))
})

test_that("glmmTMB ar1: RE contribution decays over a longer horizon (slow)", {
  skip_if_no_glmmTMB()
  skip_if_not_slow()

  dt  <- .make_panel_ar1(years = 2000:2020)
  ctx <- panel_context(unit = "gwcode", time = "year")
  m   <- glmmTMBmodel(y ~ ar1(years + 0 | gwcode),
                      data = dt[year <= 2015L], ctx = ctx)

  phi <- attr(glmmTMB::VarCorr(m$fitted)$cond$gwcode, "correlation")[1, 2]
  re  <- glmmTMB::ranef(m$fitted)$cond$gwcode
  b0  <- glmmTMB::fixef(m$fitted)$cond[[1]]
  re_last <- as.numeric(re[, ncol(re)])

  sim_ctx <- panel_context(unit = "gwcode", time = "year", sim = "sim")
  g <- data.table::copy(dt)
  g[, sim := 1L]
  p1 <- predict(m, data = g, t = 2016L, ctx = sim_ctx, what = "expectation")
  p5 <- predict(m, data = g, t = 2020L, ctx = sim_ctx, what = "expectation")
  data.table::setorder(p1, gwcode)
  data.table::setorder(p5, gwcode)

  expect_equal(p5$y - b0, phi^5 * re_last, tolerance = 1e-4)  # identity holds at k=5
  expect_lt(mean(abs(p5$y - b0)), mean(abs(p1$y - b0)))       # magnitude decays
})

# ============================================================================
# TIER 2 — Additional family response draws (t/lognormal/truncated_poisson)
# ============================================================================

test_that("glmmTMBmodel: t_family stores df and predicts finite draws", {
  skip_if_no_glmmTMB()
  dt <- .make_panel_glmmtmb(units = 8L, n_time = 25L)
  set.seed(5); dt[, y := 2 + 0.5 * x + stats::rt(.N, df = 6)]
  ctx <- panel_context(unit = "unit", time = "year")
  m <- glmmTMBmodel(y ~ lag(x), data = dt, ctx = ctx, family = glmmTMB::t_family())
  expect_true("Student-t df" %in% names(m$family_params))
  sim_ctx <- panel_context(unit = "unit", time = "year", sim = "sim")
  g <- data.table::copy(dt); g[, sim := 1L]
  p <- predict(m, data = g, t = 18L, ctx = sim_ctx, what = "pi")
  expect_false(anyNA(p$y))
})

test_that("glmmTMBmodel: lognormal predicts strictly positive draws", {
  skip_if_no_glmmTMB()
  dt <- .make_panel_glmmtmb(units = 6L, n_time = 25L)
  set.seed(6); dt[, y := stats::rlnorm(.N, meanlog = 0.2 + 0.3 * x, sdlog = 0.4)]
  ctx <- panel_context(unit = "unit", time = "year")
  m <- glmmTMBmodel(y ~ lag(x), data = dt, ctx = ctx, family = glmmTMB::lognormal())
  sim_ctx <- panel_context(unit = "unit", time = "year", sim = "sim")
  g <- data.table::copy(dt); g[, sim := 1L]
  p <- predict(m, data = g, t = 18L, ctx = sim_ctx, what = "pi")
  expect_true(all(p$y > 0)); expect_false(anyNA(p$y))
})

test_that("glmmTMBmodel: truncated_poisson predicts positive integer draws", {
  skip_if_no_glmmTMB()
  dt <- .make_panel_glmmtmb(units = 6L, n_time = 25L)
  set.seed(7); lam <- exp(0.3 + 0.3 * dt$x)
  dt[, y := stats::qpois(stats::runif(.N, stats::dpois(0L, lam), 1), lam)]
  ctx <- panel_context(unit = "unit", time = "year")
  m <- glmmTMBmodel(y ~ lag(x), data = dt, ctx = ctx,
                    family = glmmTMB::truncated_poisson())
  sim_ctx <- panel_context(unit = "unit", time = "year", sim = "sim")
  g <- data.table::copy(dt); g[, sim := 1L]
  p <- predict(m, data = g, t = 18L, ctx = sim_ctx, what = "pi")
  expect_true(all(p$y >= 1)); expect_true(all(p$y == as.integer(p$y)))
})

test_that("end-to-end: t_family simulate runs and forecasts finite values", {
  skip_if_no_glmmTMB()
  dt <- .make_panel_glmmtmb(units = 6L, n_time = 25L)
  set.seed(8); dt[, y := 1 + 0.4 * x + stats::rt(.N, df = 8)]
  system <- list(
    build_model("glmmTMB", formula = y ~ lag(x), family = glmmTMB::t_family()),
    build_model("exogen", formula = ~x)
  )
  sys <- setup_system(system, dt, train_start = 1L, test_start = 22L, horizon = 3L,
                      groupvar = "unit", timevar = "year", inner_sims = 3L)
  sim <- simulate_system(fit_system(sys, nsim = 4L))
  fc <- sim[year >= 22L]
  expect_gt(nrow(fc), 0L); expect_false(anyNA(fc$y))
})
