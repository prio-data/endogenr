# Tests for R/longhorizon.R ------------------------------------------------

# ── Helper: balanced synthetic panel with AR(1) signal ────────────────────

lh_test_data <- function() {
  set.seed(42)
  units <- 1:6
  years <- 1990:2012
  rows <- lapply(units, function(g) {
    n <- length(years)
    y <- numeric(n)
    y[1] <- stats::rnorm(1)
    for (k in 2:n) y[k] <- 0.6 * y[k - 1] + stats::rnorm(1, sd = 0.3)
    data.frame(gwcode = g, year = years, y = y, x = stats::rnorm(n))
  })
  df <- do.call(rbind, rows)
  df$gdppc <- 1000 * exp(0.3 * df$y)  # strictly positive, for asinh tests
  data.table::as.data.table(df)
}

# ── lead_horizon ──────────────────────────────────────────────────────────

test_that("lead_horizon performs a positional lead with NA padding", {
  expect_equal(lead_horizon(1:5, 2), c(3L, 4L, 5L, NA, NA))
  expect_equal(lead_horizon(c(1.5, 2.5, 3.5), 1), c(2.5, 3.5, NA))
  expect_identical(lead_horizon(1:3, 0), 1:3)
  expect_true(all(is.na(lead_horizon(1:3, 5))))
})

# ── create_lh_data ──────────────────────────────────────────────────────────

test_that("create_lh_data builds .target as lead(outcome, h) and filters leakage", {
  dt <- lh_test_data()
  aligned <- create_lh_data(dt, "y", h = 2, "gwcode", "year", test_start = 2010)

  expect_true(".target" %in% names(aligned))
  # leakage filter: t + h < test_start  <=>  t < test_start - h = 2008
  expect_true(all(aligned$year < 2008))

  # .target at (unit, year) equals the raw outcome h steps ahead
  g1 <- dt[gwcode == 1]
  expect_equal(aligned[gwcode == 1 & year == 2000, .target],
               g1[year == 2002, y])
})

# ── formula convention (lead_horizon marker) ───────────────────────────────

test_that("setup_long_horizon requires the LHS to be wrapped in lead_horizon()", {
  dt <- lh_test_data()
  expect_error(
    setup_long_horizon(dt, list(bad = y ~ x), 1, "gwcode", "year", 2010),
    "lead_horizon"
  )
})

test_that("RHS is evaluated at the origin (no implicit lag)", {
  dt <- lh_test_data()
  f <- list(lh = lead_horizon(y) ~ y + x)
  setup <- setup_long_horizon(dt, f, 1, "gwcode", "year", test_start = 2010)
  # naive formula predictors are the origin covariates, not lagged copies
  nf <- setup$models$lh[["1"]]$naive_formula
  expect_setequal(all.vars(nf), c(".target", "y", "x"))
})

# ── setup → forecast → get_lh_accuracy (native) ─────────────────────────────

test_that("forecast_long_horizon maps horizons to year_target and defaults test_start", {
  dt <- lh_test_data()
  f <- list(lh = lead_horizon(y) ~ y + x)
  setup <- setup_long_horizon(dt, f, 1:3, "gwcode", "year", test_start = 2010)
  fc <- forecast_long_horizon(setup, dt, nsim = 30, inner_sims = 4)

  expect_equal(sort(unique(fc$horizon)), 1:3)
  expect_equal(unique(fc$test_start), 2010)                  # defaults to setup test_start
  expect_equal(sort(unique(fc$year_target)), c(2010, 2011, 2012))
  expect_equal(length(fc$.draws[[1]]), 30 * 4)
  expect_false(is.null(attr(fc, "lh_formulas")))             # stamped for scoring
})

test_that("get_lh_accuracy scores per (variant, unit, horizon)", {
  dt <- lh_test_data()
  f <- list(lh = lead_horizon(y) ~ y + x)
  setup <- setup_long_horizon(dt, f, 1:3, "gwcode", "year", test_start = 2010)
  fc <- forecast_long_horizon(setup, dt, nsim = 30, inner_sims = 4)
  acc <- get_lh_accuracy(fc, dt, setup)

  expect_true(all(c("variant", "gwcode", "horizon", "crps", "mae", "winkler")
                  %in% names(acc)))
  expect_true(all(is.finite(acc$crps)))
  expect_true(all(acc$crps >= 0))
  expect_true(all(acc$winkler >= 0))
})

test_that("forecast_long_horizon warns when test_start precedes the setup test_start", {
  dt <- lh_test_data()
  f <- list(lh = lead_horizon(y) ~ y + x)
  setup <- setup_long_horizon(dt, f, 1, "gwcode", "year", test_start = 2010)
  expect_warning(
    forecast_long_horizon(setup, dt, test_start = 2005, nsim = 5, inner_sims = 2),
    "leak"
  )
})

# ── regression: NA-safe scoring (was a crash after the refactor) ────────────

test_that("get_lh_accuracy tolerates missing truth and NA draws without erroring", {
  dt <- lh_test_data()
  f <- list(lh = lead_horizon(y) ~ y + x)
  setup <- setup_long_horizon(dt, f, 1:2, "gwcode", "year", test_start = 2010)
  fc <- forecast_long_horizon(setup, dt, nsim = 20, inner_sims = 3)

  # Drop the truth row for one forecast cell (gwcode 1, year_target 2010 = h1)
  truth_missing <- data.table::copy(dt)[!(gwcode == 1 & year == 2010)]
  expect_error(acc <- get_lh_accuracy(fc, truth_missing, setup),
               NA)
  expect_true(is.na(acc[gwcode == 1 & horizon == 1, crps]))

  # NA values inside the draws must not crash crps_sample
  fc2 <- data.table::copy(fc)
  fc2$.draws[[1]][1:3] <- NA_real_
  expect_error(get_lh_accuracy(fc2, dt, setup), NA)
})

# ── scale: native vs model ──────────────────────────────────────────────────

test_that("scale='model' and scale='native' agree for an identity LHS", {
  dt <- lh_test_data()
  f <- list(lh = lead_horizon(y) ~ y + x)
  setup <- setup_long_horizon(dt, f, 1:2, "gwcode", "year", test_start = 2010)
  fc <- forecast_long_horizon(setup, dt, nsim = 20, inner_sims = 3)

  an <- get_lh_accuracy(fc, dt, setup, scale = "native")
  am <- get_lh_accuracy(fc, dt, setup, scale = "model")
  data.table::setkey(an, gwcode, horizon)
  data.table::setkey(am, gwcode, horizon)
  expect_equal(an$crps, am$crps)
})

test_that("scale handles transformed LHS (asinh)", {
  dt <- lh_test_data()
  ft <- list(lt = lead_horizon(asinh(gdppc)) ~ asinh(gdppc) + x)
  st <- setup_long_horizon(dt, ft, 1:2, "gwcode", "year", test_start = 2010)
  fct <- forecast_long_horizon(st, dt, nsim = 20, inner_sims = 3)

  am <- get_lh_accuracy(fct, dt, st, scale = "model")
  an <- get_lh_accuracy(fct, dt, st, scale = "native", inverse = sinh)
  expect_true(all(is.finite(am$crps)))
  expect_true(all(is.finite(an$crps)))
  # native is on the (large) gdppc scale; model is on the compressed asinh scale
  expect_gt(stats::median(an$crps), stats::median(am$crps))

  # native scoring of a transformed LHS needs an inverse
  expect_error(
    get_lh_accuracy(fct, dt, st, scale = "native"),
    "inverse"
  )
})

# ── compare_approaches shares (unit, horizon) ───────────────────────────────

test_that("compare_approaches stacks both approaches on a shared key", {
  dt <- lh_test_data()
  f <- list(lh = lead_horizon(y) ~ y + x)
  setup <- setup_long_horizon(dt, f, 1:2, "gwcode", "year", test_start = 2010)
  fc <- forecast_long_horizon(setup, dt, nsim = 20, inner_sims = 3)
  lh_acc <- get_lh_accuracy(fc, dt, setup)

  sim_acc <- data.table::data.table(
    gwcode  = rep(1:6, each = 2),
    horizon = rep(1:2, 6),
    crps    = stats::runif(12),
    mae     = stats::runif(12),
    winkler = stats::runif(12)
  )
  cmp <- compare_approaches(lh_acc, sim_acc)

  expect_true(all(c("approach", "gwcode", "horizon", "crps") %in% names(cmp)))
  expect_setequal(unique(cmp$approach), c("long_horizon", "simulation"))
  expect_true(all(is.na(cmp[approach == "simulation", variant])))
})

# ── cv_long_horizon ─────────────────────────────────────────────────────────

test_that("cv_long_horizon aggregates per (variant, unit, horizon) across folds", {
  dt <- lh_test_data()
  f <- list(lh = lead_horizon(y) ~ y + x)
  cv <- cv_long_horizon(dt, f, 1:2, "gwcode", "year",
                        test_starts = c(2008, 2010), nsim = 15, inner_sims = 3)

  expect_true(all(c("variant", "gwcode", "horizon", "crps", "mae", "winkler")
                  %in% names(cv)))
  expect_true(all(is.finite(cv$crps)))
})

# ── get_accuracy additive behaviour ─────────────────────────────────────────

test_that("get_accuracy keeps per-unit default and adds a horizon mode", {
  sim <- data.table::data.table(
    gwcode = rep(1:2, each = 10),
    year   = rep(rep(2010:2011, each = 5), 2),
    y      = stats::rnorm(20),
    .sim   = rep(1:5, 4)
  )
  truth <- data.table::data.table(
    gwcode = rep(1:2, each = 2),
    year   = rep(2010:2011, 2),
    y      = stats::rnorm(4)
  )
  ctx <- panel_context("gwcode", "year")

  a_default <- get_accuracy(sim, "y", truth, ctx)
  expect_false("horizon" %in% names(a_default))
  expect_true(all(c("crps", "mae", "winkler") %in% names(a_default)))

  a_horizon <- get_accuracy(sim, "y", truth, ctx,
                            test_start = 2010, by = c("gwcode", "horizon"))
  expect_true("horizon" %in% names(a_horizon))
  expect_setequal(unique(a_horizon$horizon), c(1, 2))

  expect_error(
    get_accuracy(sim, "y", truth, ctx, by = c("gwcode", "horizon")),
    "horizon"
  )
})


# ── formula list normalisation (bare formula / unnamed) ──────────────────

test_that("setup_long_horizon accepts a bare formula and auto-names unnamed lists", {
  dt <- lh_test_data()

  # bare formula (no list wrapper)
  s_bare <- setup_long_horizon(dt, lead_horizon(y) ~ y + x, 1:2,
                               "gwcode", "year", test_start = 2010)
  expect_equal(names(s_bare$models), "model1")
  fc <- forecast_long_horizon(s_bare, dt, nsim = 5, inner_sims = 2)
  expect_gt(nrow(fc), 0)
  expect_equal(unique(fc$variant), "model1")

  # unnamed list of two -> model1, model2
  s_unnamed <- setup_long_horizon(
    dt, list(lead_horizon(y) ~ y, lead_horizon(y) ~ y + x), 1,
    "gwcode", "year", test_start = 2010)
  expect_equal(names(s_unnamed$models), c("model1", "model2"))

  # partially named -> keeps the given name, fills the blank
  s_partial <- setup_long_horizon(
    dt, list(ar = lead_horizon(y) ~ y, lead_horizon(y) ~ y + x), 1,
    "gwcode", "year", test_start = 2010)
  expect_equal(names(s_partial$models), c("ar", "model2"))
})


# ── metadata inference + multi-outcome ───────────────────────────────────

test_that("get_lh_accuracy infers metadata from forecast attributes (no setup)", {
  dt <- lh_test_data()
  f <- list(lh = lead_horizon(y) ~ y + x)
  setup <- setup_long_horizon(dt, f, 1:2, "gwcode", "year", test_start = 2010)
  fc <- forecast_long_horizon(setup, dt, nsim = 15, inner_sims = 3)
  # No lh_setup passed: groupvar/timevar/formulas come from stamped attributes.
  acc <- get_lh_accuracy(fc, dt)
  expect_true(all(c("variant", "gwcode", "horizon", "crps") %in% names(acc)))
  expect_true(all(is.finite(acc$crps)))
})

test_that("get_lh_accuracy scores each variant against its own outcome", {
  dt <- lh_test_data()
  f <- list(small = lead_horizon(y) ~ y,
            big   = lead_horizon(gdppc) ~ gdppc)
  setup <- setup_long_horizon(dt, f, 1:2, "gwcode", "year", test_start = 2010)
  fc <- forecast_long_horizon(setup, dt, nsim = 15, inner_sims = 3)
  acc <- get_lh_accuracy(fc, dt, setup)

  expect_setequal(unique(acc$variant), c("small", "big"))
  expect_true(all(is.finite(acc$crps)))
  # y ~ O(1); gdppc ~ O(1000s). The big variant's CRPS must dwarf small's,
  # which only holds if each variant is scored against its own outcome column.
  expect_gt(stats::median(acc[variant == "big", crps]),
            10 * stats::median(acc[variant == "small", crps]))
})
