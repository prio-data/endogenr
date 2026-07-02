# Tests for the composed-history defect fix (defect 1) ----------------------
#
# .required_history() must COMPOSE nested time-series depths (not take their
# max, as the legacy .max_lag_depth() does), so predict.*'s per-unit history
# window is large enough that the materialised value at t equals the value the
# same expression takes over the full series. cumulative/since-event transforms
# need the whole series and report Inf.

# ── .required_history: composition rules ───────────────────────────────────

test_that(".required_history returns the lag count, composed", {
  expect_equal(.required_history(y ~ lag(x)), 1)
  expect_equal(.required_history(y ~ lag(x, 3)), 3)
  expect_equal(.required_history(y ~ lag(x, n = 3)), 3)
  expect_equal(.required_history(y ~ lag(lag(x))), 2)        # composed, not max
  expect_equal(.required_history(y ~ lag(lag(lag(x)))), 3)
})

test_that(".required_history returns 0 with no time-series terms", {
  expect_equal(.required_history(y ~ x + I(100)), 0)
  expect_equal(.required_history(y ~ factor(region) + log(gdppc)), 0)
})

test_that(".required_history composes lag over right-aligned rolling windows", {
  # rollmeanr is right-aligned: backward reach k - 1, plus the outer lag.
  expect_equal(.required_history(y ~ lag(zoo::rollsumr(x, 5))), 5)
  expect_equal(.required_history(y ~ lag(zoo::rollsumr(x, 5), 2)), 6)
  expect_equal(.required_history(y ~ lag(zoo::rollmeanr(x, k = 3), n = 2)), 4)
  expect_equal(.required_history(y ~ zoo::rollmeanr(x, 5)), 4)   # no outer lag
})

test_that(".required_history honours rolling alignment", {
  # bare rollmean defaults to center (floor(k/2)); right adds k-1; left adds 0.
  expect_equal(.required_history(y ~ lag(zoo::rollmean(x, k = 5))), 1 + 2)
  expect_equal(
    .required_history(y ~ lag(zoo::rollmean(x, k = 5, align = "center"))), 3)
  expect_equal(
    .required_history(y ~ zoo::rollsumr(x, 5, align = "center")), 2)  # override
  expect_equal(
    .required_history(y ~ zoo::rollmean(x, k = 4, align = "left")), 0)
})

test_that(".required_history treats data.table froll* as right-aligned", {
  expect_equal(.required_history(y ~ frollmean(x, 5)), 4)
  expect_equal(.required_history(y ~ data.table::frollsum(x, 3)), 2)
})

test_that(".required_history composes diff lag x differences", {
  expect_equal(.required_history(y ~ diff(x)), 1)
  expect_equal(.required_history(y ~ diff(x, 2)), 2)               # lag = 2
  expect_equal(.required_history(y ~ diff(x, differences = 2)), 2)
  expect_equal(.required_history(y ~ diff(x, lag = 2, differences = 3)), 6)
  expect_equal(.required_history(y ~ lag(diff(x))), 2)
})

test_that(".required_history is Inf for cumulative / since-event transforms", {
  expect_identical(.required_history(y ~ cumsum(x)), Inf)
  expect_identical(.required_history(y ~ lag(cumsum(x))), Inf)
  expect_identical(.required_history(y ~ cumprod(x)), Inf)
  expect_identical(.required_history(y ~ lag(decay_since_event(z, 0.5))), Inf)
  expect_identical(.required_history(y ~ time_since_event(z)), Inf)
})

test_that(".required_history maxes across terms, descends into design ops", {
  expect_equal(.required_history(y ~ lag(x) + lag(lag(z))), 2)
  expect_equal(.required_history(y ~ poly(lag(x, 3), 2)), 3)
  expect_equal(.required_history(y ~ factor(region):lag(x, 2)), 2)
  expect_equal(.required_history(y ~ I(lag(x) + lag(lag(z)))), 2)
})

# ── .history_subset: window vs full ────────────────────────────────────────

test_that(".history_subset keeps the (t - need, t] window for finite need", {
  dt <- data.table::data.table(year = 2000:2020, v = seq_along(2000:2020))
  sub <- .history_subset(dt, "year", 2010L, 2)
  expect_equal(sub$year, c(2008L, 2009L, 2010L))         # 3 rows = need + 1
})

test_that(".history_subset floors need at 1", {
  dt <- data.table::data.table(year = 2000:2020, v = seq_along(2000:2020))
  sub <- .history_subset(dt, "year", 2010L, 0)
  expect_equal(sub$year, c(2009L, 2010L))                # floored: (t-1, t]
})

test_that(".history_subset keeps the full history up to t for Inf need", {
  dt <- data.table::data.table(year = 2000:2020, v = seq_along(2000:2020))
  sub <- .history_subset(dt, "year", 2010L, Inf)
  expect_equal(sub$year, 2000:2010)
})

test_that(".history_subset is robust to a time column named like an argument", {
  # A column literally named `t` must not shadow the `t`/`need` arguments.
  dt <- data.table::data.table(t = 1:20, v = 1:20)
  sub <- .history_subset(dt, "t", 10L, 2)
  expect_equal(sub$t, c(8L, 9L, 10L))
})

# ── predict-subset materialisation == full-series materialisation ──────────

test_that("predict-subset value equals the full-series value at t", {
  set.seed(1)
  n <- 40L
  dt <- data.table::data.table(g = 1L, t = 1:n, x = rnorm(n))
  dt$y <- rnorm(n)
  ctx <- panel_context(unit = "g", time = "t")
  t0 <- 30L

  cases <- list(
    quote(lag(zoo::rollmeanr(x, k = 5, fill = NA))),
    quote(lag(zoo::rollmean(x, k = 5, fill = NA, align = "center"))),
    quote(lag(zoo::rollmeanr(x, k = 3, fill = NA), n = 2)),
    quote(lag(lag(x))),
    quote(lag(cumsum(x)))
  )

  for (rhs in cases) {
    f <- stats::reformulate(deparse(rhs), response = "y")
    environment(f) <- environment()
    need <- .required_history(f)

    sub  <- materialize_formula(f, .history_subset(dt, "t", t0, need), ctx)
    full <- materialize_formula(f, dt[dt$t <= t0], ctx)
    pc <- setdiff(names(sub), c("g", "t", "y"))

    vsub  <- sub[sub$t == t0][[pc]]
    vfull <- full[full$t == t0][[pc]]
    # The subset must reproduce the full-series VALUE at t0. For the centered
    # window the value is NA at the boundary (needs t+1): zoo returns a logical
    # all-NA vector on the short subset but a numeric NA on the full series, so
    # compare NA-aware rather than by storage mode.
    if (is.na(vfull)) {
      expect_true(is.na(vsub), info = paste("formula:", deparse(rhs)))
    } else {
      expect_equal(vsub, vfull, info = paste("formula:", deparse(rhs)))
    }
  }
})

# ── full-simulation smoke test: no NA in forecast cells ────────────────────

test_that("simulation has no NA forecast cells for composed/cumulative lags", {
  skip_on_cran()
  future::plan(future::sequential)
  on.exit(future::plan(future::sequential), add = TRUE)

  df <- expand.grid(gwcode = c(1, 2, 3), year = 2000:2014)
  set.seed(7)
  df$x <- stats::rnorm(nrow(df))
  df$y <- 0.4 * df$x + stats::rnorm(nrow(df), sd = 0.2)
  dt <- data.table::as.data.table(df)

  # y depends on x via composed and cumulative transforms that the legacy
  # max-depth window silently truncated. x is exogenous (present at all times).
  model_system <- list(
    build_model("linear",
                formula = y ~ lag(lag(x)) +
                  lag(zoo::rollmeanr(x, k = 3, fill = NA), n = 2) +
                  lag(cumsum(x))),
    build_model("exogen", formula = ~x)
  )

  setup <- setup_system(
    models = model_system, data = dt,
    train_start = 2000, test_start = 2010, horizon = 5,
    groupvar = "gwcode", timevar = "year", inner_sims = 2
  )
  set.seed(11)
  fit <- fit_system(setup, nsim = 2)
  res <- simulate_system(fit)

  forecast_cells <- res[res$year >= 2010]
  expect_false(any(is.na(forecast_cells$y)))
  # x (exogenous) is also fully populated
  expect_false(any(is.na(forecast_cells$x)))
})
