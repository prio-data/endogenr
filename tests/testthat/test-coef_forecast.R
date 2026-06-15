# Tests for R/coef_forecast.R ----------------------------------------------

# ── Helper: balanced synthetic panel with a linear y ~ x signal ───────────
# When `trend = TRUE` the slope beta_t rises over time, so expanding/rolling
# window estimates of the coefficient genuinely trend (a real TVP signal).

cf_test_data <- function(trend = FALSE) {
  set.seed(123)
  units <- 1:10
  years <- 1990:2015
  rows <- lapply(units, function(g) {
    n <- length(years)
    x <- stats::rnorm(n, mean = 5, sd = 1)
    beta <- if (trend) 0.2 + 0.04 * (years - 1990) else rep(0.5, n)
    y <- beta * x + stats::rnorm(n, sd = 0.2)
    data.frame(gwcode = g, year = years, x = x, y = y)
  })
  data.table::as.data.table(do.call(rbind, rows))
}

cf_setup <- function(data = cf_test_data(), min_window = 8L) {
  setup_system(
    models = list(
      build_model("linear", formula = y ~ x, boot = "resid"),
      build_model("exogen", formula = ~x)
    ),
    data = data, train_start = 1990L, test_start = 2010L, horizon = 5L,
    groupvar = "gwcode", timevar = "year", inner_sims = 2L,
    min_window = min_window
  )
}

# ── .cf_anchor_grid ─────────────────────────────────────────────────────────

test_that(".cf_anchor_grid is leakage-free and well-formed", {
  g <- .cf_anchor_grid(1990L, 2010L, "expanding", min_train = 8L,
                       width = 8L, step = 1L)
  expect_equal(min(g), 1997L)            # train_start + min_train - 1
  expect_equal(max(g), 2009L)            # test_start - 1 (origin, never crossed)
  expect_true(all(diff(g) > 0))

  gr <- .cf_anchor_grid(1990L, 2010L, "rolling", min_train = 8L,
                        width = 12L, step = 5L)
  expect_equal(min(gr), 2001L)           # train_start + width - 1
  expect_equal(max(gr), 2009L)           # origin force-included despite step
  expect_true(2009L %in% gr)

  expect_error(
    .cf_anchor_grid(1990L, 2010L, "expanding", min_train = 50L,
                    width = 50L, step = 1L),
    "anchor"
  )
})

# ── .cf_rmv ─────────────────────────────────────────────────────────────────

test_that(".cf_rmv reproduces the target mean and covariance", {
  set.seed(1)
  S <- matrix(c(1, 0.6, 0.6, 2), 2L)
  X <- .cf_rmv(50000L, c(3, -1), S)
  expect_equal(colMeans(X), c(3, -1), tolerance = 0.05)
  expect_equal(stats::cov(X), S, tolerance = 0.05)
})

test_that(".cf_rmv tolerates a non-PSD covariance", {
  Sbad <- matrix(c(1, 2, 2, 1), 2L)      # eigenvalues -1, 3
  expect_silent(.cf_rmv(100L, c(0, 0), Sbad))
})

# ── forecast_coefficients: contract ──────────────────────────────────────────

test_that("forecast_coefficients returns a well-formed coef_forecast", {
  sys <- cf_setup()
  set.seed(2)
  fc <- forecast_coefficients(sys, horizon = 5L, method = "rw")

  expect_s3_class(fc, "coef_forecast")
  expect_true(all(c("outcome", "term", ".time", ".type", ".h", ".mean",
                    ".sd", ".q05", ".q50", ".q95") %in% names(fc)))
  expect_false(".draws" %in% names(fc))
  expect_setequal(unique(fc$.type), c("observed", "forecast"))
  expect_equal(attr(fc, "panel_unit"), "gwcode")
  expect_equal(attr(fc, "panel_time"), "year")
  expect_equal(attr(fc, "cf_method"), "rw")

  fcst <- fc[fc$.type == "forecast", ]
  expect_equal(sort(unique(fcst$.h)), 1:5)
  expect_setequal(unique(fcst$.time), 2010:2014)
  expect_true(all(c("(Intercept)", "x") %in% unique(fc$term)))
  # forecast band is ordered
  expect_true(all(fcst$.q05 <= fcst$.q50 & fcst$.q50 <= fcst$.q95))
})

test_that("draws attaches per-time draws; fitted systems are accepted", {
  sys <- cf_setup()
  set.seed(3)
  fc <- forecast_coefficients(sys, horizon = 3L, nsim = 200L, draws = TRUE)
  expect_true(".draws" %in% names(fc))
  fr <- which(fc$.type == "forecast")[1L]
  expect_length(fc$.draws[[fr]], 200L)
  orow <- which(fc$.type == "observed")[1L]
  expect_null(fc$.draws[[orow]])

  set.seed(123)
  fit <- fit_system(sys, nsim = 2L)
  fc2 <- forecast_coefficients(fit, horizon = 3L)
  expect_s3_class(fc2, "coef_forecast")
})

# ── forecast_coefficients: random-walk behaviour ─────────────────────────────

test_that("random-walk forecast fans out and stays centered on the origin", {
  sys <- cf_setup()
  set.seed(4)
  fc <- forecast_coefficients(sys, horizon = 6L, method = "rw", nsim = 4000L)
  fcst <- fc[fc$.type == "forecast", ]

  for (tm in unique(fcst$term)) {
    s <- fcst[fcst$term == tm, ]
    s <- s[order(s$.h), ]
    expect_gt(s$.sd[nrow(s)], s$.sd[1L])           # variance grows with horizon
    origin_mean <- fc$.mean[fc$.type == "observed" &
                              fc$term == tm & fc$.time == 2009L]
    expect_lt(max(abs(s$.mean - origin_mean)), 0.5 * max(s$.sd))  # no drift
  }
})

# ── drift recovers the path slope (internal forecaster) ──────────────────────

test_that("drift continues the path slope while rw stays flat", {
  set.seed(99)
  od <- data.frame(x = stats::rnorm(50), y = stats::rnorm(50))
  origin_fit <- stats::lm(y ~ x, od)     # supplies coef names / vcov / df only
  taus <- 1:12
  path <- data.table::rbindlist(lapply(taus, function(t) {
    data.table::data.table(
      outcome = "y", term = c("(Intercept)", "x"), .tau = t,
      estimate = c(2, 0.5 + 0.1 * t),    # 'x' trends at +0.1 per step
      std.error = c(0.01, 0.01)
    )
  }))

  set.seed(5)
  d <- .cf_forecast_one("y", path, origin_fit, test_start = 13L, horizon = 5L,
                        method = "drift", nsim = 4000L, student_t = FALSE, step = 1L)
  fx_d <- d$fc[d$fc$term == "x", ]
  fx_d <- fx_d[order(fx_d$.h), ]
  expect_equal(mean(diff(fx_d$.mean)), 0.1, tolerance = 0.03)   # slope recovered

  r <- .cf_forecast_one("y", path, origin_fit, test_start = 13L, horizon = 5L,
                        method = "rw", nsim = 4000L, student_t = FALSE, step = 1L)
  fx_r <- r$fc[r$fc$term == "x", ]
  fx_r <- fx_r[order(fx_r$.h), ]
  expect_lt(abs(mean(diff(fx_r$.mean))), 0.03)                  # rw is flat
})

# ── errors ───────────────────────────────────────────────────────────────────

test_that("forecast_coefficients errors clearly on bad inputs", {
  sys <- cf_setup()
  expect_error(forecast_coefficients(sys, min_train = 100L), "anchor")  # no windows
  expect_error(forecast_coefficients(sys, min_train = 19L), "anchor")   # only 2 anchors

  no_lin <- setup_system(
    models = list(
      build_model("deterministic", formula = y ~ I(lag(y) + x)),
      build_model("exogen", formula = ~x)
    ),
    data = cf_test_data(), train_start = 1990L, test_start = 2010L, horizon = 5L,
    groupvar = "gwcode", timevar = "year", inner_sims = 2L
  )
  expect_error(forecast_coefficients(no_lin), "linear")

  expect_error(forecast_coefficients(list(a = 1)),
               "endogenr_system_setup")
})

# ── rolling windows + plotting ───────────────────────────────────────────────

test_that("rolling windows and plotting produce expected objects", {
  sys <- cf_setup()
  set.seed(6)
  fr <- forecast_coefficients(sys, horizon = 4L, window = "rolling", width = 12L)
  expect_s3_class(fr, "coef_forecast")
  expect_equal(attr(fr, "cf_window"), "rolling")

  p <- plot_coefficient_forecast(fr)
  expect_s3_class(p, "ggplot")

  p2 <- plot_coefficient_forecast(sys, method = "rw", horizon = 3L)
  expect_s3_class(p2, "ggplot")
})
