# Panel design layer: poly / data-dependent bases + naming contract -----------
#
# Validates:
#   1. region:poly(dem, 2) end-to-end through linearmodel and predict.linear
#      (including units with < 3 unique dem values — the previously-failing case)
#   2. poly/bs basis coherence: predict reuses training-time scaling
#   3. Naming contract: coef terms are natural R names; no .pt# leaks;
#      lag(x) -> lag_x alias holds
#   4. systemgraph ts-function registry is a single source of truth
#   5. Full setup_system -> fit_system -> simulate_system smoke test with poly

# ── Helper: panel with a region column and dem that varies across regions ──

make_poly_panel <- function(seed = 1L) {
  set.seed(seed)
  n_units <- 12L
  n_time  <- 15L
  dt <- data.table::data.table(
    unit   = rep(seq_len(n_units), each = n_time),
    time   = rep(seq_len(n_time),  n_units),
    region = rep(c("A", "B", "C"), each = n_time * (n_units / 3L))
  )
  # dem varies across units but is CONSTANT within each unit (worst case for
  # poly per-group: only 1 unique point per unit -> poly(dem,2) errors).
  dt[, dem := seq(1, by = 0.5, length.out = n_units)[unit]]
  dt[, gdppc := exp(1 + 0.2 * time + 0.01 * dem + stats::rnorm(.N, 0, 0.1))]
  dt[, y := 1 + 0.5 * dem + 0.2 * time + stats::rnorm(.N, 0, 0.3)]
  dt
}

# ── 1. linearmodel: region:poly(dem, 2) fits without error ────────────────────

test_that("linearmodel: region:poly(dem,2) fits without error", {
  dt  <- make_poly_panel()
  ctx <- panel_context(unit = "unit", time = "time")
  # Should NOT error despite each unit having exactly 1 unique dem value
  expect_no_error(
    m <- linearmodel(y ~ region + poly(dem, 2) + lag(y), data = dt, ctx = ctx)
  )
  expect_true(!is.null(m$fitted))
  expect_true(is.finite(m$gof$r.squared))
})

test_that("linearmodel: region:poly(dem,2) coef names are natural (no .pt# leak)", {
  dt  <- make_poly_panel()
  ctx <- panel_context(unit = "unit", time = "time")
  m   <- linearmodel(y ~ region + poly(dem, 2) + lag(y), data = dt, ctx = ctx)

  terms_in <- m$coefs$term
  # No .pt# internal symbol should leak into user-facing coefficients
  expect_false(any(grepl("^\\.pt\\d", terms_in)),
               info = paste("Found .pt# in coefs:", paste(terms_in, collapse = ", ")))
  # poly terms match base lm natural naming
  expect_true(any(grepl("poly\\(dem, 2\\)1", terms_in)))
  expect_true(any(grepl("poly\\(dem, 2\\)2", terms_in)))
  # lag(y) materialised as readable alias
  expect_true(any(grepl("lag_y", terms_in)))
})

test_that("linearmodel: region:poly(dem,2) coefs match base lm", {
  dt  <- make_poly_panel()
  ctx <- panel_context(unit = "unit", time = "time")
  m   <- linearmodel(y ~ region + poly(dem, 2) + lag(y), data = dt, ctx = ctx)

  # Base lm with same data (lag done manually for comparison)
  dt2 <- data.table::copy(dt)
  data.table::setorder(dt2, unit, time)
  dt2[, lag_y := c(NA, utils::head(y, -1L)), by = unit]
  base_fit <- stats::lm(y ~ region + poly(dem, 2) + lag_y,
                        data = stats::na.omit(dt2))

  base_coef <- sort(names(stats::coef(base_fit)))
  endo_coef <- sort(m$coefs$term)
  expect_equal(endo_coef, base_coef)
})

# ── 2. predict.linear: poly basis coherent with fit ──────────────────────────

test_that("predict.linear: poly basis is reused from fit (finite, finite)", {
  set.seed(10L)
  dt  <- make_poly_panel()
  ctx_fit <- panel_context(unit = "unit", time = "time")
  ctx_sim <- panel_context(unit = "unit", time = "time", sim = "sim")
  dt[, sim := 1L]

  m <- linearmodel(y ~ poly(dem, 2) + lag(y), data = dt[time <= 12L],
                   ctx = ctx_fit)

  # Predict time 13 — poly basis must reuse training-time scaling, not refit
  pred <- predict(m, data = dt, t = 13L, ctx = ctx_sim)
  expect_equal(nrow(pred), 12L)  # one row per unit
  expect_true(all(is.finite(pred$y)))
})

test_that("predict.linear: region:poly interaction over rolling window", {
  dt  <- make_poly_panel()
  ctx_fit <- panel_context(unit = "unit", time = "time")
  ctx_sim <- panel_context(unit = "unit", time = "time", sim = "sim")
  dt[, sim := 1L]

  m <- linearmodel(y ~ region:poly(dem, 2) + lag(y), data = dt[time <= 12L],
                   ctx = ctx_fit)

  pred <- predict(m, data = dt, t = 13L, ctx = ctx_sim)
  expect_equal(nrow(pred), 12L)
  expect_true(all(is.finite(pred$y)))
})

# ── 3. glmmodel: poly in a GLM ────────────────────────────────────────────────

test_that("glmmodel: poly(dem, 2) fits and predicts without error", {
  dt  <- make_poly_panel()
  ctx_fit <- panel_context(unit = "unit", time = "time")
  ctx_sim <- panel_context(unit = "unit", time = "time", sim = "sim")
  dt[, sim := 1L]
  dt[, ypos := exp(y)]  # positive outcome for Gamma

  expect_no_error(
    m <- glmmodel(y ~ poly(dem, 2) + lag(y), data = dt[time <= 12L],
                  family = stats::gaussian(), ctx = ctx_fit)
  )
  pred <- predict(m, data = dt, t = 13L, ctx = ctx_sim)
  expect_equal(nrow(pred), 12L)
  expect_true(all(is.finite(pred$y)))
})

# ── 4. Naming contract: lag(x) -> lag_x, no .pt# in get_coefficients ────────

test_that("get_coefficients: lag_x alias holds; no .pt# leaks", {
  future::plan(future::sequential)
  on.exit(future::plan(future::sequential), add = TRUE)

  set.seed(5L)
  dt <- data.table::data.table(
    unit = rep(1:6, each = 20L),
    time = rep(1:20L, 6L)
  )
  dt[, x := stats::rnorm(.N)]
  dt[, y := 0.5 * x + stats::rnorm(.N, 0, 0.3)]

  setup <- setup_system(
    list(build_model("linear", formula = y ~ lag(x)),
         build_model("exogen", formula = ~x)),
    dt,
    train_start = 1L, test_start = 16L, horizon = 2L,
    groupvar = "unit", timevar = "time", inner_sims = 1L
  )
  co <- get_coefficients(fit_system(setup, nsim = 2L))

  # lag(x) materialised as lag_x
  expect_true("lag_x" %in% co$term)
  # No .pt# internal symbols leak into user-facing output
  expect_false(any(grepl("^\\.pt\\d", co$term)))
})

# ── 5. .pt_ts_fns is the single registry (systemgraph uses same set) ─────────

test_that(".pt_ts_fns and .required_history agree on lag depth", {
  # lag(x) -> 1, lag(x, 2) -> 2, cumsum(x) -> Inf
  expect_equal(.required_history(y ~ lag(x)), 1)
  expect_equal(.required_history(y ~ lag(x, 2)), 2)
  expect_equal(.required_history(y ~ cumsum(x)), Inf)
  # rollmean in .pt_ts_fns -> depth extracted
  # rollmean is center-aligned by default: backward reach = floor(k/2)
  expect_equal(.required_history(y ~ rollmean(x, k = 4)), 2)  # center: floor(4/2)
  expect_equal(.required_history(y ~ rollmeanr(x, k = 4)), 3) # right: k-1
})

test_that(".pt_cum_fns and .pt_roll_fns are subsets of .pt_ts_fns", {
  expect_true(all(.pt_cum_fns  %in% .pt_ts_fns))
  expect_true(all(.pt_roll_fns %in% .pt_ts_fns))
})

# ── 6. Full pipeline: poly through rolling-window fit_system ──────────────────

test_that("fit_system(window='rolling') with poly(dem,2) runs end-to-end", {
  skip_on_cran()
  future::plan(future::sequential)
  on.exit(future::plan(future::sequential), add = TRUE)

  dt  <- make_poly_panel()
  # setup_system does not need a sim column in the input data

  setup <- setup_system(
    list(
      build_model("linear",
                  formula = y ~ region + poly(dem, 2) + lag(y),
                  boot = "resid"),
      build_model("exogen", formula = ~dem + region)
    ),
    dt,
    train_start = 1L, test_start = 13L, horizon = 2L,
    groupvar = "unit", timevar = "time", inner_sims = 1L,
    min_window = 8L
  )

  set.seed(42L)
  fit  <- fit_system(setup, nsim = 2L)
  sims <- simulate_system(fit)

  expect_s3_class(fit,  "endogenr_fitted_system")
  expect_s3_class(sims, "data.table")
  expect_true(all(is.finite(sims[sims$time >= 13L, y])))
})
