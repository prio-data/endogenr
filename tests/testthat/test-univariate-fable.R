# Tests for the fable path-coherence defect fix (defect 2) ------------------
#
# predict.univariate_fable must draw COHERENT sample paths via
# fabletools::generate() on the fitted mable (correlated innovations simulated
# forward), not stitch independent per-horizon draws from the marginal
# forecast() distributions. The latter destroys the temporal structure of the
# ETS/ARIMA path before it feeds (via lag()) into the endogenous equations.

# Strongly-trending balanced panel: a fitted ETS(A,A,N) should yield smooth,
# highly autocorrelated coherent paths.
fable_panel <- function(units = c(1, 2, 3), t_n = 40, seed = 42) {
  set.seed(seed)
  rows <- lapply(units, function(u) {
    y <- cumsum(0.8 + stats::rnorm(t_n, sd = 0.4))
    data.table::data.table(gwcode = u, year = seq_len(t_n), dem = y)
  })
  data.table::rbindlist(rows)
}

fit_fable <- function(dt) {
  ctx <- panel_context(unit = "gwcode", time = "year", sim = "sim")
  list(
    mod = univariate_fable_model(
      dem ~ error("A") + trend("A") + season("N"),
      data = dt, method = "ets", ctx = ctx),
    ctx = ctx
  )
}

# Mean within-path lag-1 autocorrelation of DEVIATIONS from the per-horizon
# ensemble mean (removes the shared deterministic trend, isolating innovation
# coherence). ~0 for independent per-horizon draws, strongly positive for
# coherent paths.
.dev_acf <- function(pred) {
  p <- data.table::as.data.table(pred)
  p[, dv := dem - mean(dem), by = c("gwcode", "year")]
  data.table::setorderv(p, c("gwcode", "sim", "year"))
  a <- p[, {
    v <- dv
    if (length(v) >= 3 && stats::sd(v) > 0) {
      list(a = stats::cor(v[-length(v)], v[-1]))
    } else list(a = NA_real_)
  }, by = c("gwcode", "sim")]
  mean(a$a, na.rm = TRUE)
}

.raw_acf <- function(pred) {
  p <- data.table::as.data.table(pred)
  data.table::setorderv(p, c("gwcode", "sim", "year"))
  a <- p[, {
    v <- dem
    if (length(v) >= 3 && stats::sd(v) > 0) {
      list(a = stats::cor(v[-length(v)], v[-1]))
    } else list(a = NA_real_)
  }, by = c("gwcode", "sim")]
  mean(a$a, na.rm = TRUE)
}

# ── deterministic mapping: predict() IS a relabelling of generate() ─────────

test_that("predict.univariate_fable returns generate() paths, mapped correctly", {
  skip_if_no_fable()
  dt <- fable_panel()
  ff <- fit_fable(dt)

  # The fable model forecasts from the last observed period (year 40 here);
  # `test_start` is accepted but unused by the fable predict.
  set.seed(99)
  pred <- predict(ff$mod, data = dt, ctx = ff$ctx,
                  test_start = 35, horizon = 6, inner_sims = 20)
  origin <- max(dt$year)

  # Column contract: unit, time, sim, outcome — and the sim/time grid is full.
  expect_setequal(names(pred), c("gwcode", "year", "sim", "dem"))
  expect_setequal(sort(unique(pred$sim)), 1:20)
  expect_setequal(sort(unique(pred$year)), (origin + 1):(origin + 6))
  expect_equal(nrow(pred), 3L * 6L * 20L)
  expect_false(anyNA(pred$dem))

  # The values ARE fabletools::generate() under the same seed (this is what
  # the old forecast()+distributional::generate() stitch did NOT produce).
  set.seed(99)
  manual <- fabletools::generate(ff$mod$fitted, h = 6, times = 20) |>
    dplyr::as_tibble() |>
    dplyr::transmute(gwcode, year,
                     sim = as.integer(.rep), dem = .sim) |>
    data.table::as.data.table()
  data.table::setkeyv(pred, c("gwcode", "sim", "year"))
  data.table::setkeyv(manual, c("gwcode", "sim", "year"))
  expect_equal(pred[, c("gwcode", "sim", "year", "dem")],
               manual[, c("gwcode", "sim", "year", "dem")])
})

test_that("predict.univariate_fable produces one coherent path per (unit, sim)", {
  skip_if_no_fable()
  # Within a single path, consecutive horizon values are highly autocorrelated
  # (a single ETS draw is smooth). Independent per-horizon stitching is not.
  dt <- fable_panel()
  ff <- fit_fable(dt)
  set.seed(3)
  pred <- predict(ff$mod, data = dt, ctx = ff$ctx,
                  test_start = 29, horizon = 12, inner_sims = 60)
  expect_gt(.raw_acf(pred), 0.6)
})

# ── statistical: coherence flips vs the old index-stitching ────────────────

test_that("coherent paths beat independent-draw stitching on innovation ACF", {
  skip_on_cran()
  skip_if_no_fable()
  skip_if_not_slow()

  dt <- fable_panel()
  ff <- fit_fable(dt)
  h <- 12L
  s <- 300L

  # OLD behaviour: per-horizon marginal forecast() then independent generate().
  set.seed(1)
  old <- ff$mod$fitted |>
    fabletools::forecast(h = h) |>
    dplyr::mutate(.s = distributional::generate(dem, s)) |>
    dplyr::as_tibble() |>
    dplyr::select(gwcode, year, .s) |>
    tidyr::unnest(.s) |>
    dplyr::mutate(sim = rep(1:s, dplyr::n() / s)) |>
    dplyr::rename(dem = .s) |>
    data.table::as.data.table()

  # NEW behaviour: coherent paths from the fitted mable.
  set.seed(1)
  new <- predict(ff$mod, data = dt, ctx = ff$ctx,
                 test_start = 29, horizon = h, inner_sims = s)

  old_acf <- .dev_acf(old)
  new_acf <- .dev_acf(new)

  expect_lt(old_acf, 0.15)             # old stitch ~ independent
  expect_gt(new_acf, 0.4)              # new paths carry innovation structure
  expect_gt(new_acf - old_acf, 0.3)    # the flip, with margin
})

test_that("generated paths preserve the per-horizon marginal moments", {
  skip_on_cran()
  skip_if_no_fable()
  skip_if_not_slow()

  dt <- fable_panel()
  ff <- fit_fable(dt)
  h <- 12L
  s <- 400L

  set.seed(5)
  new <- predict(ff$mod, data = dt, ctx = ff$ctx,
                 test_start = 29, horizon = h, inner_sims = s)

  fc <- ff$mod$fitted |>
    fabletools::forecast(h = h) |>
    dplyr::as_tibble()
  fc$.sd <- sqrt(distributional::variance(fc$dem))
  fc <- data.table::as.data.table(fc)[, c("gwcode", "year", ".mean", ".sd")]

  nm <- data.table::as.data.table(new)[, list(m = mean(dem), s = stats::sd(dem)),
                                       by = c("gwcode", "year")]
  mg <- merge(nm, fc, by = c("gwcode", "year"))

  expect_lt(max(abs((mg$m - mg$.mean) / mg$.mean)), 0.05)
  expect_lt(max(abs((mg$s - mg$.sd) / mg$.sd)), 0.20)
})
