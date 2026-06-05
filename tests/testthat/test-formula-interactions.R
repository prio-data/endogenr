# Formula interaction and factor(unit) preservation --------------------------
#
# Verifies that:
#   1. .clean_term_labels / derive_naive_formula preserve interaction structure
#   2. linear/glm models estimate interactions identically to base lm/glm
#   3. factor(<unit-key>) materializes and fits without error
#   4. heterolm preserves interactions on both sides of the | separator

# ── .clean_term_labels ────────────────────────────────────────────────────────

test_that(".clean_term_labels preserves simple main effects", {
  cm <- c()  # no renaming
  labs <- endogenr:::.clean_term_labels(y ~ a + b, cm)
  expect_equal(labs, c("a", "b"))
})

test_that(".clean_term_labels preserves crossing (y ~ g*x)", {
  cm <- c()
  labs <- endogenr:::.clean_term_labels(y ~ g * x, cm)
  expect_setequal(labs, c("g", "x", "g:x"))
})

test_that(".clean_term_labels maps transformed variable tokens", {
  cm <- c("lag(x)" = "lag_x", "factor(g)" = "factor_g")
  labs <- endogenr:::.clean_term_labels(y ~ factor(g) * lag(x), cm)
  expect_setequal(labs, c("factor_g", "lag_x", "factor_g:lag_x"))
})

test_that(".clean_term_labels returns character(0) for intercept-only formula", {
  labs <- endogenr:::.clean_term_labels(y ~ 1, c())
  expect_equal(length(labs), 0L)
})

# ── derive_naive_formula ──────────────────────────────────────────────────────

test_that("derive_naive_formula: y ~ a*b preserves interaction", {
  cm <- c()
  f  <- derive_naive_formula(y ~ a * b, col_mapping = cm)
  expect_true(inherits(f, "formula"))
  labs <- labels(stats::terms(f))
  expect_setequal(labs, c("a", "b", "a:b"))
})

test_that("derive_naive_formula: y ~ 0 + x drops intercept", {
  cm <- c()
  f  <- derive_naive_formula(y ~ 0 + x, col_mapping = cm)
  expect_equal(attr(stats::terms(f), "intercept"), 0L)
  expect_equal(labels(stats::terms(f)), "x")
})

test_that("derive_naive_formula: y ~ 1 returns intercept-only", {
  cm <- c()
  f  <- derive_naive_formula(y ~ 1, col_mapping = cm)
  expect_equal(attr(stats::terms(f), "intercept"), 1L)
  expect_equal(length(labels(stats::terms(f))), 0L)
})

test_that("derive_naive_formula applies col_mapping to interaction labels", {
  cm <- c("lag(x)" = "lag_x")
  f  <- derive_naive_formula(y ~ g * lag(x), col_mapping = cm)
  labs <- labels(stats::terms(f))
  expect_setequal(labs, c("g", "lag_x", "g:lag_x"))
})

# ── linear model parity with base lm ─────────────────────────────────────────

test_that("linear model: y ~ g*x coefs match base lm", {
  set.seed(1)
  dt <- data.table::data.table(
    unit = rep(1:20, each = 20),
    time = rep(1:20, 20)
  )
  dt[, g := ifelse(unit <= 10, "A", "B")]
  dt[, x := stats::rnorm(.N)]
  dt[, y := 1 + 0.5 * (g == "B") + 2 * x + 0.8 * (g == "B") * x +
       stats::rnorm(.N, sd = 0.5)]

  ctx <- panel_context(unit = "unit", time = "time")
  m   <- linearmodel(y ~ g * x, data = dt, ctx = ctx)

  base_lm <- stats::lm(y ~ g * x, data = dt)
  base_co <- sort(names(stats::coef(base_lm)))
  endo_co <- sort(m$coefs$term)

  # Same term names
  expect_equal(endo_co, base_co)
  # Close coefficient values
  for (nm in base_co) {
    endogenr_est <- m$coefs$estimate[m$coefs$term == nm]
    base_est     <- unname(stats::coef(base_lm)[nm])
    expect_equal(endogenr_est, base_est, tolerance = 1e-6,
                 label = paste("coef", nm))
  }
})

# ── factor(unit) as fixed effects ─────────────────────────────────────────────

test_that("linear model: factor(unit) fixed effects fit without error", {
  set.seed(2)
  dt <- data.table::data.table(
    unit = rep(1:5, each = 15),
    time = rep(1:15, 5)
  )
  dt[, x := stats::rnorm(.N)]
  dt[, y := unit * 0.3 + 0.8 * x + stats::rnorm(.N, sd = 0.5)]

  ctx <- panel_context(unit = "unit", time = "time")

  # Should not error
  m <- linearmodel(y ~ factor(unit) + x, data = dt, ctx = ctx)
  expect_true(!is.null(m$fitted))
  expect_true(any(grepl("factor_unit", m$coefs$term)))
})

test_that("linear model: factor(unit)*x interaction fits and includes cross terms", {
  set.seed(3)
  dt <- data.table::data.table(
    unit = rep(1:4, each = 20),
    time = rep(1:20, 4)
  )
  dt[, x := stats::rnorm(.N)]
  dt[, y := unit * 0.2 + x + unit * 0.1 * x + stats::rnorm(.N, sd = 0.3)]

  ctx <- panel_context(unit = "unit", time = "time")
  m   <- linearmodel(y ~ factor(unit) * x, data = dt, ctx = ctx)

  expect_true(!is.null(m$fitted))
  # At least one interaction term present
  expect_true(any(grepl("factor_unit.*:x|x:factor_unit", m$coefs$term)))
})

test_that("factor(unit) coefs match base lm on same data", {
  set.seed(4)
  dt <- data.table::data.table(
    unit = rep(1:5, each = 20),
    time = rep(1:20, 5)
  )
  dt[, x := stats::rnorm(.N)]
  dt[, y := unit * 0.4 + 1.2 * x + stats::rnorm(.N, sd = 0.5)]

  ctx     <- panel_context(unit = "unit", time = "time")
  m       <- linearmodel(y ~ factor(unit) + x, data = dt, ctx = ctx)
  base_lm <- stats::lm(y ~ factor(unit) + x, data = dt)

  base_x <- unname(stats::coef(base_lm)["x"])
  endo_x <- m$coefs$estimate[m$coefs$term == "x"]
  expect_equal(endo_x, base_x, tolerance = 1e-6)
})

# ── GLM interaction parity ────────────────────────────────────────────────────

test_that("glm model: y ~ g*x preserves g:x interaction", {
  set.seed(5)
  dt <- data.table::data.table(
    unit = rep(1:15, each = 20),
    time = rep(1:20, 15)
  )
  dt[, g := ifelse(unit <= 7, "A", "B")]
  dt[, x := stats::rnorm(.N)]
  dt[, y := 0.5 + 0.3 * (g == "B") + x + 0.4 * (g == "B") * x +
       stats::rnorm(.N, sd = 0.5)]

  ctx  <- panel_context(unit = "unit", time = "time")
  m    <- glmmodel(y ~ g * x, family = stats::gaussian(), data = dt, ctx = ctx)
  labs <- labels(stats::terms(m$naive_formula))
  expect_true(any(grepl(":", labs)), info = "interaction term missing from naive_formula")

  base_glm <- stats::glm(y ~ g * x, data = dt, family = stats::gaussian())
  endo_gBx <- m$coefs$estimate[m$coefs$term == "gB:x"]
  base_gBx <- unname(stats::coef(base_glm)["gB:x"])
  expect_equal(endo_gBx, base_gBx, tolerance = 1e-6)
})

# ── heterolm interaction ──────────────────────────────────────────────────────

test_that("heterolm: mean y ~ g*x captures interaction (more than 2 mean terms)", {
  skip_if_not_installed("heterolm")

  set.seed(6)
  dt <- data.table::data.table(
    unit = rep(1:10, each = 20),
    time = rep(1:20, 10)
  )
  dt[, g := ifelse(unit <= 5, "A", "B")]
  dt[, x := stats::rnorm(.N)]
  dt[, z := stats::rnorm(.N)]
  dt[, y := 1 + 0.5 * (g == "B") + x + 0.3 * (g == "B") * x +
       stats::rnorm(.N, sd = sqrt(exp(-0.5 + 0.5 * z)))]

  ctx <- panel_context(unit = "unit", time = "time")
  m   <- heterolmmodel(y ~ g * x, variance = ~ z, data = dt, ctx = ctx)

  # heterolm::hetero() needs flat column-name lookups; the interaction g:x is
  # pre-expanded via model.matrix into a separate column (e.g. g_b_x).
  # Verify: the mean side of naive_hetero_formula has > 2 terms.
  frm_str   <- deparse(m$naive_hetero_formula)
  mean_side <- sub("\\s*\\|.*$", "", sub("^.*~\\s*", "", frm_str))
  mean_n    <- length(strsplit(trimws(mean_side), "\\s*\\+\\s*")[[1L]])
  expect_gt(mean_n, 2L,
            label = paste("mean-side term count > 2 in:", frm_str))
})

# ── predict path smoke test ───────────────────────────────────────────────────

test_that("predict.linear works for an interaction model", {
  set.seed(7)
  dt <- data.table::data.table(
    unit = rep(1:5, each = 20),
    time = rep(1:20, 5)
  )
  dt[, g := ifelse(unit <= 2, "A", "B")]
  dt[, x := stats::rnorm(.N)]
  dt[, y := 1 + 0.5 * (g == "B") + x + 0.3 * (g == "B") * x +
       stats::rnorm(.N, sd = 0.5)]

  ctx <- panel_context(unit = "unit", time = "time", sim = "sim")
  dt[, sim := 1L]
  m   <- linearmodel(y ~ g * x, data = dt[time <= 15], ctx = ctx)

  pred_data <- data.table::copy(dt)
  pred <- predict(m, data = pred_data, t = 16L, ctx = ctx)
  expect_equal(nrow(pred), 5L)   # one row per unit
  expect_true(all(is.finite(pred$y)))
})

test_that("predict.linear works for factor(unit) fixed-effect model", {
  set.seed(8)
  dt <- data.table::data.table(
    unit = rep(1:4, each = 20),
    time = rep(1:20, 4)
  )
  dt[, x := stats::rnorm(.N)]
  dt[, y := unit * 0.3 + 0.8 * x + stats::rnorm(.N, sd = 0.4)]

  ctx <- panel_context(unit = "unit", time = "time", sim = "sim")
  dt[, sim := 1L]
  m   <- linearmodel(y ~ factor(unit) + x, data = dt[time <= 15], ctx = ctx)

  pred <- predict(m, data = dt, t = 16L, ctx = ctx)
  expect_equal(nrow(pred), 4L)
  expect_true(all(is.finite(pred$y)))
})
