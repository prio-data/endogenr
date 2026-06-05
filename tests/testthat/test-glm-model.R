# Tests for the GLM response-scale dispersion fix (defect 3) ----------------
#
# getpi_glm() must add the family's RESPONSE-scale dispersion, not just the
# link-scale parameter uncertainty. A gaussian GLM then matches the equivalent
# lm PI width; count/positive/proportion families produce realistic draws.

glm_data <- function(seed = 1, n = 400) {
  set.seed(seed)
  x <- stats::rnorm(n)
  data.table::data.table(
    g = 1L, t = seq_len(n), x = x,
    y      = 2 + 1.5 * x + stats::rnorm(n, sd = 2),                 # gaussian
    ycount = stats::rpois(n, lambda = exp(0.5 + 0.3 * x)),          # poisson
    ypos   = stats::rgamma(n, shape = 2, scale = exp(0.2 + 0.1 * x) / 2), # Gamma
    yprop  = stats::rbinom(n, 1, stats::plogis(-0.2 + 0.5 * x))     # proportion
  )
}

q_width <- function(v, lo = 0.05, hi = 0.95) {
  unname(diff(stats::quantile(v, c(lo, hi))))
}

# ── Gaussian GLM PI width matches lm; link-only under-disperses (the flip) ──

test_that("gaussian GLM PI width matches the lm PI width", {
  dt <- glm_data()
  nd <- dt[t == nrow(dt)]
  m  <- 5000L

  lm_fit  <- stats::lm(y ~ x, dt)
  lm_pred <- stats::predict(lm_fit, newdata = nd, se.fit = TRUE)
  set.seed(2)
  lm_w <- q_width(getpi(lm_pred, nsamples = m))

  glm_g <- stats::glm(y ~ x, dt, family = stats::gaussian())
  disp  <- summary(glm_g)$dispersion
  gp    <- stats::predict(glm_g, newdata = nd, type = "link", se.fit = TRUE)
  set.seed(2)
  glm_w <- q_width(getpi_glm(gp, glm_g$family, glm_g$df.residual, disp, nsamples = m))

  # Pre-fix behaviour: link-scale parameter draw only (no residual dispersion).
  set.seed(2)
  link_only <- as.vector(glm_g$family$linkinv(
    gp$fit + outer(gp$se.fit, stats::rt(m, glm_g$df.residual))))
  old_w <- q_width(link_only)

  expect_gt(glm_w / lm_w, 0.8)         # parity restored
  expect_lt(glm_w / lm_w, 1.25)
  expect_lt(old_w / lm_w, 0.5)         # the defect: old PI far too narrow
  expect_lt(old_w, glm_w)              # response draw widens the interval
})

# ── Structural per-family contracts (deterministic, always-on) ─────────────

test_that("poisson PI draws are non-negative integers", {
  dt <- glm_data()
  nd <- dt[t == nrow(dt)]
  glm_p <- stats::glm(ycount ~ x, dt, family = stats::poisson())
  pp <- stats::predict(glm_p, newdata = nd, type = "link", se.fit = TRUE)
  set.seed(3)
  draws <- getpi_glm(pp, glm_p$family, glm_p$df.residual, 1, nsamples = 500)
  expect_true(all(draws == round(draws)))
  expect_true(all(draws >= 0))
})

test_that("Gamma PI draws are strictly positive", {
  dt <- glm_data()
  nd <- dt[t == nrow(dt)]
  glm_gam <- stats::glm(ypos ~ x, dt, family = stats::Gamma(link = "log"))
  disp <- summary(glm_gam)$dispersion
  gmp <- stats::predict(glm_gam, newdata = nd, type = "link", se.fit = TRUE)
  set.seed(4)
  draws <- getpi_glm(gmp, glm_gam$family, glm_gam$df.residual, disp, nsamples = 500)
  expect_true(all(draws > 0))
})

test_that("binomial (proportion) PI draws stay in [0, 1]", {
  dt <- glm_data()
  nd <- dt[t == nrow(dt)]
  glm_b <- stats::glm(yprop ~ x, dt, family = stats::binomial())
  bp <- stats::predict(glm_b, newdata = nd, type = "link", se.fit = TRUE)
  set.seed(5)
  draws <- getpi_glm(bp, glm_b$family, glm_b$df.residual, 1, nsamples = 500)
  expect_true(all(draws >= 0 & draws <= 1))
})

test_that("unsupported families warn once and fall back to a link-only draw", {
  dt <- glm_data()
  nd <- dt[t == nrow(dt)]
  glm_ig <- stats::glm(ypos ~ x, dt,
                       family = stats::inverse.gaussian(link = "log"))
  ig <- stats::predict(glm_ig, newdata = nd, type = "link", se.fit = TRUE)
  disp <- summary(glm_ig)$dispersion

  options(.endogenr_glm_warned_inverse.gaussian = NULL)
  set.seed(6)
  expect_warning(
    d1 <- getpi_glm(ig, glm_ig$family, glm_ig$df.residual, disp, nsamples = 50),
    "no response-scale predictive draw")
  # one-time-per-session: a second call is silent
  expect_silent(getpi_glm(ig, glm_ig$family, glm_ig$df.residual, disp, nsamples = 50))
  # fallback equals the link-scale mean (parameter uncertainty only)
  expect_true(all(d1 > 0))
})

test_that("glmmodel caches the response-scale dispersion", {
  dt <- glm_data()
  fit_ctx <- panel_context(unit = "g", time = "t")
  mg <- glmmodel(y ~ x, family = stats::gaussian(), data = dt, ctx = fit_ctx)
  mp <- glmmodel(ycount ~ x, family = stats::poisson(), data = dt, ctx = fit_ctx)
  expect_gt(mg$dispersion, 0)                 # estimated for gaussian
  expect_equal(mp$dispersion, 1)              # fixed for poisson
})

# ── what = "expectation" is unchanged (pinned) ─────────────────────────────

test_that("predict.glm_endogenr expectation is the deterministic mean", {
  df <- data.table::as.data.table(expand.grid(g = c(1, 2), year = 1:40))
  set.seed(1)
  df$x <- stats::rnorm(nrow(df))
  df$y <- stats::rpois(nrow(df), exp(0.3 + 0.4 * df$x))
  fit_ctx <- panel_context(unit = "g", time = "year")
  sim_ctx <- panel_context(unit = "g", time = "year", sim = "sim")
  mod <- glmmodel(y ~ x, family = stats::poisson(), data = df, ctx = fit_ctx)

  grid <- data.table::copy(df)
  grid[, sim := 1L]

  e1 <- predict(mod, data = grid, t = 40L, ctx = sim_ctx, what = "expectation")
  e2 <- predict(mod, data = grid, t = 40L, ctx = sim_ctx, what = "expectation")
  expect_equal(e1$y, e2$y)                    # deterministic, no draw

  # equals linkinv(linear predictor)
  nd <- materialize_formula(mod$formula, grid[year == 40L], sim_ctx,
                            mod$mat_formula, mod$col_mapping)
  mu <- stats::predict(mod$fitted, newdata = nd, type = "response")
  expect_equal(unname(e1$y), unname(as.numeric(mu)))

  # pi draws differ from the expectation (carry randomness)
  set.seed(9)
  p1 <- predict(mod, data = grid, t = 40L, ctx = sim_ctx, what = "pi")
  expect_false(isTRUE(all.equal(p1$y, e1$y)))
  expect_true(all(p1$y == round(p1$y)))       # integer counts end-to-end
})

# ── statistical moment matching (gated) ────────────────────────────────────

test_that("poisson draws have mean and variance ~ mu", {
  skip_on_cran()
  skip_if_not_slow()
  dt <- glm_data()
  nd <- dt[t == nrow(dt)]
  glm_p <- stats::glm(ycount ~ x, dt, family = stats::poisson())
  pp <- stats::predict(glm_p, newdata = nd, type = "link", se.fit = TRUE)
  mu <- as.numeric(stats::predict(glm_p, nd, type = "response"))
  set.seed(3)
  draws <- as.vector(getpi_glm(pp, glm_p$family, glm_p$df.residual, 1,
                               nsamples = 40000))
  expect_equal(mean(draws), mu, tolerance = 0.1)
  expect_equal(stats::var(draws), mu, tolerance = 0.25)   # Poisson var = mean
})

test_that("Gamma draws have mean ~ mu and variance ~ dispersion * mu^2", {
  skip_on_cran()
  skip_if_not_slow()
  dt <- glm_data()
  nd <- dt[t == nrow(dt)]
  glm_gam <- stats::glm(ypos ~ x, dt, family = stats::Gamma(link = "log"))
  disp <- summary(glm_gam)$dispersion
  gmp <- stats::predict(glm_gam, newdata = nd, type = "link", se.fit = TRUE)
  mu <- as.numeric(stats::predict(glm_gam, nd, type = "response"))
  set.seed(4)
  draws <- as.vector(getpi_glm(gmp, glm_gam$family, glm_gam$df.residual, disp,
                               nsamples = 40000))
  expect_equal(mean(draws), mu, tolerance = 0.1 * mu)
  expect_equal(stats::var(draws), disp * mu^2, tolerance = 0.3 * disp * mu^2)
})
