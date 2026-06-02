# Tests for R/systemsim.R --------------------------------------------------

# ── Helper: minimal test data ──────────────────────────────────────────────

make_test_data <- function() {
  df <- expand.grid(gwcode = c(1, 2), year = 2000:2015)
  set.seed(1)
  df$x <- rnorm(nrow(df))
  df$y <- 0.5 * df$x + rnorm(nrow(df), sd = 0.2)
  data.table::as.data.table(df)
}

# Fit + simulate in one go (the 3-step pipeline, minus setup).
make_results <- function(setup, nsim = 1L, seed = 42L) {
  set.seed(seed)
  fit <- fit_system(setup, nsim = nsim)
  simulate_system(fit)
}

# ── prepare_simulation_data ──────────────────────────────────────────────

test_that("prepare_simulation_data returns a data.table with correct dimensions", {
  dt <- make_test_data()
  ctx <- panel_context(unit = "gwcode", time = "year")
  ctx$sim <- "sim"

  sim_data <- prepare_simulation_data(
    data       = dt,
    ctx        = ctx,
    train_start = 2000,
    test_start  = 2010,
    horizon     = 3,
    inner_sims  = 2
  )

  expect_s3_class(sim_data, "data.table")

  # Should span 2000–2012 (train_start to test_start + horizon - 1)
  expect_equal(range(sim_data$year), c(2000, 2012))

  # 2 groups × 13 years × 2 inner_sims = 52 rows
  expect_equal(nrow(sim_data), 2 * 13 * 2)

  # Training rows should have non-NA y; forecast rows should be NA
  train_rows <- sim_data[sim_data$year < 2010]
  expect_false(all(is.na(train_rows$y)))

  forecast_rows <- sim_data[sim_data$year >= 2010]
  expect_true(all(is.na(forecast_rows$y)))
})

test_that("prepare_simulation_data includes sim column", {
  dt <- make_test_data()
  ctx <- panel_context(unit = "gwcode", time = "year")
  ctx$sim <- "sim"

  sim_data <- prepare_simulation_data(
    data       = dt,
    ctx        = ctx,
    train_start = 2000,
    test_start  = 2010,
    horizon     = 1,
    inner_sims  = 5
  )

  expect_true("sim" %in% names(sim_data))
  expect_equal(sort(unique(sim_data$sim)), 1:5)
})

# ── process_independent_models ───────────────────────────────────────────

test_that("process_independent_models validates test_start", {
  dt <- make_test_data()
  ctx <- panel_context(unit = "gwcode", time = "year")
  ctx$sim <- "sim"

  sim_data <- prepare_simulation_data(
    data = dt, ctx = ctx,
    train_start = 2000, test_start = 2010, horizon = 1, inner_sims = 1
  )

  # Should error because test_start is character
  expect_error(
    process_independent_models(sim_data, list(NULL), ctx, "2010", 1, 1),
    "numeric"
  )
})

# ── setup_system (integration) ────────────────────────────────────────────

test_that("setup_system returns expected structure and does not fit", {
  dt <- make_test_data()

  model_system <- list(
    build_model("linear", formula = y ~ lag(x)),
    build_model("exogen", formula = ~x)
  )

  result <- setup_system(
    models      = model_system,
    data        = dt,
    train_start = 2000,
    test_start  = 2010,
    horizon     = 3,
    groupvar    = "gwcode",
    timevar     = "year",
    inner_sims  = 2
  )

  expect_s3_class(result, "endogenr_system_setup")
  expect_true(all(c(
    "simulation_data", "specs", "test_start", "horizon", "execution_order",
    "groupvar", "timevar", "inner_sims", "ctx", "fit_ctx", "min_window"
  ) %in% names(result)))

  # No fitting happens in setup: there are no stored model objects.
  expect_false("fitted_models" %in% names(result))
  expect_false("fitted_draws" %in% names(result))

  expect_equal(result$test_start, 2010)
  expect_equal(result$horizon, 3)
  expect_equal(result$inner_sims, 2)
  expect_s3_class(result$simulation_data, "data.table")
  expect_s3_class(result$ctx, "panel_context")
})

test_that("fit_system stores coefs and gof for linear models", {
  dt <- make_test_data()

  model_system <- list(
    build_model("linear", formula = y ~ lag(x)),
    build_model("exogen", formula = ~x)
  )

  setup <- setup_system(
    models      = model_system,
    data        = dt,
    train_start = 2000,
    test_start  = 2010,
    horizon     = 1,
    groupvar    = "gwcode",
    timevar     = "year",
    inner_sims  = 1
  )

  fit <- fit_system(setup, nsim = 1)
  expect_s3_class(fit, "endogenr_fitted_system")

  fm <- fit$fitted_models[[1]]
  expect_true(!is.null(fm$coefs))
  expect_true(!is.null(fm$gof))
  expect_true("estimate" %in% names(fm$coefs))
  expect_true("r.squared" %in% names(fm$gof))
})

test_that("setup_system builds execution_order from specs (no fitting)", {
  dt <- make_test_data()

  model_system <- list(
    build_model("linear", formula = y ~ lag(x)),
    build_model("exogen", formula = ~x)
  )

  result <- setup_system(
    models      = model_system,
    data        = dt,
    train_start = 2000,
    test_start  = 2010,
    horizon     = 1,
    groupvar    = "gwcode",
    timevar     = "year",
    inner_sims  = 1
  )

  expect_true("y" %in% result$execution_order)
})

test_that("setup_system works with deterministic model", {
  dt <- make_test_data()

  # Self-lag deterministic: a valid model with no unmodeled same-period predictor
  # (an unlagged data column such as `x` would now be rejected by closure checks).
  model_system <- list(
    build_model("deterministic", formula = y ~ I(lag(y) * 1.02))
  )

  result <- setup_system(
    models      = model_system,
    data        = dt,
    train_start = 2000,
    test_start  = 2010,
    horizon     = 1,
    groupvar    = "gwcode",
    timevar     = "year",
    inner_sims  = 1
  )

  expect_s3_class(result, "endogenr_system_setup")
  expect_true("y" %in% result$execution_order)
})

test_that("setup_system handles multiple models with dependency ordering", {
  df <- expand.grid(gwcode = c(1, 2), year = 2000:2015)
  set.seed(1)
  df$x <- rnorm(nrow(df))
  df$y <- 0.5 * df$x + rnorm(nrow(df), sd = 0.2)
  df$z <- 0.3 * df$y + rnorm(nrow(df), sd = 0.1)
  dt <- data.table::as.data.table(df)

  model_system <- list(
    build_model("linear", formula = y ~ lag(x)),
    build_model("linear", formula = z ~ lag(y)),
    build_model("exogen", formula = ~x)
  )

  result <- setup_system(
    models      = model_system,
    data        = dt,
    train_start = 2000,
    test_start  = 2010,
    horizon     = 1,
    groupvar    = "gwcode",
    timevar     = "year",
    inner_sims  = 1
  )

  # y must come before z in execution order
  ord <- result$execution_order
  expect_true(which(ord == "y") < which(ord == "z"))
})

test_that("setup_system accepts tsibble input", {
  skip_if_not_installed("tsibble")
  df <- expand.grid(gwcode = c(1, 2), year = 2000:2015)
  set.seed(1)
  df$x <- rnorm(nrow(df))
  df$y <- 0.5 * df$x + rnorm(nrow(df), sd = 0.2)
  ts <- tsibble::tsibble(df, key = "gwcode", index = "year")

  model_system <- list(
    build_model("linear", formula = y ~ lag(x)),
    build_model("exogen", formula = ~x)
  )

  result <- setup_system(
    models      = model_system,
    data        = ts,
    train_start = 2000,
    test_start  = 2010,
    horizon     = 1,
    groupvar    = "gwcode",
    timevar     = "year",
    inner_sims  = 1
  )

  # Should still work — tsibble is converted to data.table internally
  expect_s3_class(result$simulation_data, "data.table")
})

# ── fit_system ─────────────────────────────────────────────────────────────

test_that("fit_system produces nsim refit draws with differing estimates", {
  df <- expand.grid(gwcode = c(1, 2), year = 2000:2030)
  set.seed(7)
  df$x <- rnorm(nrow(df))
  df$y <- 0.5 * df$x + rnorm(nrow(df), sd = 0.3)
  dt <- data.table::as.data.table(df)

  setup <- setup_system(
    models      = list(build_model("linear", formula = y ~ lag(x), boot = "resid"),
                       build_model("exogen", formula = ~x)),
    data        = dt,
    train_start = 2000,
    test_start  = 2025,
    horizon     = 2,
    groupvar    = "gwcode",
    timevar     = "year",
    inner_sims  = 1,
    min_window  = 10
  )

  set.seed(123)
  fit <- fit_system(setup, nsim = 4)

  expect_length(fit$fitted_draws, 4L)
  expect_equal(fit$nsim, 4L)

  # The refit linear model's slope must vary from draw to draw.
  slopes <- vapply(fit$fitted_draws, function(models) {
    unname(stats::coef(models[[1]]$fitted)[["lag_x"]])
  }, numeric(1))
  expect_gt(length(unique(slopes)), 1L)
})

test_that("fit_system shares a single fit across draws when min_window is NULL", {
  dt <- make_test_data()

  setup <- setup_system(
    models      = list(build_model("linear", formula = y ~ lag(x)),
                       build_model("exogen", formula = ~x)),
    data        = dt,
    train_start = 2000,
    test_start  = 2010,
    horizon     = 1,
    groupvar    = "gwcode",
    timevar     = "year",
    inner_sims  = 1
  )

  fit <- fit_system(setup, nsim = 3)
  expect_length(fit$fitted_draws, 3L)

  slopes <- vapply(fit$fitted_draws, function(models) {
    unname(stats::coef(models[[1]]$fitted)[["lag_x"]])
  }, numeric(1))
  expect_equal(length(unique(slopes)), 1L)

  # No random window: get_coefficients reports the full [train_start,
  # test_start - 1] range for every draw.
  gc <- get_coefficients(fit)
  expect_true(all(gc$.window_start == 2000))
  expect_true(all(gc$.window_end == 2009))
})

test_that("fit_system strips the cached training frame but predict still works", {
  dt <- make_test_data()
  set.seed(2)
  dt$z <- 0.4 * dt$x + rnorm(nrow(dt), sd = 0.2)

  setup <- setup_system(
    models      = list(build_model("linear", formula = y ~ lag(x)),
                       build_model("glm", formula = z ~ lag(x)),
                       build_model("exogen", formula = ~x)),
    data        = dt,
    train_start = 2000,
    test_start  = 2010,
    horizon     = 1,
    groupvar    = "gwcode",
    timevar     = "year",
    inner_sims  = 1
  )

  fit <- fit_system(setup, nsim = 1)

  lin <- fit$fitted_models[[1]]
  gl  <- fit$fitted_models[[2]]
  expect_null(lin$data)
  expect_null(gl$data)

  # predict()/model.matrix() must still work on the stripped objects.
  pred_names <- setdiff(names(stats::coef(lin$fitted)), "(Intercept)")
  nd <- as.data.frame(stats::setNames(
    rep(list(c(0.1, -0.2)), length(pred_names)), pred_names
  ))
  pp_lin <- predict(lin$fitted, newdata = nd, se.fit = TRUE)
  expect_length(pp_lin$fit, 2L)

  pp_glm <- predict(gl$fitted, newdata = nd, type = "link", se.fit = TRUE)
  expect_length(pp_glm$fit, 2L)
})

# ── get_coefficients / plot_coefficients ──────────────────────────────────

test_that("get_coefficients returns a tidy across-draw table", {
  df <- expand.grid(gwcode = c(1, 2), year = 2000:2030)
  set.seed(7)
  df$x <- rnorm(nrow(df))
  df$y <- 0.5 * df$x + rnorm(nrow(df), sd = 0.3)
  dt <- data.table::as.data.table(df)

  setup <- setup_system(
    models      = list(build_model("linear", formula = y ~ lag(x), boot = "resid"),
                       build_model("exogen", formula = ~x)),
    data        = dt,
    train_start = 2000,
    test_start  = 2025,
    horizon     = 2,
    groupvar    = "gwcode",
    timevar     = "year",
    inner_sims  = 1,
    min_window  = 10
  )

  set.seed(123)
  fit <- fit_system(setup, nsim = 5)
  gc <- get_coefficients(fit)

  expect_s3_class(gc, "data.table")
  expect_true(all(c(".draw", ".window_start", ".window_end", "outcome", "term",
                    "estimate", "std.error", "statistic", "p.value") %in% names(gc)))

  # One coefficient-bearing model (the linear), tidied to 2 terms per draw.
  n_terms <- nrow(fit$fitted_models[[1]]$coefs)
  expect_equal(nrow(gc), 5L * n_terms)
  expect_setequal(unique(gc$.draw), 1:5)
  expect_setequal(unique(gc$outcome), "y")

  # The refit window is recorded per draw, lies inside the trainable range,
  # and varies across draws (random windows).
  expect_true(all(gc$.window_start >= 2000 & gc$.window_end <= 2025))
  win_by_draw <- unique(gc[, .(.draw, .window_start, .window_end)])
  expect_gt(nrow(unique(win_by_draw[, .(.window_start, .window_end)])), 1L)

  # Across-draw spread for the refit slope.
  slope <- gc[gc$term != "(Intercept)"]
  expect_gt(stats::sd(slope$estimate), 0)
})

test_that("plot_coefficients returns a ggplot", {
  df <- expand.grid(gwcode = c(1, 2), year = 2000:2030)
  set.seed(7)
  df$x <- rnorm(nrow(df))
  df$y <- 0.5 * df$x + rnorm(nrow(df), sd = 0.3)
  dt <- data.table::as.data.table(df)

  setup <- setup_system(
    models      = list(build_model("linear", formula = y ~ lag(x), boot = "resid"),
                       build_model("exogen", formula = ~x)),
    data        = dt,
    train_start = 2000,
    test_start  = 2025,
    horizon     = 2,
    groupvar    = "gwcode",
    timevar     = "year",
    inner_sims  = 1,
    min_window  = 10
  )

  set.seed(123)
  fit <- fit_system(setup, nsim = 4)
  p <- plot_coefficients(fit)
  expect_s3_class(p, "ggplot")
  expect_silent(ggplot2::ggplot_build(p))
})

# ── simulate_system (minimal integration) ─────────────────────────────────

test_that("simulate_system runs and returns expected output", {
  dt <- make_test_data()

  model_system <- list(
    build_model("linear", formula = y ~ lag(x)),
    build_model("exogen", formula = ~x)
  )

  setup <- setup_system(
    models      = model_system,
    data        = dt,
    train_start = 2000,
    test_start  = 2013,
    horizon     = 2,
    groupvar    = "gwcode",
    timevar     = "year",
    inner_sims  = 2
  )

  res <- make_results(setup, nsim = 1)

  expect_s3_class(res, "data.table")
  expect_true(".sim" %in% names(res))
  expect_true("y" %in% names(res))
  expect_true("gwcode" %in% names(res))
  expect_true("year" %in% names(res))

  # Should have rows for both training and forecast periods
  expect_true(all(c(2013, 2014) %in% res$year))
})

test_that("simulate_system derives .sim count from nsim * inner_sims", {
  dt <- make_test_data()

  setup <- setup_system(
    models      = list(build_model("linear", formula = y ~ lag(x)),
                       build_model("exogen", formula = ~x)),
    data        = dt,
    train_start = 2000,
    test_start  = 2013,
    horizon     = 2,
    groupvar    = "gwcode",
    timevar     = "year",
    inner_sims  = 2
  )

  res <- make_results(setup, nsim = 3)
  expect_equal(length(unique(res$.sim)), 3L * 2L)
})

# ── sim_to_dist ──────────────────────────────────────────────────────────

test_that("sim_to_dist nests simulations into list-columns", {
  dt <- make_test_data()

  setup <- setup_system(
    models      = list(build_model("linear", formula = y ~ lag(x)), build_model("exogen", formula = ~x)),
    data        = dt,
    train_start = 2000,
    test_start  = 2013,
    horizon     = 2,
    groupvar    = "gwcode",
    timevar     = "year",
    inner_sims  = 2
  )

  res <- make_results(setup, nsim = 1)
  forecast_res <- res[res$year >= 2013]

  # Remove .sim from ctx for sim_to_dist (it groups by unit+time)
  ctx <- setup$ctx

  dist_result <- sim_to_dist(forecast_res, "y", ctx)

  expect_s3_class(dist_result, "data.table")
  # y column should be a list-column of numeric vectors
  expect_true(is.list(dist_result$y))
  expect_true(is.numeric(dist_result$y[[1]]))
})

# ── get_accuracy ─────────────────────────────────────────────────────────

test_that("get_accuracy computes CRPS, MAE, and Winkler scores", {
  dt <- make_test_data()

  setup <- setup_system(
    models      = list(build_model("linear", formula = y ~ lag(x)), build_model("exogen", formula = ~x)),
    data        = dt,
    train_start = 2000,
    test_start  = 2013,
    horizon     = 2,
    groupvar    = "gwcode",
    timevar     = "year",
    inner_sims  = 2
  )

  res <- make_results(setup, nsim = 2)
  forecast_res <- res[res$year >= 2013]

  ctx <- setup$ctx
  truth <- dt[dt$year >= 2013 & dt$year <= 2014]

  acc <- get_accuracy(forecast_res, "y", truth, ctx)

  expect_s3_class(acc, "data.table")
  expect_true(all(c("crps", "mae", "winkler") %in% names(acc)))
  expect_true(all(acc$crps >= 0, na.rm = TRUE))
  expect_true(all(acc$mae >= 0, na.rm = TRUE))
})
# ── ctx inference from panel attributes ──────────────────────────────────

test_that("simulate_system stamps panel_unit and panel_time attributes", {
  dt <- make_test_data()

  setup <- setup_system(
    models      = list(build_model("linear", formula = y ~ lag(x)), build_model("exogen", formula = ~x)),
    data        = dt,
    train_start = 2000,
    test_start  = 2013,
    horizon     = 2,
    groupvar    = "gwcode",
    timevar     = "year",
    inner_sims  = 1
  )

  res <- make_results(setup, nsim = 1)

  expect_identical(attr(res, "panel_unit"), "gwcode")
  expect_identical(attr(res, "panel_time"), "year")
})

test_that("get_accuracy infers ctx from simulation_results attributes", {
  dt <- make_test_data()

  setup <- setup_system(
    models      = list(build_model("linear", formula = y ~ lag(x)), build_model("exogen", formula = ~x)),
    data        = dt,
    train_start = 2000,
    test_start  = 2013,
    horizon     = 2,
    groupvar    = "gwcode",
    timevar     = "year",
    inner_sims  = 2
  )

  res <- make_results(setup, nsim = 2)
  forecast_res <- res[res$year >= 2013]
  truth <- dt[dt$year >= 2013 & dt$year <= 2014]

  # No ctx supplied — should infer from res attributes
  acc <- get_accuracy(forecast_res, "y", truth)

  expect_s3_class(acc, "data.table")
  expect_true(all(c("crps", "mae", "winkler") %in% names(acc)))
})

test_that("get_accuracy infers ctx from truth attributes when results lack them", {
  dt <- make_test_data()

  setup <- setup_system(
    models      = list(build_model("linear", formula = y ~ lag(x)), build_model("exogen", formula = ~x)),
    data        = dt,
    train_start = 2000,
    test_start  = 2013,
    horizon     = 2,
    groupvar    = "gwcode",
    timevar     = "year",
    inner_sims  = 2
  )

  res <- make_results(setup, nsim = 2)
  forecast_res <- res[res$year >= 2013]

  # Strip attributes from simulation results
  data.table::setattr(forecast_res, "panel_unit", NULL)
  data.table::setattr(forecast_res, "panel_time", NULL)

  # Stamp attributes on truth (mimicking paneltools::as_panel)
  truth <- dt[dt$year >= 2013 & dt$year <= 2014]
  data.table::setattr(truth, "panel_unit", "gwcode")
  data.table::setattr(truth, "panel_time", "year")

  acc <- get_accuracy(forecast_res, "y", truth)

  expect_s3_class(acc, "data.table")
  expect_true(all(c("crps", "mae", "winkler") %in% names(acc)))
})

test_that("plotsim infers ctx from attributes", {
  dt <- make_test_data()

  setup <- setup_system(
    models      = list(build_model("linear", formula = y ~ lag(x)), build_model("exogen", formula = ~x)),
    data        = dt,
    train_start = 2000,
    test_start  = 2013,
    horizon     = 2,
    groupvar    = "gwcode",
    timevar     = "year",
    inner_sims  = 2
  )

  res <- make_results(setup, nsim = 1)

  p <- plotsim(res, "y", units = c(1, 2), true_data = dt)
  expect_s3_class(p, "ggplot")
})

test_that("sim_to_dist infers ctx from simulation_results attributes", {
  dt <- make_test_data()

  setup <- setup_system(
    models      = list(build_model("linear", formula = y ~ lag(x)), build_model("exogen", formula = ~x)),
    data        = dt,
    train_start = 2000,
    test_start  = 2013,
    horizon     = 2,
    groupvar    = "gwcode",
    timevar     = "year",
    inner_sims  = 2
  )

  res <- make_results(setup, nsim = 1)
  forecast_res <- res[res$year >= 2013]

  dist_result <- sim_to_dist(forecast_res, "y")

  expect_s3_class(dist_result, "data.table")
  expect_true(is.list(dist_result$y))
})

test_that("get_accuracy errors helpfully when ctx cannot be inferred", {
  # Plain data.frames with no panel attributes
  sim <- data.table::data.table(gwcode = 1, year = 2013, y = 1.0, .sim = 1L)
  truth <- data.table::data.table(gwcode = 1, year = 2013, y = 1.0)

  expect_error(
    get_accuracy(sim, "y", truth),
    "Could not infer panel context"
  )
})
