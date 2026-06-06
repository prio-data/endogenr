# Tests for R/experiments.R ------------------------------------------------

# ── Helpers ────────────────────────────────────────────────────────────────

make_exp_data <- function() {
  df <- expand.grid(gwcode = c(1, 2), year = 2000:2020)
  set.seed(1)
  df$x <- rnorm(nrow(df))
  df$y <- 0.5 * df$x + rnorm(nrow(df), sd = 0.2)
  data.table::as.data.table(df)
}

# A base system: linear y ~ lag(x), exogen x. The linear spec is the swap target.
base_system <- function() {
  list(
    build_model("linear", formula = y ~ lag(x), boot = "resid"),
    build_model("exogen", formula = ~x)
  )
}

# ── vary_model ───────────────────────────────────────────────────────────--

test_that("vary_model swaps the matched-outcome spec and preserves the rest", {
  base <- base_system()
  systems <- vary_model(base, list(
    e1 = build_model("linear", formula = y ~ lag(x), boot = "resid"),
    e2 = build_model("linear", formula = y ~ lag(x) + lag(y), boot = "resid")
  ))

  expect_named(systems, c("e1", "e2"))
  # Each variant is a full 2-spec system in the original order.
  expect_length(systems$e1, 2L)
  expect_identical(systems$e1[[2]]$type, "exogen")            # untouched
  expect_identical(systems$e2[[2]]$type, "exogen")            # untouched
  # The swapped slot carries the variant's formula.
  expect_identical(deparse(systems$e2[[1]]$formula),
                   deparse(y ~ lag(x) + lag(y)))
})

test_that("vary_model auto-names blank variants and rejects unmatched outcomes", {
  base <- base_system()

  # Unnamed entries become variant1, variant2, ...
  systems <- vary_model(base, list(
    build_model("linear", formula = y ~ lag(x), boot = "resid")
  ))
  expect_named(systems, "variant1")

  # An override whose outcome is absent from base is an error (swap, not append).
  expect_error(
    vary_model(base, list(bad = build_model("linear", formula = z ~ lag(x)))),
    "not produced"
  )
})

# ── run_experiments: model axis ─────────────────────────────────────────────

test_that("run_experiments stacks a model axis with readable labels + attrs", {
  dt <- make_exp_data()
  systems <- vary_model(base_system(), list(
    e1 = build_model("linear", formula = y ~ lag(x), boot = "resid"),
    e2 = build_model("linear", formula = y ~ lag(x) + lag(y), boot = "resid")
  ))

  set.seed(42)
  res <- run_experiments(
    data = dt, models = systems,
    train_start = 2000, test_start = 2010, horizon = 3,
    groupvar = "gwcode", timevar = "year",
    inner_sims = 2, nsim = 2, min_window = 5, parallel = "draws"
  )

  expect_s3_class(res, "data.table")
  expect_setequal(unique(res$.experiment), c("e1", "e2"))
  expect_true(all(c(".experiment", "model", "train_start", "test_start", ".sim")
                  %in% names(res)))
  # Leading metadata columns.
  expect_identical(names(res)[1:4],
                   c(".experiment", "model", "train_start", "test_start"))
  # Panel + experiment metadata stamped.
  expect_identical(attr(res, "panel_unit"), "gwcode")
  expect_identical(attr(res, "panel_time"), "year")
  meta <- attr(res, "experiments")
  expect_s3_class(meta, "data.table")
  expect_equal(nrow(meta), 2L)
  expect_setequal(meta$.experiment, c("e1", "e2"))
  expect_true(all(meta$test_start == 2010))
})

# ── run_experiments: varying test_start ─────────────────────────────────────

test_that("run_experiments varies test_start: labels + per-experiment spans", {
  dt <- make_exp_data()

  set.seed(42)
  res <- run_experiments(
    data = dt, models = base_system(),
    train_start = 2000, test_start = c(2010, 2015), horizon = 3,
    groupvar = "gwcode", timevar = "year",
    inner_sims = 2, nsim = 1, min_window = 5, parallel = "draws"
  )

  expect_setequal(unique(res$.experiment), c("test2010", "test2015"))

  # Each experiment spans [test_start, test_start + horizon - 1] (forecast-only).
  span <- res[, .(lo = min(year), hi = max(year)), by = .experiment]
  expect_equal(span[.experiment == "test2010", hi], 2012)  # 2010 + 3 - 1
  expect_equal(span[.experiment == "test2015", hi], 2017)  # 2015 + 3 - 1
  expect_equal(span[.experiment == "test2010", lo], 2010)
  expect_equal(span[.experiment == "test2015", lo], 2015)

  meta <- attr(res, "experiments")
  expect_setequal(meta$test_start, c(2010, 2015))
})

# ── run_experiments: Cartesian crossing ─────────────────────────────────────

test_that("run_experiments crosses model x test_start (Cartesian)", {
  dt <- make_exp_data()
  systems <- vary_model(base_system(), list(
    e1 = build_model("linear", formula = y ~ lag(x), boot = "resid"),
    e2 = build_model("linear", formula = y ~ lag(x) + lag(y), boot = "resid")
  ))

  set.seed(42)
  res <- run_experiments(
    data = dt, models = systems,
    train_start = 2000, test_start = c(2010, 2015), horizon = 2,
    groupvar = "gwcode", timevar = "year",
    inner_sims = 1, nsim = 1, min_window = 5, parallel = "draws"
  )

  meta <- attr(res, "experiments")
  expect_equal(nrow(meta), 4L)
  expect_equal(length(unique(meta$.experiment)), 4L)
  expect_setequal(unique(res$.experiment),
                  c("e1_test2010", "e1_test2015", "e2_test2010", "e2_test2015"))
})

# ── Reproducibility ─────────────────────────────────────────────────────────

test_that("run_experiments is reproducible within a scheme", {
  dt <- make_exp_data()
  systems <- vary_model(base_system(), list(
    e1 = build_model("linear", formula = y ~ lag(x), boot = "resid"),
    e2 = build_model("linear", formula = y ~ lag(x) + lag(y), boot = "resid")
  ))
  run <- function(scheme) {
    set.seed(7)
    run_experiments(
      data = dt, models = systems,
      train_start = 2000, test_start = 2010, horizon = 3,
      groupvar = "gwcode", timevar = "year",
      inner_sims = 2, nsim = 2, min_window = 5, parallel = scheme
    )
  }

  future::plan(future::sequential)
  a <- run("draws")
  b <- run("draws")
  expect_equal(a$y, b$y)
  expect_equal(a$.sim, b$.sim)

  c1 <- run("experiments")
  c2 <- run("experiments")
  expect_equal(c1$y, c2$y)
})

test_that("single-experiment 'draws' run equals the manual pipeline", {
  dt <- make_exp_data()
  future::plan(future::sequential)

  setup <- setup_system(
    models = base_system(), data = dt,
    train_start = 2000, test_start = 2010, horizon = 3,
    groupvar = "gwcode", timevar = "year", inner_sims = 2, min_window = 5
  )
  set.seed(99)
  manual <- simulate_system(fit_system(setup, nsim = 3))

  set.seed(99)
  res <- run_experiments(
    data = dt, models = base_system(),
    train_start = 2000, test_start = 2010, horizon = 3,
    groupvar = "gwcode", timevar = "year",
    inner_sims = 2, nsim = 3, min_window = 5, parallel = "draws"
  )

  # run_experiments only adds metadata columns; the simulated values match.
  expect_equal(res$y, manual$y)
  expect_equal(res$.sim, manual$.sim)
})

# ── keep_fits ───────────────────────────────────────────────────────────────

test_that("keep_fits attaches per-experiment fits usable by get_coefficients", {
  dt <- make_exp_data()
  systems <- vary_model(base_system(), list(
    e1 = build_model("linear", formula = y ~ lag(x), boot = "resid"),
    e2 = build_model("linear", formula = y ~ lag(x) + lag(y), boot = "resid")
  ))

  set.seed(42)
  res <- run_experiments(
    data = dt, models = systems,
    train_start = 2000, test_start = 2010, horizon = 2,
    groupvar = "gwcode", timevar = "year",
    inner_sims = 1, nsim = 2, min_window = 5, parallel = "draws",
    keep_fits = TRUE
  )

  fits <- attr(res, "fits")
  expect_named(fits, c("e1", "e2"))
  expect_s3_class(fits$e1, "endogenr_fitted_system")

  co <- get_coefficients(fits$e1)
  expect_s3_class(co, "data.table")
  expect_true(all(c(".draw", "outcome", "term", "estimate") %in% names(co)))
})

# ── get_experiment_accuracy ─────────────────────────────────────────────────

test_that("get_experiment_accuracy scores per experiment with its own test_start", {
  dt <- make_exp_data()

  set.seed(42)
  res <- run_experiments(
    data = dt, models = base_system(),
    train_start = 2000, test_start = c(2010, 2015), horizon = 3,
    groupvar = "gwcode", timevar = "year",
    inner_sims = 3, nsim = 2, min_window = 5, parallel = "draws"
  )

  acc <- get_experiment_accuracy(res, "y", dt)

  expect_s3_class(acc, "data.table")
  expect_true(all(c(".experiment", "gwcode", "horizon",
                    "crps", "mae", "winkler") %in% names(acc)))
  expect_setequal(unique(acc$.experiment), c("test2010", "test2015"))

  # Horizon is relative to each experiment's own test_start: both experiments
  # cover horizons 1..3, despite different forecast years.
  for (lab in c("test2010", "test2015")) {
    expect_setequal(acc[.experiment == lab, sort(unique(horizon))], 1:3)
  }
  # Scores are finite.
  expect_true(all(is.finite(acc$crps)))
})

test_that("get_experiment_accuracy honours a custom `by`", {
  dt <- make_exp_data()
  set.seed(42)
  res <- run_experiments(
    data = dt, models = base_system(),
    train_start = 2000, test_start = 2010, horizon = 3,
    groupvar = "gwcode", timevar = "year",
    inner_sims = 3, nsim = 1, min_window = 5, parallel = "draws"
  )

  acc <- get_experiment_accuracy(res, "y", dt, by = "gwcode")
  expect_true("gwcode" %in% names(acc))
  expect_false("horizon" %in% names(acc))
  # One row per (experiment, unit).
  expect_equal(nrow(acc), data.table::uniqueN(dt$gwcode))
})

# ── Parallel across experiments with user globals in a formula ──────────────

test_that("parallel = 'experiments' exports user globals for formula NSE", {
  skip_on_cran()
  dt <- make_exp_data()
  bump <- function(z) z + 0          # referenced via NSE inside a formula

  systems <- list(
    e1 = list(
      build_model("linear", formula = y ~ lag(bump(x)), boot = "resid"),
      build_model("exogen", formula = ~x)
    )
  )

  future::plan(future::multisession, workers = 2)
  on.exit(future::plan(future::sequential), add = TRUE)

  set.seed(42)
  res <- run_experiments(
    data = dt, models = systems,
    train_start = 2000, test_start = 2010, horizon = 2,
    groupvar = "gwcode", timevar = "year",
    inner_sims = 1, nsim = 2, min_window = 5,
    globals = list(bump = bump), parallel = "experiments"
  )

  expect_s3_class(res, "data.table")
  expect_setequal(unique(res$.experiment), "e1")
  expect_true(any(!is.na(res[year >= 2010, y])))
})

# ── window_config constructor ────────────────────────────────────────────────

test_that("window_config() creates a valid endogenr_window with correct fields", {
  wc <- window_config()
  expect_s3_class(wc, "endogenr_window")
  expect_equal(wc$window, "random")
  expect_null(wc$width)
  expect_equal(wc$step, 1L)
  expect_equal(wc$window_policy, "latest")
  expect_equal(wc$decay, 0.5)
  expect_null(wc$weights)

  # Rolling config captures all fields
  wc2 <- window_config("rolling", width = 10L, step = 2L,
                        window_policy = "decay", decay = 0.8)
  expect_equal(wc2$window, "rolling")
  expect_equal(wc2$width, 10L)
  expect_equal(wc2$step, 2L)
  expect_equal(wc2$window_policy, "decay")
  expect_equal(wc2$decay, 0.8)
})

test_that("window_config() validates arguments", {
  # Bad window enum
  expect_error(window_config("sliding"), "should be one of")
  # Bad window_policy enum
  expect_error(window_config(window_policy = "newest"), "should be one of")
  # Non-positive step
  expect_error(window_config(step = 0L), "positive integer")
  # Non-positive width
  expect_error(window_config("rolling", width = 0L), "positive integer")
  # decay out of range
  expect_error(window_config(decay = 1.5), "\\(0, 1\\]")
  expect_error(window_config(decay = 0), "\\(0, 1\\]")
  # bad weights type
  expect_error(window_config(weights = "bad"), "function")
})

test_that("window_config() warns when random + non-default policy", {
  expect_warning(window_config(window = "random", window_policy = "equal"),
                 "ignored")
  expect_warning(window_config(window = "random", weights = matrix(1, 1, 1)),
                 "ignored")
  # no warning for random + latest (the default)
  expect_no_warning(window_config(window = "random", window_policy = "latest"))
})

# ── run_experiments: windows axis ───────────────────────────────────────────

test_that("run_experiments windows axis: 2 configs produce labelled experiments with correct fit_mode", {
  dt <- make_exp_data()
  future::plan(future::sequential)

  set.seed(42)
  res <- run_experiments(
    data = dt, models = base_system(),
    train_start = 2000, test_start = 2010, horizon = 3,
    groupvar = "gwcode", timevar = "year",
    inner_sims = 2, nsim = 2, min_window = 5, parallel = "draws",
    windows = list(
      rand = window_config(),
      roll = window_config("rolling", width = 5, step = 2)
    ),
    keep_fits = TRUE
  )

  expect_s3_class(res, "data.table")
  # Both window labels present
  expect_setequal(unique(res$window_cfg), c("rand", "roll"))
  # Experiments named by window (single model, single test_start)
  expect_setequal(unique(res$.experiment), c("rand", "roll"))
  # window_cfg is the 5th leading column
  expect_equal(names(res)[5], "window_cfg")
  # experiments metadata carries window columns
  meta <- attr(res, "experiments")
  expect_true(all(c("window_cfg", "window", "window_policy") %in% names(meta)))
  expect_setequal(meta$window_cfg, c("rand", "roll"))
  expect_setequal(meta$window, c("random", "rolling"))

  # fits reflect the correct fit_mode
  fits <- attr(res, "fits")
  expect_equal(fits$rand$fit_mode, "random")
  expect_equal(fits$roll$fit_mode, "sliding")
  expect_equal(fits$roll$window, "rolling")
})

test_that("run_experiments crosses models x windows (Cartesian)", {
  dt <- make_exp_data()
  future::plan(future::sequential)
  systems <- vary_model(base_system(), list(
    e1 = build_model("linear", formula = y ~ lag(x), boot = "resid"),
    e2 = build_model("linear", formula = y ~ lag(x) + lag(y), boot = "resid")
  ))

  set.seed(7)
  res <- run_experiments(
    data = dt, models = systems,
    train_start = 2000, test_start = 2010, horizon = 2,
    groupvar = "gwcode", timevar = "year",
    inner_sims = 1, nsim = 1, min_window = 5, parallel = "draws",
    windows = list(
      rand = window_config(),
      roll = window_config("rolling", width = 5)
    )
  )

  meta <- attr(res, "experiments")
  expect_equal(nrow(meta), 4L)  # 2 models x 2 windows
  expect_setequal(unique(res$.experiment),
                  c("e1_rand", "e1_roll", "e2_rand", "e2_roll"))
})

test_that("single window config applies uniformly without adding a label component", {
  dt <- make_exp_data()
  future::plan(future::sequential)

  set.seed(42)
  res <- run_experiments(
    data = dt, models = base_system(),
    train_start = 2000, test_start = 2010, horizon = 2,
    groupvar = "gwcode", timevar = "year",
    inner_sims = 2, nsim = 2, min_window = 5, parallel = "draws",
    windows = window_config("expanding"),
    keep_fits = TRUE
  )

  # Only 1 window config → label doesn't include window component
  expect_equal(unique(res$.experiment), "model1")
  # fit_mode is sliding (expanding)
  fits <- attr(res, "fits")
  expect_equal(fits$model1$fit_mode, "sliding")
  expect_equal(fits$model1$window, "expanding")
})

test_that("default windows = NULL reproduces existing manual-pipeline behaviour", {
  # This mirrors 'single-experiment draws run equals manual pipeline' but now
  # with windows = NULL (default), confirming backward compatibility.
  dt <- make_exp_data()
  future::plan(future::sequential)

  setup <- setup_system(
    models = base_system(), data = dt,
    train_start = 2000, test_start = 2010, horizon = 3,
    groupvar = "gwcode", timevar = "year", inner_sims = 2, min_window = 5
  )
  set.seed(99)
  manual <- simulate_system(fit_system(setup, nsim = 3))

  set.seed(99)
  res <- run_experiments(
    data = dt, models = base_system(),
    train_start = 2000, test_start = 2010, horizon = 3,
    groupvar = "gwcode", timevar = "year",
    inner_sims = 2, nsim = 3, min_window = 5, parallel = "draws"
    # windows = NULL by default
  )
  expect_equal(res$y, manual$y)
  expect_equal(res$.sim, manual$.sim)
})

test_that("get_experiment_accuracy carries window_cfg column through", {
  dt <- make_exp_data()
  future::plan(future::sequential)

  set.seed(42)
  res <- run_experiments(
    data = dt, models = base_system(),
    train_start = 2000, test_start = 2010, horizon = 3,
    groupvar = "gwcode", timevar = "year",
    inner_sims = 3, nsim = 2, min_window = 5, parallel = "draws",
    windows = list(
      rand = window_config(),
      roll = window_config("rolling", width = 5)
    )
  )

  acc <- get_experiment_accuracy(res, "y", dt)
  expect_true("window_cfg" %in% names(acc))
  expect_setequal(unique(acc$window_cfg), c("rand", "roll"))
  expect_true(all(is.finite(acc$crps)))
})
