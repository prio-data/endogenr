# parametric_distribution (fitdistrplus) models ------------------------------
#
# Pins the t_ls auto-start rule, the .with_dist_visible() shim that makes the
# package's d/p/q/r functions resolvable by fitdist(), start-as-function
# plumbing, the unattached (endogenr::-only) fitting path, and the
# param_uncertainty engine opt-in.

# ── t_ls auto-start (no start argument) ─────────────────────────────────────

test_that("t_ls fits without start values and recovers location/scale", {
  skip_if_not_installed("fitdistrplus")

  set.seed(1)
  x <- 3 + 2 * stats::rt(4000, df = 6)
  m <- parametric_distribution_model(
    ~v, distribution = "t_ls",
    data = data.table::data.table(v = x, unit = 1L, time = seq_along(x)),
    ctx = panel_context("unit", "time", "sim"))

  expect_s3_class(m, "parametric_distribution")
  est <- m$fitted$estimate
  expect_true(all(c("df", "mu", "sigma") %in% names(est)))
  expect_lt(abs(est[["mu"]] - 3), 0.3)
  expect_lt(abs(est[["sigma"]] - 2), 0.4)
  expect_gt(est[["df"]], 2)          # heavy tails detected, finite df
  expect_true(all(is.finite(est)))
})

# ── .dist_start contract ────────────────────────────────────────────────────

test_that(".dist_start returns moment-matched t_ls starts and NULL otherwise", {
  set.seed(2)
  x <- 3 + 2 * stats::rt(4000, df = 6)

  st <- .dist_start("t_ls", x)
  expect_type(st, "list")
  expect_named(st, c("df", "mu", "sigma"))
  expect_gte(st$df, 2.5)
  expect_lte(st$df, 100)
  expect_equal(st$mu, stats::median(x))  # location start = robust median
  expect_gt(st$sigma, 0)

  # fitdist's own defaults cover the standard families -> no override.
  expect_null(.dist_start("norm", x))

  # Near-normal input (excess kurtosis <= 0.1) -> large-df fallback of 30.
  set.seed(3)
  z <- stats::rnorm(4000)
  expect_equal(.dist_start("t_ls", z)$df, 30)
})

# ── start as function(x), evaluated on the training window ─────────────────

test_that("start = function(x) flows through the pipeline and fits", {
  skip_if_not_installed("fitdistrplus")
  future::plan(future::sequential)
  on.exit(future::plan(future::sequential), add = TRUE)

  set.seed(4)
  dt <- data.table::as.data.table(expand.grid(unit = 1:3, time = 1:40))
  dt[, g := stats::rnorm(.N, mean = 0.02, sd = 0.05)]

  models <- list(
    build_model("parametric_distribution", formula = ~g, distribution = "t_ls",
                start = function(x) list(df = 10, mu = mean(x),
                                         sigma = stats::sd(x)))
  )
  setup <- setup_system(models, dt, train_start = 1, test_start = 35,
                        horizon = 6, groupvar = "unit", timevar = "time",
                        inner_sims = 5)
  set.seed(5)
  res <- simulate_system(fit_system(setup, nsim = 1))

  fc <- res[time >= 35]
  expect_gt(nrow(fc), 0)
  expect_true(all(is.finite(fc$g)))
})

# ── .with_dist_visible mechanics ────────────────────────────────────────────

test_that(".with_dist_visible exposes dt_ls during the call and cleans up", {
  inside <- .with_dist_visible(
    "t_ls", exists("dt_ls", envir = globalenv(), mode = "function"))
  expect_true(inside)
  # No binding is left behind in globalenv itself.
  expect_false(exists("dt_ls", envir = globalenv(), inherits = FALSE))
})

test_that(".with_dist_visible is a no-op for already-visible distributions", {
  expect_false(exists("dnorm", envir = globalenv(), inherits = FALSE))
  out <- .with_dist_visible("norm", stats::dnorm(0))
  expect_equal(out, stats::dnorm(0))
  # dnorm is visible via package:stats -> nothing was injected.
  expect_false(exists("dnorm", envir = globalenv(), inherits = FALSE))
})

test_that(".with_dist_visible cleans up when the expression errors", {
  expect_error(.with_dist_visible("t_ls", stop("boom")), "boom")
  expect_false(exists("dt_ls", envir = globalenv(), inherits = FALSE))
})

test_that(".with_dist_visible restores a shadowed non-function binding", {
  marker <- "user data, not a function"
  assign("dt_ls", marker, envir = globalenv())
  on.exit(if (exists("dt_ls", envir = globalenv(), inherits = FALSE))
            rm("dt_ls", envir = globalenv()),
          add = TRUE)

  inside <- .with_dist_visible(
    "t_ls", exists("dt_ls", envir = globalenv(), mode = "function"))
  expect_true(inside)
  # The user's non-function binding survives the round trip untouched.
  expect_identical(get("dt_ls", envir = globalenv(), inherits = FALSE), marker)
})

# ── unattached usage: endogenr:: without library(endogenr) ─────────────────

test_that("t_ls fits end-to-end in a subprocess that never attaches endogenr", {
  skip_on_cran()
  skip_if_not_slow()
  skip_if_not_installed("fitdistrplus")
  skip_if_not_installed("pkgload")
  # Needs the package source tree (fitdist resolves d/p/q/r by name from the
  # search path; this reproduces the reported endogenr::-only failure mode).
  skip_if_not(file.exists(testthat::test_path("..", "..", "DESCRIPTION")),
              "package source root not available")

  pkg_root <- normalizePath(testthat::test_path("..", ".."))
  script <- tempfile(fileext = ".R")
  writeLines(c(
    # Subprocess runs --vanilla: re-point it at this session's libraries.
    sprintf(".libPaths(unique(c(%s, .libPaths())))",
            paste(deparse(.libPaths()), collapse = "")),
    # Prefer the in-repo source over any (possibly stale) installed endogenr;
    # fall back to requireNamespace only if load_all is unusable.
    "loaded <- tryCatch({",
    sprintf("  pkgload::load_all(%s, attach = FALSE, export_all = FALSE, quiet = TRUE)",
            deparse(pkg_root)),
    "  TRUE",
    "}, error = function(e) FALSE)",
    "if (!loaded) stopifnot(requireNamespace('endogenr', quietly = TRUE))",
    "stopifnot(!'package:endogenr' %in% search())",
    "set.seed(42)",
    "d <- data.frame(unit = rep(1:2, each = 250), time = rep(1:250, 2),",
    "                v = stats::rt(500, 5) * 2 + 1)",
    "m <- endogenr::build_model('parametric_distribution', formula = ~v,",
    "                           distribution = 't_ls')",
    "s <- endogenr::setup_system(list(m), d, train_start = 1, test_start = 240,",
    "                            horizon = 5, groupvar = 'unit', timevar = 'time',",
    "                            inner_sims = 5)",
    "r <- endogenr::simulate_system(endogenr::fit_system(s, nsim = 1))",
    "stopifnot(all(is.finite(r$v[r$time >= 240])))"
  ), script)

  out <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
                                  c("--vanilla", script),
                                  stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  if (is.null(status)) status <- 0L
  expect_identical(status, 0L,
                   label = paste(c("subprocess exit status", out),
                                 collapse = "\n"))
})

# ── param_uncertainty engine opt-in ─────────────────────────────────────────

test_that("param_uncertainty = TRUE draws one parameter vector per simulation", {
  skip_if_not_installed("fitdistrplus")
  future::plan(future::sequential)
  on.exit(future::plan(future::sequential), add = TRUE)

  set.seed(6)
  dt <- data.table::as.data.table(expand.grid(unit = 1:4, time = 1:30))
  dt[, g := stats::rnorm(.N, mean = 1, sd = 0.5)]

  # The endogenr-level option lands on the model and is NOT forwarded to
  # fitdist() (an unused 'param_uncertainty' argument would abort the fit).
  m <- parametric_distribution_model(
    ~g, distribution = "norm", data = dt,
    ctx = panel_context("unit", "time", "sim"), param_uncertainty = TRUE)
  expect_true(m$param_uncertainty)
  expect_named(m$fitted$estimate, c("mean", "sd"))

  models <- list(
    build_model("parametric_distribution", formula = ~g, distribution = "norm",
                param_uncertainty = TRUE)
  )
  setup <- setup_system(models, dt, train_start = 1, test_start = 26,
                        horizon = 5, groupvar = "unit", timevar = "time",
                        inner_sims = 100)
  set.seed(7)
  res <- simulate_system(fit_system(setup, nsim = 1))

  fc <- res[time >= 26]
  expect_true(all(is.finite(fc$g)))
  per_sim <- fc[, .(m = mean(g)), by = ".sim"]
  expect_equal(nrow(per_sim), 100L)
  # Parameter draws differ across sims -> per-sim means genuinely vary.
  expect_gt(stats::sd(per_sim$m), 0)
})
