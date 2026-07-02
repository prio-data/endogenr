# Test gating helpers ---------------------------------------------------------
#
# Deterministic unit tests run always-on. Replication-heavy statistical tests
# (coverage, calibration, parameter recovery, true-vs-wrong model ordering) are
# gated so CI stays fast; enable them locally/nightly with
#   ENDOGENR_SLOW_TESTS=true
# and run via a sequential future plan for reproducibility.

# TRUE when the slow / statistical battery is opted in via the env var.
slow_tests_enabled <- function() {
  val <- tolower(trimws(Sys.getenv("ENDOGENR_SLOW_TESTS", "")))
  val %in% c("true", "1", "yes", "on")
}

# Skip the calling test unless the slow battery is enabled.
skip_if_not_slow <- function() {
  testthat::skip_if_not(slow_tests_enabled(),
                        "ENDOGENR_SLOW_TESTS not enabled (slow/statistical test)")
}

# The fable family needs several Suggests packages; skip the whole stack at once.
skip_if_no_fable <- function() {
  testthat::skip_if_not_installed("fable")
  testthat::skip_if_not_installed("fabletools")
  testthat::skip_if_not_installed("tsibble")
  testthat::skip_if_not_installed("distributional")
}

# glmmTMB needs only the glmmTMB package itself.
skip_if_no_glmmTMB <- function() {
  testthat::skip_if_not_installed("glmmTMB")
}

# gamlss needs both gamlss and gamlss.dist.
skip_if_no_gamlss <- function() {
  testthat::skip_if_not_installed("gamlss")
  testthat::skip_if_not_installed("gamlss.dist")
}
