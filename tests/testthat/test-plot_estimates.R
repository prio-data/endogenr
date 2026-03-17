# Helpers ------------------------------------------------------------------

#' Build a minimal fitted-model list element matching the structure
#' that plot_estimates() expects: $coefs (broom::tidy-like), $gof
#' (broom::glance-like), and $outcome.
make_mock_model <- function(outcome, terms, n = 100) {
  coefs <- data.frame(
    term      = c("(Intercept)", terms),
    estimate  = stats::rnorm(length(terms) + 1),
    std.error = abs(stats::rnorm(length(terms) + 1, sd = 0.1)),
    statistic = stats::rnorm(length(terms) + 1),
    p.value   = stats::runif(length(terms) + 1),
    stringsAsFactors = FALSE
  )
  gof <- data.frame(
    r.squared = stats::runif(1),
    nobs      = n
  )
  list(outcome = outcome, coefs = coefs, gof = gof)
}

#' Wrap a set of fitted-model lists into the named-list-by-test-start
#' structure that plot_estimates() consumes.
make_mock_models <- function(test_starts, model_specs) {
  stats::setNames(
    lapply(test_starts, function(ts) {
      list(fitted_models = lapply(model_specs, function(spec) {
        make_mock_model(spec$outcome, spec$terms)
      }))
    }),
    as.character(test_starts)
  )
}

# Tests --------------------------------------------------------------------

test_that("plot_estimates returns a ggplot when show_gof = FALSE", {
  models <- make_mock_models(
    test_starts = c(2010, 2012),
    model_specs = list(
      list(outcome = "gdp",      terms = c("lag_gdp", "pop")),
      list(outcome = "conflict", terms = c("lag_conflict", "gdp"))
    )
  )

  p <- plot_estimates(models, show_gof = FALSE)
  expect_s3_class(p, "ggplot")
})

test_that("plot_estimates returns a patchwork when show_gof = TRUE", {
  models <- make_mock_models(
    test_starts = c(2010, 2012),
    model_specs = list(
      list(outcome = "gdp",      terms = c("lag_gdp", "pop")),
      list(outcome = "conflict", terms = c("lag_conflict", "gdp"))
    )
  )

  p <- plot_estimates(models, show_gof = TRUE)
  expect_s3_class(p, "patchwork")
})

test_that("plot_estimates handles models with varying number of terms", {
  models <- make_mock_models(
    test_starts = c(2010, 2012),
    model_specs = list(
      list(outcome = "gdp",      terms = c("lag_gdp", "pop", "trade")),
      list(outcome = "conflict", terms = c("lag_conflict"))
    )
  )

  p <- plot_estimates(models, show_gof = FALSE)
  expect_s3_class(p, "ggplot")

  # The facet should use ncol = 3 (max terms) with empty cells for conflict
  facet <- p$facet
  expect_equal(facet$params$ncol, 3L)
})

test_that("outcome_labels are applied correctly", {
  models <- make_mock_models(
    test_starts = c(2010),
    model_specs = list(
      list(outcome = "gdp",      terms = c("lag_gdp")),
      list(outcome = "conflict", terms = c("lag_conflict"))
    )
  )

  p <- plot_estimates(
    models,
    outcome_labels = c(gdp = "GDP per capita", conflict = "Armed conflict"),
    show_gof = FALSE
  )

  # Check that outcome levels use the labels
  pdata <- ggplot2::ggplot_build(p)$layout$layout
  expect_true("GDP per capita" %in% as.character(pdata$outcome))
  expect_true("Armed conflict" %in% as.character(pdata$outcome))
})
