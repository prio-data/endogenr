# Tests for R/panel_transform.R -------------------------------------------

# Rewrite an expression string; return the rewritten deparse and the ts map.
.rw <- function(expr_str, ts_fns = NULL) {
  st <- new.env(parent = emptyenv())
  st$ts_fns  <- union(.pt_ts_fns, ts_fns)
  st$map     <- list()
  st$keys    <- list()
  st$counter <- 0L
  out <- .rewrite_panel_formula(str2lang(expr_str), st)
  list(expr = deparse(out),
       map  = vapply(st$map, deparse, character(1)))
}

# ── .rewrite_panel_formula ─────────────────────────────────────────────────

test_that("design operators are kept and time-series terms are extracted", {
  expect_equal(.rw("factor(region)")$expr, "factor(region)")
  expect_equal(.rw("lag(x)")$expr, ".pt1")
  expect_equal(.rw("poly(lag(x), 2)")$expr, "poly(.pt1, 2)")
  expect_equal(.rw("factor(region):lag(x)")$expr, "factor(region):.pt1")
  expect_equal(.rw("log(lag(x))")$expr, "log(.pt1)")
  expect_equal(.rw("-1 + factor(region) + log(gdppc)")$expr,
               "-1 + factor(region) + log(gdppc)")
})

test_that("identical time-series sub-expressions are deduplicated", {
  r <- .rw("lag(x) + lag(x)")
  expect_equal(r$expr, ".pt1 + .pt1")
  expect_length(r$map, 1L)
})

test_that("pkg::-qualified time-series functions are matched", {
  r <- .rw("data.table::shift(x, 2)")
  expect_equal(r$expr, ".pt1")
  expect_equal(unname(r$map[[".pt1"]]), "data.table::shift(x, 2)")
})

test_that("ts_fns registers extra within-unit functions", {
  expect_equal(.rw("myroll(x, 3)")$expr, "myroll(x, 3)")            # not registered
  expect_equal(.rw("myroll(x, 3)", ts_fns = "myroll")$expr, ".pt1") # registered
})

# ── .apply_ts_map ──────────────────────────────────────────────────────────

test_that(".apply_ts_map lags within unit, time-ordered, with no cross-unit bleed", {
  dt <- data.table::data.table(
    g = c(1, 1, 1, 2, 2, 2),
    t = c(2, 1, 3, 3, 1, 2),          # deliberately shuffled within unit
    x = c(20, 10, 30, 300, 100, 200)
  )
  map <- list(.pt1 = quote(lag(log(x))))
  out <- .apply_ts_map(map, dt, "g", "t", environment())

  # earliest row of each unit is NA; values shift within unit, not across units
  expect_true(is.na(out[g == 1 & t == 1, .pt1]))
  expect_true(is.na(out[g == 2 & t == 1, .pt1]))
  expect_equal(out[g == 1 & t == 2, .pt1], log(10))
  expect_equal(out[g == 1 & t == 3, .pt1], log(20))
  expect_equal(out[g == 2 & t == 2, .pt1], log(100))
  expect_equal(out[g == 2 & t == 3, .pt1], log(200))

  # input is not mutated
  expect_false(".pt1" %in% names(dt))
})

test_that(".apply_ts_map with an empty map returns a keyed copy", {
  dt <- data.table::data.table(g = c(2, 1), t = c(1, 1), x = c(1, 2))
  out <- .apply_ts_map(list(), dt, "g", "t", environment())
  expect_equal(data.table::key(out), c("g", "t"))
  expect_equal(names(out), c("g", "t", "x"))
})

# ── panel_materialize ──────────────────────────────────────────────────────

test_that("panel_materialize rewrites both sides and materialises ts columns", {
  dt <- data.table::data.table(g = c(1, 1, 2, 2), t = c(1, 2, 1, 2),
                               x = c(1, 2, 3, 4), y = c(10, 20, 30, 40))
  pm <- panel_materialize(y ~ factor(g) + lag(x), dt, "g", "t")

  expect_equal(deparse(pm$formula), "y ~ factor(g) + .pt1")
  expect_true(".pt1" %in% names(pm$data))
  expect_true(is.na(pm$data[g == 1 & t == 1, .pt1]))
  expect_equal(pm$data[g == 1 & t == 2, .pt1], 1)
  expect_equal(pm$data[g == 2 & t == 2, .pt1], 3)
  # original data untouched
  expect_false(".pt1" %in% names(dt))
})
