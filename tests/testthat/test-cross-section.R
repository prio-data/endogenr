# cross_section models ---------------------------------------------------------
#
# A cross_section model evaluates `outcome ~ I(expr)` across all units at one
# simulated time step, separately per simulation draw.

cs_panel <- function() {
  dt <- data.table::as.data.table(expand.grid(gwcode = 1:5, year = 2000:2014))
  data.table::setkeyv(dt, c("gwcode", "year"))
  set.seed(1)
  dt[, g := rnorm(.N, 0.02, 0.05)]
  dt[, x := cumsum(g), by = gwcode]
  dt[, x_l1 := data.table::shift(x), by = gwcode]
  dt <- dt[year > 2000]
  dt[, fr := quantile(x_l1, 0.9, na.rm = TRUE), by = year]
  dt[, rk := rank(x_l1), by = year]
  dt
}

cs_simulate <- function(fr_formula) {
  future::plan(future::sequential)
  on.exit(future::plan(future::sequential), add = TRUE)
  models <- list(
    # Stochastic growth driver so draws diverge.
    build_model("linear", formula = g ~ 1, boot = "resid"),
    build_model("deterministic", formula = x ~ I(lag(x) + g)),
    build_model("deterministic", formula = x_l1 ~ I(lag(x))),
    build_model("cross_section", formula = fr_formula),
    build_model("cross_section", formula = rk ~ I(rank(x_l1)))
  )
  setup <- setup_system(models, data = cs_panel(), train_start = 2001,
                        test_start = 2010, horizon = 3, groupvar = "gwcode",
                        timevar = "year", inner_sims = 2)
  set.seed(2)
  data.table::as.data.table(simulate_system(fit_system(setup, nsim = 2)))
}

test_that("cross_section broadcasts a per-draw cross-sectional quantile", {
  res <- cs_simulate(fr ~ I(quantile(x_l1, 0.9, na.rm = TRUE)))

  chk <- res[, .(ok = isTRUE(all.equal(
    fr, rep(unname(quantile(x_l1, 0.9, na.rm = TRUE)), .N)))), by = .(.sim, year)]
  expect_true(all(chk$ok))
  # Computed within each draw, not pooled across draws.
  expect_gt(data.table::uniqueN(res[year == 2012, fr[1], by = .sim]$V1), 1)
})

test_that("cross_section keeps a per-unit-length result", {
  res <- cs_simulate(fr ~ I(quantile(x_l1, 0.9, na.rm = TRUE)))

  chk <- res[, .(ok = isTRUE(all.equal(rk, rank(x_l1)))), by = .(.sim, year)]
  expect_true(all(chk$ok))
})

test_that("cross_section errors on a result of the wrong length", {
  # future.apply warns "Canceling all iterations" when a worker errors.
  expect_error(suppressWarnings(cs_simulate(fr ~ I(head(x_l1, 2)))),
               "returned length 2; expected 1 or 5")
})

test_that("build_model rejects time-series functions and non-I() RHS in cross_section", {
  expect_error(build_model("cross_section", fr ~ I(quantile(lag(x), 0.9))),
               "time-series functions \\(lag\\)")
  expect_error(build_model("cross_section", fr ~ quantile(x_l1, 0.9)),
               "single I\\(\\) term")
})
