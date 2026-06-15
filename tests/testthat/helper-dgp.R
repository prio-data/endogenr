# Known-truth data-generating processes for behavioural tests ----------------
#
# All produce balanced panels with contiguous integer time and a complete
# initial state, so they can feed setup_system() directly. Column names are
# configurable; defaults are `unit` / `time`.

# Panel AR(1)-X DGP: y_{u,t} = rho * y_{u,t-1} + beta * x_{u,t-1} + e_{u,t},
# with an exogenous AR(1) covariate x and Gaussian innovations. The workhorse
# for parameter recovery, error propagation, and true-vs-wrong-model tests.
sim_panel_ar1 <- function(units = 30L, n_time = 40L, rho = 0.6, beta = 0.4,
                          sd = 1, x_rho = 0.5, x_sd = 1, seed = NULL,
                          groupvar = "unit", timevar = "time") {
  if (!is.null(seed)) set.seed(seed)
  rows <- vector("list", units)
  for (u in seq_len(units)) {
    x <- numeric(n_time)
    x[1] <- stats::rnorm(1, sd = x_sd)
    for (tt in 2:n_time) {
      x[tt] <- x_rho * x[tt - 1] + stats::rnorm(1, sd = x_sd)
    }
    y <- numeric(n_time)
    y[1] <- stats::rnorm(1, sd = sd / sqrt(max(1 - rho^2, 0.05)))
    for (tt in 2:n_time) {
      y[tt] <- rho * y[tt - 1] + beta * x[tt - 1] + stats::rnorm(1, sd = sd)
    }
    rows[[u]] <- data.table::data.table(
      unit = u, time = seq_len(n_time), x = x, y = y)
  }
  dt <- data.table::rbindlist(rows)
  data.table::setnames(dt, c("unit", "time"), c(groupvar, timevar))
  dt[]
}

# Ragged-panel wrapper around sim_panel_ar1: trims per-unit spans so units may
# enter late or exit early. `enter`/`exit` are named vectors mapping unit ->
# first/last retained time step (names are coerced with as.integer()).
sim_panel_ragged <- function(units = 6L, n_time = 40L, enter = NULL,
                             exit = NULL, groupvar = "unit",
                             timevar = "time", ...) {
  dt <- sim_panel_ar1(units = units, n_time = n_time,
                      groupvar = groupvar, timevar = timevar, ...)
  u  <- dt[[groupvar]]
  tt <- dt[[timevar]]
  keep <- rep(TRUE, nrow(dt))
  for (nm in names(enter)) keep <- keep & !(u == as.integer(nm) & tt < enter[[nm]])
  for (nm in names(exit))  keep <- keep & !(u == as.integer(nm) & tt > exit[[nm]])
  dt <- dt[keep]
  data.table::setkeyv(dt, c(groupvar, timevar))
  dt[]
}

# Pooled cross-sectional DGP: y = b0 + b1 * x1 + b2 * x2 + e (homoscedastic).
# Independent rows, correct functional form -> a setting where the per-equation
# estimator's assumptions actually hold (used as the calibration baseline).
sim_panel_pooled <- function(units = 40L, n_time = 30L,
                             b0 = 1, b1 = 2, b2 = -1, sd = 1, seed = NULL,
                             groupvar = "unit", timevar = "time") {
  if (!is.null(seed)) set.seed(seed)
  grid <- data.table::CJ(unit = seq_len(units), time = seq_len(n_time))
  n <- nrow(grid)
  grid[, x1 := stats::rnorm(n)]
  grid[, x2 := stats::rnorm(n)]
  grid[, y := b0 + b1 * x1 + b2 * x2 + stats::rnorm(n, sd = sd)]
  data.table::setnames(grid, c("unit", "time"), c(groupvar, timevar))
  grid[]
}

# Heteroscedastic DGP: y = b0 + b1 * x + e, e ~ N(0, sigma(z)) with
# log sigma^2 = g0 + g1 * z. Variance rises with z when g1 > 0.
sim_panel_hetero <- function(units = 40L, n_time = 30L,
                            b0 = 0, b1 = 1, g0 = -1, g1 = 1, seed = NULL,
                            groupvar = "unit", timevar = "time") {
  if (!is.null(seed)) set.seed(seed)
  grid <- data.table::CJ(unit = seq_len(units), time = seq_len(n_time))
  n <- nrow(grid)
  grid[, x := stats::rnorm(n)]
  grid[, z := stats::rnorm(n)]
  sigma <- sqrt(exp(g0 + g1 * grid$z))
  grid[, y := b0 + b1 * x + stats::rnorm(n, sd = sigma)]
  data.table::setnames(grid, c("unit", "time"), c(groupvar, timevar))
  grid[]
}

# Common-shock DGP: y_{u,t} = rho * y_{u,t-1} + s_t + e_{u,t}, where s_t is a
# shock SHARED across all units at time t (cross-unit residual correlation).
# Used to document the SUR / common-shock gap (innovations drawn independently).
sim_panel_common_shock <- function(units = 25L, n_time = 30L, rho = 0.5,
                                   common_sd = 1, idio_sd = 0.3, seed = NULL,
                                   groupvar = "unit", timevar = "time") {
  if (!is.null(seed)) set.seed(seed)
  shock <- stats::rnorm(n_time, sd = common_sd)        # shared across units
  rows <- vector("list", units)
  for (u in seq_len(units)) {
    y <- numeric(n_time)
    y[1] <- stats::rnorm(1, sd = idio_sd)
    for (tt in 2:n_time) {
      y[tt] <- rho * y[tt - 1] + shock[tt] + stats::rnorm(1, sd = idio_sd)
    }
    rows[[u]] <- data.table::data.table(unit = u, time = seq_len(n_time), y = y)
  }
  dt <- data.table::rbindlist(rows)
  data.table::setnames(dt, c("unit", "time"), c(groupvar, timevar))
  dt[]
}
