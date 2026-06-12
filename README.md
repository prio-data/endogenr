
<!-- README.md is generated from README.Rmd. Please edit that file -->

# endogenr

<!-- badges: start -->

<!-- badges: end -->

The goal of `endogenr` is to make it easy to simulate dynamic systems
from regression models, mathematical equations, and exogenous inputs
(either based on a stochastic distribution, or given by some data). It
assumes a panel-data structure with two columns identifying the time and
the unit dimensions.

The simulator identifies the dependency graph of the models added to the
system and derives the order of calculation from that graph. Parallel
execution is opt-in via the `future` package.

## Installation

You can install the development version of `endogenr` from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("prio-data/endogenr")
# alternatively with
renv::install("prio-data/endogenr")
```

You can also clone the repository, open it as a project in RStudio, find
the “Build” tab, and press “Install”.

## Example

`setup_system()` accepts a plain `data.frame` or `data.table` — no
`tsibble` conversion is required. Formula RHS terms can use `lag()`,
`zoo::rollmean()`, or any function you define yourself (pass user
functions through the `globals` argument).

``` r
library(endogenr)
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
df <- endogenr::example_data

# Drop units with any NA in modelled outcomes over the training window 1970–2009.
required <- c("gdppc", "gdppc_grwt", "best", "v2x_polyarchy", "psecprop", "population")
train_window <- df[df$year >= 1970 & df$year <= 2009, ]
ok <- aggregate(train_window[, required],
                by = list(gwcode = train_window$gwcode),
                FUN = function(x) all(!is.na(x)))
keep <- ok$gwcode[apply(ok[, required], 1, all)]
df   <- df[df$gwcode %in% keep, ]
df$gdp <- df$gdppc * df$population  # derived outcome; pre-computed so the initial state is non-NA

c1 <- best ~ lag(best) + lag(log(gdppc)) + lag(log(population)) +
  lag(psecprop) + lag(v2x_polyarchy) + lag(gdppc_grwt) +
  lag(zoo::rollmean(best, k = 5, fill = NA, align = "right"))

model_system <- list(
  build_model("deterministic", formula = gdppc ~ I(abs(lag(gdppc) * (1 + gdppc_grwt)))),
  build_model("deterministic", formula = gdp ~ I(abs(gdppc * population))),
  build_model("parametric_distribution", formula = ~gdppc_grwt, distribution = "norm"),
  build_model("linear", formula = c1, boot = "resid"),
  build_model("univariate_fable",
              formula = v2x_polyarchy ~ error("A") + trend("N") + season("N"),
              method = "ets"),
  build_model("exogen", formula = ~psecprop),
  build_model("exogen", formula = ~population)
)
```

### Validation

`setup_system()` runs `validate_panel()` and `validate_system_closure()`
on your inputs and builds the dependency graph (no models are fit yet).
These check that time is contiguous and integer-valued within each
unit, the initial state at `test_start - 1` has no NAs in any modelled
outcome for units present there, and every variable referenced by a
formula is either modelled or supplied as a column. Panels may be
unbalanced: units may enter late or exit early. Units without a row at
`test_start - 1` are used for training only and are excluded from the
simulation. With `factor()` terms, prefer pre-converted factor columns
so window fits keep all levels. See `?validate_panel` and
`?validate_system_closure` if you want to run them ad-hoc on a
candidate panel.

``` r
sys <- setup_system(
  models      = model_system,
  data        = df,
  train_start = 1970,
  test_start  = 2010,
  horizon     = 12,
  groupvar    = "gwcode",
  timevar     = "year",
  inner_sims  = 2,
  min_window  = 10
)
```

### Fit the system

`fit_system()` estimates the models and **stores** the fitted objects,
so the coefficients that drive the simulation can be inspected with
`get_coefficients()` or plotted with `plot_coefficients()`. `nsim` lives
here: it is the number of coefficient draws. Because this example sets
`min_window`, the bootstrapped linear model is refit on a random
training window for each draw, so the stored coefficients vary across
draws.

``` r
future::plan(future::multisession, workers = 2)

set.seed(42)
fit <- fit_system(sys, nsim = 2)

# The coefficients actually used across draws
get_coefficients(fit)
#>     .draw .window_start .window_end outcome
#>     <int>         <num>       <num>  <char>
#>  1:     1          1984        2010    best
#>  2:     1          1984        2010    best
#>  3:     1          1984        2010    best
#>  4:     1          1984        2010    best
#>  5:     1          1984        2010    best
#>  6:     1          1984        2010    best
#>  7:     1          1984        2010    best
#>  8:     1          1984        2010    best
#>  9:     2          1970        1987    best
#> 10:     2          1970        1987    best
#> 11:     2          1970        1987    best
#> 12:     2          1970        1987    best
#> 13:     2          1970        1987    best
#> 14:     2          1970        1987    best
#> 15:     2          1970        1987    best
#> 16:     2          1970        1987    best
#>                                              term      estimate    std.error
#>                                            <char>         <num>        <num>
#>  1:                                   (Intercept)  3.426502e+02 4.876389e+02
#>  2:                                      lag_best  3.767427e-01 2.035664e-02
#>  3:                                 lag_log_gdppc -4.406651e+01 5.986201e+01
#>  4:                            lag_log_population  3.539008e+01 3.115072e+01
#>  5:                                  lag_psecprop  1.036039e+03 9.693822e+02
#>  6:                             lag_v2x_polyarchy -1.779253e+02 2.050505e+02
#>  7:                                lag_gdppc_grwt  2.463622e+03 7.655523e+02
#>  8: lag_zoo_rollmean_best_k_5_fill_na_align_right  3.489623e-01 2.186694e-02
#>  9:                                   (Intercept) -8.430086e+02 8.902369e+02
#> 10:                                      lag_best  6.619216e-01 2.339521e-02
#> 11:                                 lag_log_gdppc  1.302106e+02 1.078525e+02
#> 12:                            lag_log_population  1.280010e+02 5.984988e+01
#> 13:                                  lag_psecprop  6.686443e+02 2.639668e+03
#> 14:                             lag_v2x_polyarchy -1.283131e+03 3.964505e+02
#> 15:                                lag_gdppc_grwt -1.387291e+03 1.508337e+03
#> 16: lag_zoo_rollmean_best_k_5_fill_na_align_right  6.765428e-02 2.654208e-02
#>      statistic       p.value
#>          <num>         <num>
#>  1:  0.7026719  4.823053e-01
#>  2: 18.5071116  3.785206e-73
#>  3: -0.7361349  4.616961e-01
#>  4:  1.1360920  2.559930e-01
#>  5:  1.0687621  2.852479e-01
#>  6: -0.8677147  3.856078e-01
#>  7:  3.2180975  1.301786e-03
#>  8: 15.9584422  1.833862e-55
#>  9: -0.9469485  3.437913e-01
#> 10: 28.2930395 3.698411e-146
#> 11:  1.2073030  2.274731e-01
#> 12:  2.1387009  3.259326e-02
#> 13:  0.2533062  8.000603e-01
#> 14: -3.2365470  1.231689e-03
#> 15: -0.9197488  3.578265e-01
#> 16:  2.5489444  1.088676e-02
```

### Checking the fit against a plain regression

The fit path is designed to match a plain pooled regression exactly.
With `boot = NULL` on a spec and no `min_window`, `fit_system()` fits
that spec once on the full training window
`[train_start, test_start - 1]`, and its coefficients equal `lm()` on
the materialized training data. With `boot = "resid"` the bootstrap
draws centre on that fit and use the same estimation sample (complete
cases on the model's own columns). If a spec ever looks off, compare it
against the reference regression directly:

``` r
# endogenr's fit (one shared draw, full window)
sys <- setup_system(
  list(build_model("linear", formula = y ~ lag(y) + lag(x)),
       build_model("exogen", formula = ~x)),
  dt, train_start = 1965, test_start = 2010, horizon = 5,
  groupvar = "unit", timevar = "time", inner_sims = 1
)
fit <- fit_system(sys, nsim = 1)
get_coefficients(fit)

# the same regression by hand
dt[, `:=`(lag_y = shift(y), lag_x = shift(x)), by = unit]
coef(lm(y ~ lag_y + lag_x, dt[time < 2010]))
```

### Parallel execution and progress

The simulator no longer manages a `future` plan internally and no longer
refits anything: `simulate_system()` predicts using the stored draws and
derives `nsim` from them. Set a plan yourself before calling, and wrap
the call in `progressr::with_progress()` if you want a progress bar:

``` r
set.seed(42)
progressr::with_progress({
  res <- simulate_system(fit)
})

future::plan(future::sequential)
```

### Scoring and plotting

`simulate_system()` stamps `panel_unit` / `panel_time` attributes on the
result so `get_accuracy()` and `plotsim()` can infer the panel context
without any further configuration. Filter to the forecast window using a
plain `data.table` predicate.

``` r
res <- res[res$year >= sys$test_start, ]

acc <- get_accuracy(res, "gdppc_grwt", df)
acc |>
  dplyr::summarize(dplyr::across(crps:winkler, ~ mean(.x))) |>
  dplyr::arrange(crps) |>
  knitr::kable()
```

|      crps |       mae |   winkler |
|----------:|----------:|----------:|
| 0.0354269 | 0.0441132 | 0.1456342 |

``` r

plotsim(res, "gdppc", c(2, 20, 530), df)
#> Warning: Removed 4 rows containing missing values or values outside the scale range
#> (`geom_line()`).
```

<img src="man/figures/README-postprocess-1.png" alt="" width="100%" />

``` r
plotsim(res, "gdppc_grwt", c(2, 20, 530), df)
#> Warning: Removed 7 rows containing missing values or values outside the scale range
#> (`geom_line()`).
```

<img src="man/figures/README-postprocess-2.png" alt="" width="100%" />

``` r
plotsim(res, "v2x_polyarchy", c(2, 20, 530), df)
```

<img src="man/figures/README-postprocess-3.png" alt="" width="100%" />

## Spatial lag

`endogenr` supports spatial-lag variables in the simulation. The setup
is two steps: (1) compute the spatial lag for the historical data, and
(2) register a `spatial_lag` model in the system so the lag is
recomputed at each simulated time step.

### Step 1: prepare spatial weights and the historical lag

`st_weights_from_sf()` builds a neighbourhood list and weights from an
`sf` object. The default is queen contiguity. For other schemes, call
`sfdep` directly and supply `nb`, `wt`, and `unit_ids` yourself (in the
order they appear in the spatial object).

Computing the historical lag can be fiddly when the panel is unbalanced
or has missing observations; the example below filters to units present
in the neighbourhood structure before computing the lag.

``` r
library(endogenr)
library(dplyr)

df <- endogenr::example_data

# Load a map and filter to units present in the data
map <- poldat::cshp_gw_modifications(france_overseas = FALSE) |>
  dplyr::filter(end == as.Date("2019-12-31"))
map <- map |> dplyr::filter(gwcode %in% unique(df$gwcode))

# Build spatial weights (queen contiguity by default)
sf::sf_use_s2(FALSE)
neigh <- st_weights_from_sf(map, "gwcode", weights_args = list(allow_zero = TRUE))

# Compute the spatial lag of `best` for each year in the historical data
df <- df |>
  dplyr::filter(gwcode %in% neigh$unit_ids) |>
  dplyr::group_by(year) |>
  dplyr::mutate(
    sl_best = sfdep::st_lag(best, neigh$nb, neigh$wt, allow_zero = TRUE)
  ) |>
  dplyr::ungroup()
```

### Step 2: add a `spatial_lag` model to the system

A `spatial_lag` model is a cross-sectional transformation applied at
every simulated `t`. Reference the lag in other formulas as `lag(sl_y)`
— using the same-period value `sl_y` would create a circular dependency.

``` r
c1 <- best ~ lag(best) + lag(sl_best) + lag(log(gdppc)) +
  lag(log(population)) + lag(psecprop) + lag(v2x_polyarchy) + lag(gdppc_grwt) +
  lag(zoo::rollmean(best, k = 5, fill = NA, align = "right"))

model_system <- list(
  # spatial_lag recomputes sl_best from `best` at each simulated t
  build_model("spatial_lag", formula = sl_best ~ best,
              nb = neigh$nb, wt = neigh$wt, unit_ids = neigh$unit_ids,
              island_default = 0),
  build_model("deterministic", formula = gdppc ~ I(abs(lag(gdppc) * (1 + gdppc_grwt)))),
  build_model("deterministic", formula = gdp ~ I(abs(gdppc * population))),
  build_model("parametric_distribution", formula = ~gdppc_grwt, distribution = "norm"),
  build_model("linear", formula = c1, boot = "resid"),
  build_model("univariate_fable",
              formula = v2x_polyarchy ~ error("A") + trend("N") + season("N"),
              method = "ets"),
  build_model("exogen", formula = ~psecprop),
  build_model("exogen", formula = ~population)
)

sys <- setup_system(
  models      = model_system,
  data        = df,
  train_start = 1970,
  test_start  = 2010,
  horizon     = 12,
  groupvar    = "gwcode",
  timevar     = "year",
  inner_sims  = 2,
  min_window  = 10
)

future::plan(future::multisession, workers = 2)
set.seed(42)
fit <- fit_system(sys, nsim = 2)
set.seed(42)
progressr::with_progress({
  res <- simulate_system(fit)
})
future::plan(future::sequential)
```

## Long-horizon comparison

For each horizon `h`, the long-horizon API fits a single *direct*
regression of the h-step-ahead outcome on covariates observed at the
forecast origin (`test_start - 1`, the last observed period — the same
information the dynamic simulator conditions on). It is a reduced-form
benchmark for the dynamic simulator. See `?cv_long_horizon` for the
cross-validated entry point.

The outcome on the formula LHS must be wrapped in `lead_horizon()`; the
horizon `h` is supplied internally per horizon. The RHS is evaluated at
the origin, so write the covariates *unlagged* for the standard
benchmark (use `lag()` only when you deliberately want history older
than the origin).

``` r
formulas <- list(
  lh_linear = lead_horizon(gdppc_grwt) ~ gdppc_grwt + log(gdppc) + best
)

lh_setup <- setup_long_horizon(
  data       = df,
  formulas   = formulas,
  horizons   = 1:12,
  groupvar   = "gwcode",
  timevar    = "year",
  test_start = 2010
)

lh_forecasts <- forecast_long_horizon(
  lh_setup,
  data       = df,
  test_start = 2010,   # defaults to lh_setup$test_start
  nsim       = 100,
  inner_sims = 10
)

# Score both approaches on the same (unit, horizon) grid, then stack them.
lh_acc  <- get_lh_accuracy(lh_forecasts, df, lh_setup)
sim_acc <- get_accuracy(res, "gdppc_grwt", df,
                        test_start = 2010, by = c("gwcode", "horizon"))

compare_approaches(lh_acc, sim_acc) |>
  dplyr::group_by(approach, horizon) |>
  dplyr::summarize(crps = mean(crps), .groups = "drop") |>
  dplyr::arrange(horizon, approach)
```

For a transformed outcome (e.g. `lead_horizon(asinh(gdppc)) ~ ...`),
score on the modelled scale with
`get_lh_accuracy(..., scale = "model")`, or on the native scale with
`scale = "native", inverse = sinh`.
